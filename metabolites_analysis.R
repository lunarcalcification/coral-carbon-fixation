# ========= 基础设置 =========
setwd("D:/metabolites")
options(scipen = 999, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ropls)
  library(ggplot2)
  library(cowplot)
  library(readxl)
  library(scales)
  if (!requireNamespace("sysfonts", quietly = TRUE)) install.packages("sysfonts")
  library(sysfonts)
  library(showtext)
  library(dplyr)
  library(stringr)
})
set.seed(2025)

# （可选）Levene
.have_car <- requireNamespace("car", quietly = TRUE)
if (!.have_car) message("提示：建议安装 car 包以执行 Levene 检验：install.packages('car')")

# ---- 全局常量 ----
SCALE_C      <- "standard"
PERM_I       <- 200
CV_I         <- 7
SD_THRESH    <- 1e-6
COR_CUTOFF   <- 0.95
MIN_VIP      <- 1
FDR_THR      <- 0.05
L2FC_THR     <- 1
FORCE_LOG10  <- TRUE
P_THR        <- 0.05
RUN_LABEL    <- "pval0.05"
VOLCANO_COLS <- c("Down"="#9AB6FF","NS"="#BFBFBF","Up"="#FF9AA2")
AXIS_TEXT_SIZE   <- 16
AXIS_TITLE_SIZE  <- 18
LEGEND_TEXT_SIZE <- 14
LEGEND_TITLE_SIZE<- 14
TITLE_SIZE       <- 18
TICK_LEN_PT      <- 4

# —— 自定义分组配色（仅对 AH_* 三组指定）——
GROUP_COLORS <- c(
  "AH_27" = "#E64B35", "PD_27" = "#E64B35", "PL_27" = "#E64B35",
  "AH_31" = "#4DBBD5", "PD_31" = "#4DBBD5", "PL_31" = "#4DBBD5",
  "AH_34" = "#00A087", "PD_34" = "#00A087", "PL_34" = "#00A087"
)

# —— 变量筛选开关（为避免未定义报错，这里显式设为 FALSE）——
FILTER_VARS <- FALSE

# ========= 9个配对（固定顺序）=========
desired_pairs <- list(
  c("AH_27","AH_31"), c("AH_27","AH_34"), c("AH_31","AH_34"),
  c("PD_27","PD_31"), c("PD_27","PD_34"), c("PD_31","PD_34"),
  c("PL_27","PL_31"), c("PL_27","PL_34"), c("PL_31","PL_34")
)
pair_order <- vapply(desired_pairs, function(x) sprintf("%s_vs_%s", x[1], x[2]), character(1))

# ========= 工具函数 =========
.safe_log10 <- function(M){
  M <- as.matrix(M); mode(M) <- "numeric"
  pos <- M[is.finite(M) & M > 0]
  pseudo <- if (length(pos)) max(1e-9, 0.5*stats::quantile(pos,0.01,na.rm=TRUE)) else 1e-9
  M <- pmax(M, pseudo)
  base::log10(M)
}

.safe_pick_num <- function(df, rownm, colnm){
  if (is.null(df) || !is.data.frame(df)) return(NA_real_)
  if (!(rownm %in% rownames(df)) || !(colnm %in% colnames(df))) return(NA_real_)
  v <- suppressWarnings(as.numeric(df[rownm, colnm])); if (length(v)==0) NA_real_ else v
}

.extract_ropls_metrics <- function(mdl){
  met <- list(strategy=NA_character_, R2Y=NA_real_, Q2=NA_real_,
              pR2Y=NA_real_, pQ2=NA_real_, pre=NA_integer_, ort=NA_integer_, ok=FALSE)
  if (is.null(mdl)) return(met)
  type_guess <- tryCatch(mdl@typeC, error=function(e) NA_character_)
  sumdf <- tryCatch(mdl@summaryDF, error=function(e) NULL)
  if (!is.null(sumdf) && nrow(sumdf)>0){
    rn <- rownames(sumdf); rp <- if ("Total" %in% rn) "Total" else rn[1]
    met$R2Y  <- suppressWarnings(as.numeric(sumdf[rp,"R2Y(cum)"]))
    met$Q2   <- suppressWarnings(as.numeric(sumdf[rp,"Q2(cum)"]))
    met$pR2Y <- suppressWarnings(as.numeric(sumdf[rp,"pR2Y"]))
    met$pQ2  <- suppressWarnings(as.numeric(sumdf[rp,"pQ2"]))
    if ("pre" %in% colnames(sumdf)) met$pre <- suppressWarnings(as.integer(sumdf[rp,"pre"]))
    if ("ort" %in% colnames(sumdf)) met$ort <- suppressWarnings(as.integer(sumdf[rp,"ort"]))
  } else {
    s <- tryCatch(summary(mdl), silent=TRUE)
    if (!inherits(s,"try-error") && is.data.frame(s) && nrow(s)>0){
      rn <- rownames(s); rp <- if ("Total" %in% rn) "Total" else rn[1]
      g <- function(cn) if (cn %in% colnames(s)) suppressWarnings(as.numeric(s[rp,cn])) else NA_real_
      met$R2Y <- g("R2Y(cum)"); met$Q2 <- g("Q2(cum)")
      met$pR2Y <- g("pR2Y"); met$pQ2 <- g("pQ2")
      if ("pre" %in% colnames(s)) met$pre <- suppressWarnings(as.integer(s$pre[1]))
      if ("ort" %in% colnames(s)) met$ort <- suppressWarnings(as.integer(s$ort[1]))
    }
  }
  met$ok <- any(is.finite(c(met$R2Y,met$Q2,met$pR2Y,met$pQ2)))
  if (met$ok && is.character(type_guess) && !is.na(type_guess)){
    if (grepl("opls",type_guess,TRUE) && grepl("da",type_guess,TRUE)) met$strategy <- "OPLS-DA"
    else if (grepl("pls",type_guess,TRUE) && grepl("da",type_guess,TRUE)) met$strategy <- "PLS-DA"
  }
  met
}

any_finite <- function(...) isTRUE(any(is.finite(c(...))))

get_pair <- function(X, info, g1, g2){
  sel <- info$Group %in% c(g1,g2)
  Xp  <- X[sel,,drop=FALSE]
  grp <- droplevels(factor(info$Group[sel], levels=c(g1,g2)))
  list(X=Xp, group=grp, info=info[sel,,drop=FALSE])
}

# （保留但默认不用）
filter_vars <- function(Xp, sd_thresh=SD_THRESH, cor_cutoff=COR_CUTOFF){
  keep <- apply(Xp, 2, function(v) sd(v,na.rm=TRUE) > sd_thresh)
  Xf <- Xp[, keep, drop=FALSE]
  if (ncol(Xf) >= 2){
    cm <- suppressWarnings(cor(Xf, use="pairwise.complete.obs")); diag(cm) <- 0
    drop <- rep(FALSE, ncol(Xf))
    for (j in seq_len(ncol(Xf))){
      if (drop[j]) next
      hi <- which(abs(cm[j,]) > cor_cutoff & !drop); drop[hi] <- TRUE; drop[j] <- FALSE
    }
    Xf <- Xf[, !drop, drop=FALSE]
  }
  Xf
}

# ========== 统一美化画图（仅保存 PDF，标题含 R2Y/Q2） ==========
.get_slot <- function(obj, slot) tryCatch(slot(obj, slot), error=function(e) NULL)
.cols2 <- c("#F59E0B","#60A5FA")
.pt    <- 2.4

.plotnum <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "NA")

plot_scores_pretty <- function(mdl, group, tag, outprefix, r2y=NA_real_, q2=NA_real_){
  tMN <- .get_slot(mdl, "scoreMN"); if (is.null(tMN)) tMN <- .get_slot(mdl, "xScoreMN")
  oMN <- .get_slot(mdl, "orthoScoreMN")
  if (is.null(tMN)) return(invisible(NULL))
  t1 <- as.numeric(tMN[,1]); t1[!is.finite(t1)] <- 0
  o1 <- if (!is.null(oMN)) as.numeric(oMN[,1]) else rep(0, length(t1)); o1[!is.finite(o1)] <- 0
  
  df <- data.frame(t1=t1, o1=o1, grp=factor(group))
  xlab <- "T score[1]"
  ylab <- if (!isTRUE(all(o1==0, na.rm=TRUE))) "Orthogonal T score[1]" else "T(ortho) score[1]"
  ttl  <- sprintf("%s | R2Y=%s, Q2=%s", tag, .plotnum(r2y), .plotnum(q2))
  
  # 配色：AH_* 用指定色，其它组自动分配互异色
  lvl <- levels(df$grp)
  vals <- GROUP_COLORS[intersect(names(GROUP_COLORS), lvl)]
  missing <- setdiff(lvl, names(vals))
  if (length(missing) > 0) {
    extra <- scales::hue_pal()(length(missing))
    names(extra) <- missing
    vals <- c(vals, extra)
  }
  vals <- vals[lvl]
  
  p <- ggplot(df, aes(t1, o1, color=grp)) +
    geom_hline(yintercept=0, linetype=2, linewidth=0.4, color="grey40") +
    geom_vline(xintercept=0, linetype=2, linewidth=0.4, color="grey40") +
    stat_ellipse(type="norm", linewidth=0.9, alpha=0.25, show.legend=FALSE) +
    geom_point(size=.pt) +
    scale_color_manual(values = vals, name = NULL, drop = FALSE) +
    labs(title = ttl, x = xlab, y = ylab) +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(face="bold", hjust=0.5), legend.position="right")
  
  ggsave(paste0(outprefix, "_scores.pdf"), p, width=7, height=5.2)
  invisible(p)
}

# ========= 智能建模（不再进行变量筛选；仅输出得分图 PDF）=========
run_pair_smart <- function(pr, g1, g2,
                           scaleC=SCALE_C, permI=PERM_I, cvI=CV_I,
                           outdir="OPLS_pairwise_SMART", bad_frac_thresh=Inf,
                           cor_cutoff=COR_CUTOFF){
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
  comp_tag <- sprintf("%s_vs_%s", g1, g2)
  
  if (isTRUE(FILTER_VARS)) Xf <- filter_vars(pr$X, sd_thresh=SD_THRESH, cor_cutoff=cor_cutoff) else Xf <- pr$X
  kept_vars <- colnames(Xf)
  
  if (length(kept_vars) < 2){
    summary_row <- data.frame(
      Comparison=comp_tag, Strategy=NA, Transform="log10", Scale=scaleC,
      Vars_Kept=length(kept_vars), R2Y=NA, Q2=NA, pR2Y=NA, pQ2=NA, Samples=nrow(pr$X)
    )
    return(list(summary=summary_row, kept=kept_vars, mdl=NULL, vip=NULL, trans="log10"))
  }
  
  mdl <- tryCatch(ropls::opls(Xf, pr$group, predI=1, orthoI=NA, scaleC=scaleC, permI=permI, crossvalI=cvI), error=function(e) NULL)
  met <- .extract_ropls_metrics(mdl)
  if (is.null(mdl) || !met$ok){
    mdl <- tryCatch(ropls::opls(Xf, pr$group, predI=1, orthoI=0, scaleC=scaleC, permI=permI, crossvalI=cvI), error=function(e) NULL)
    met <- .extract_ropls_metrics(mdl)
  }
  vip_vec <- tryCatch(ropls::getVipVn(mdl), error=function(e) NULL)
  
  if (!is.null(mdl) && met$ok && !is.na(met$strategy)){
    outprefix <- file.path(outdir, sprintf("%s_%s_%s_%s", met$strategy, comp_tag, "log10", scaleC))
    plot_scores_pretty(mdl, pr$group, comp_tag, outprefix, r2y = met$R2Y, q2 = met$Q2)
  }
  
  summary_row <- data.frame(
    Comparison=comp_tag,
    Strategy=if (met$ok) met$strategy else NA_character_,
    Transform="log10", Scale=scaleC, Vars_Kept=length(kept_vars),
    R2Y=met$R2Y, Q2=met$Q2, pR2Y=met$pR2Y, pQ2=met$pQ2, Samples=nrow(pr$X)
  )
  list(summary=summary_row, kept=kept_vars, mdl=if (met$ok) mdl else NULL,
       vip=if (met$ok) vip_vec else NULL, trans="log10")
}

# ========= 单次（混合极性）完整流程 =========
run_pipeline_mixed <- function(infile){
  message(sprintf("===== 开始处理：%s（MIX） =====", infile))
  stopifnot(file.exists(infile))
  set.seed(123)
  
  dat <- readxl::read_excel(infile, sheet = 1)
  dat <- as.data.frame(dat, check.names = FALSE)
  names(dat) <- trimws(names(dat))
  
  # ---- 智能识别 PubChem 列并统一名为 PubChemID ----
  pubchem_candidates <- c("PubChemID", "PubChem", "PubChem_ID", "PubChem CID")
  PUBCHEM_COL <- intersect(pubchem_candidates, names(dat))
  PUBCHEM_COL <- if (length(PUBCHEM_COL)) PUBCHEM_COL[1] else character(0)
  
  # 样本/注释列
  sample_cols <- grep("^(AH|PD|PL)(27|31|34)-\\d+$", names(dat), value=TRUE)
  subdata <- dat[, c("Compounds", sample_cols), drop=FALSE]
  
  # 基础注释列：Compounds/HMDB/C_ID +（若存在）PubChemID
  anno_basic_cols <- c("Compounds","HMDB","C_ID","Pathway")
  if (length(PUBCHEM_COL)) anno_basic_cols <- c(anno_basic_cols, PUBCHEM_COL)
  anno_basic_cols <- intersect(anno_basic_cols, names(dat))
  anno_basic <- dat[, anno_basic_cols, drop=FALSE]
  if (length(PUBCHEM_COL) && PUBCHEM_COL != "PubChemID") {
    names(anno_basic)[names(anno_basic) == PUBCHEM_COL] <- "PubChemID"
  }
  
  samples <- setdiff(names(subdata), "Compounds")
  species <- sub("^(AH|PD|PL).*", "\\1", samples)
  temp    <- sub("^(?:AH|PD|PL)(\\d{2}).*", "\\1", samples)
  rep_id  <- sub(".*-(\\d+)$", "\\1", samples)
  sample_info <- data.frame(
    Sample=samples, Species=species, Temp=temp,
    Group=paste(species, temp, sep="_"), Replicate=rep_id
  )
  
  missing_groups <- setdiff(unique(unlist(desired_pairs)), unique(sample_info$Group))
  if (length(missing_groups)>0) stop(sprintf("以下分组在数据中不存在：%s", paste(missing_groups, collapse=", ")))
  
  stopifnot(!anyDuplicated(subdata$Compounds))
  X0 <- as.matrix(subdata[, samples, drop=FALSE]); mode(X0) <- "numeric"
  rownames(X0) <- subdata$Compounds; colnames(X0) <- samples
  X_raw <- t(X0)
  ord <- match(sample_info$Sample, rownames(X_raw)); stopifnot(!any(is.na(ord)))
  X_raw <- X_raw[ord,,drop=FALSE]
  stopifnot(identical(rownames(X_raw), sample_info$Sample))
  cat(sprintf("[MIX] X_raw 维度：%d × %d\n", nrow(X_raw), ncol(X_raw)))
  
  X_log <- if (FORCE_LOG10) .safe_log10(X_raw) else X_raw
  cat(sprintf("[MIX] X_log 维度：%d × %d（%s）\n",
              nrow(X_log), ncol(X_log), if (FORCE_LOG10) "log10" else "raw"))
  
  outdir_main <- sprintf("OPLS_pairwise_SMART_MIX_%s", RUN_LABEL)
  if (!dir.exists(outdir_main)) dir.create(outdir_main, recursive=TRUE)
  diff_outdir <- file.path(outdir_main, "DEMs"); if (!dir.exists(diff_outdir)) dir.create(diff_outdir, TRUE)
  
  results_list <- vector("list", length(desired_pairs))
  pairwise_summary <- vector("list", length(desired_pairs))
  kept_union <- character(0)
  
  for (i in seq_along(desired_pairs)){
    g1 <- desired_pairs[[i]][1]; g2 <- desired_pairs[[i]][2]
    message(sprintf("===> SMART: %s vs %s (MIX)", g1, g2))
    pr_log <- get_pair(X_log, sample_info, g1, g2)
    pr_raw <- get_pair(X_raw, sample_info, g1, g2)
    res <- run_pair_smart(pr_log, g1, g2, scaleC=SCALE_C, permI=PERM_I, cvI=CV_I,
                          outdir=outdir_main, bad_frac_thresh=Inf, cor_cutoff=COR_CUTOFF)
    results_list[[i]] <- list(pair=c(g1,g2), pr_log=pr_log, pr_raw=pr_raw, res=res)
    pairwise_summary[[i]] <- res$summary
    kept_union <- union(kept_union, res$kept)
  }
  
  pairwise_summary <- do.call(rbind, pairwise_summary)
  
  get_comp_counts <- function(mdl){
    pre <- NA_integer_; ort <- NA_integer_
    if (!is.null(mdl)){
      s <- try(summary(mdl), silent=TRUE)
      if (!inherits(s,"try-error") && is.data.frame(s) && all(c("pre","ort") %in% colnames(s))){
        pre <- suppressWarnings(as.integer(s$pre[1])); ort <- suppressWarnings(as.integer(s$ort[1]))
      } else if (!is.null(mdl@summaryDF) && all(c("pre","ort") %in% colnames(mdl@summaryDF))){
        pre <- suppressWarnings(as.integer(mdl@summaryDF$pre[1])); ort <- suppressWarnings(as.integer(mdl@summaryDF$ort[1]))
      }
    }
    c(pre=pre, ort=ort)
  }
  pred_list <- orth_list <- integer(length(results_list))
  for (i in seq_along(results_list)){
    cc <- get_comp_counts(results_list[[i]]$res$mdl)
    pred_list[i] <- cc[["pre"]]; orth_list[i] <- cc[["ort"]]
  }
  pairwise_summary$PredComp <- pred_list
  pairwise_summary$OrthComp <- orth_list
  pairwise_summary$Transform <- "log10"
  
  out_sum <- file.path(outdir_main, "OPLSDA_PLSDA_SMART_pairwise_summary_selected9.csv")
  write.csv(pairwise_summary, out_sum, row.names=FALSE)
  print(pairwise_summary)
  
  # 背景并集（无筛选时等于全变量并集）
  bg_df <- data.frame(Compounds=sort(unique(kept_union)))
  write.csv(bg_df, file.path(diff_outdir, "background_metabolites_union.csv"), row.names=FALSE)
  
  make_diff_table_for_pair <- function(pr_raw, g1, g2, res_obj,
                                       min_vip=MIN_VIP, p_thr=P_THR, l2fc_thr=L2FC_THR){
    if (is.null(res_obj$mdl) || is.null(res_obj$vip) || length(res_obj$vip)==0) return(NULL)
    vars <- names(res_obj$vip)
    idx1 <- pr_raw$group==g1; idx2 <- pr_raw$group==g2
    Xraw <- pr_raw$X[, vars, drop=FALSE]
    mean1 <- colMeans(Xraw[idx1,,drop=FALSE], na.rm=TRUE)
    mean2 <- colMeans(Xraw[idx2,,drop=FALSE], na.rm=TRUE)
    log2FC <- log2((mean2+1e-9)/(mean1+1e-9))
    pvals <- sapply(vars, function(v){
      tryCatch(stats::t.test(Xraw[idx2, v], Xraw[idx1, v], var.equal=FALSE)$p.value, error=function(e) NA_real_)
    })
    fdr <- p.adjust(pvals, method="BH")
    
    df <- data.frame(
      Compounds=vars, VIP=as.numeric(res_obj$vip[vars]),
      mean_g1=as.numeric(mean1[vars]), mean_g2=as.numeric(mean2[vars]),
      log2FC=as.numeric(log2FC[vars]), pvalue=as.numeric(pvals[vars]), FDR=as.numeric(fdr[vars]),
      Strategy=res_obj$summary$Strategy[1], Transform="log10", stringsAsFactors=FALSE
    )
    df$Sig <- ifelse(!is.na(df$pvalue)&!is.na(df$VIP)&!is.na(df$log2FC) &
                       (df$VIP>=min_vip & df$pvalue < P_THR & abs(df$log2FC)>l2fc_thr),
                     ifelse(df$log2FC>0, "Up", "Down"), "NS")
    df$Direction <- ifelse(df$Sig=="Up","Up", ifelse(df$Sig=="Down","Down", NA))
    
    df_sig <- df[df$Sig!="NS", , drop=FALSE]
    df_sig <- df_sig[order(-df_sig$VIP, -abs(df_sig$log2FC)), , drop=FALSE]
    
    tag <- sprintf("%s_vs_%s", g1, g2)
    # 合并基础注释（含 PubChemID）
    df_out     <- if (exists("anno_basic")) merge(df,     anno_basic, by="Compounds", all.x=TRUE) else df
    df_sig_out <- if (exists("anno_basic")) merge(df_sig, anno_basic, by="Compounds", all.x=TRUE) else df_sig
    
    write.csv(df_out,     file.path(diff_outdir, sprintf("DEMs_stats_%s_all.csv", tag)), row.names=FALSE)
    write.csv(df_sig_out, file.path(diff_outdir, sprintf("DEMs_%s.csv",       tag)), row.names=FALSE)
    
    # 额外导出上/下调 C_ID 供后续分析
    if ("C_ID" %in% names(df_sig_out)){
      up_ids   <- na.omit(unique(df_sig_out$C_ID[df_sig_out$Direction=="Up"]))
      down_ids <- na.omit(unique(df_sig_out$C_ID[df_sig_out$Direction=="Down"]))
      write.table(up_ids,   file.path(diff_outdir, paste0("C_ID_up_",   tag, ".txt")), row.names=FALSE, col.names=FALSE, quote=FALSE)
      write.table(down_ids, file.path(diff_outdir, paste0("C_ID_down_", tag, ".txt")), row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
    if (nrow(df_sig_out)>0){
      df_sig_out$Pair <- tag
      df_sig_out <- df_sig_out[, c("Pair", setdiff(names(df_sig_out), "Pair"))]
    }
    df_sig_out
  }
  
  dems_all <- NULL
  for (tag in pair_order){
    sp <- strsplit(tag, "_vs_")[[1]]; g1 <- sp[1]; g2 <- sp[2]
    i <- which(vapply(results_list,function(z) paste(z$pair, collapse="_vs_"), character(1))==tag)
    if (length(i)==1){
      pr_raw_i <- results_list[[i]]$pr_raw; res_obj <- results_list[[i]]$res
      d <- make_diff_table_for_pair(pr_raw_i, g1, g2, res_obj)
      if (!is.null(d) && nrow(d)>0) dems_all <- if (is.null(dems_all)) d else rbind(dems_all, d)
    }
  }
  if (!is.null(dems_all) && nrow(dems_all)>0){
    write.csv(dems_all, file.path(diff_outdir, "DEMs_all_pairs.csv"), row.names=FALSE)
  }
  
  # 差异代谢物计数
  counts_list <- list()
  rs_index <- setNames(seq_along(results_list),
                       vapply(results_list,function(z) paste(z$pair, collapse="_vs_"), character(1)))
  for (tag in pair_order){
    f_sig <- file.path(diff_outdir, sprintf("DEMs_%s.csv", tag))
    if (!file.exists(f_sig)) next
    dem <- try(read.csv(f_sig), silent=TRUE); if (inherits(dem,"try-error") || nrow(dem)==0) next
    up_n <- sum(dem$Direction=="Up", na.rm=TRUE)
    down_n <- sum(dem$Direction=="Down", na.rm=TRUE)
    tot_n <- nrow(dem)
    i <- rs_index[[tag]]
    strat <- results_list[[i]]$res$summary$Strategy[1]
    r2y <- suppressWarnings(as.numeric(results_list[[i]]$res$summary$R2Y[1]))
    q2  <- suppressWarnings(as.numeric(results_list[[i]]$res$summary$Q2[1]))
    prc <- pairwise_summary$PredComp[i]; orc <- pairwise_summary$OrthComp[i]
    counts_list[[tag]] <- data.frame(
      Pair=tag, Up=up_n, Down=down_n, Total=tot_n,
      Strategy=strat, R2Y=r2y, Q2=q2, PredComp=prc, OrthComp=orc, stringsAsFactors=FALSE
    )
  }
  dem_counts <- do.call(rbind, counts_list[pair_order[pair_order %in% names(counts_list)]])
  if (!is.null(dem_counts) && nrow(dem_counts)>0){
    write.csv(dem_counts, file.path(diff_outdir,"DEMs_counts_by_pair.csv"), row.names=FALSE)
  } else message("[计数表] 没有可统计的配对。")
  
  # 火山图（仅保存 PDF）
  for (tag in pair_order) {
    f_all <- file.path(diff_outdir, sprintf("DEMs_stats_%s_all.csv", tag))
    if (!file.exists(f_all)) next
    dem <- read.csv(f_all); if (nrow(dem) == 0) next
    
    dem$Sig <- factor(dem$Sig, levels = c("Down", "NS", "Up"))
    up_n   <- sum(dem$Sig == "Up", na.rm = TRUE)
    down_n <- sum(dem$Sig == "Down", na.rm = TRUE)
    ns_n   <- sum(dem$Sig == "NS", na.rm = TRUE)
    
    p <- ggplot(dem, aes(x = log2FC, y = -log10(pvalue), color = Sig)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_vline(xintercept = c(-L2FC_THR, L2FC_THR), linetype = 2, color = "black") +
      geom_hline(yintercept = -log10(P_THR), linetype = 2, color = "black") +
      scale_color_manual(
        values = VOLCANO_COLS,
        labels = c(
          paste0("Down (n = ", down_n, ")"),
          paste0("Not Significant (n = ", ns_n, ")"),
          paste0("Up (n = ", up_n, ")")
        ),
        name = NULL
      ) +
      labs(x = expression(log[2]~Fold~Change),
           y = expression(-log[10](P~value)),
           title = paste("Volcano Plot of", tag),
           color = NULL) +
      theme_classic(base_size = 16) +
      theme(
        axis.line         = element_line(color = "black", linewidth = 0.6),
        panel.border      = element_blank(),
        axis.ticks        = element_line(color = "black"),
        axis.ticks.length = grid::unit(TICK_LEN_PT, "pt"),
        panel.grid        = element_blank(),
        legend.position   = c(0.8, 0.7),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        axis.text         = element_text(size = AXIS_TEXT_SIZE,  color = "black"),
        axis.title        = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
        legend.text       = element_text(size = LEGEND_TEXT_SIZE),
        legend.title      = element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
        plot.title        = element_text(size = TITLE_SIZE, hjust = 0.5, face = "bold")
      ) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    
    ggsave(file.path(diff_outdir, paste0("Volcano_", tag, ".pdf")), p, width = 6, height = 5)
  }
  
  # 富集准备：把 PubChemID 一并带上
  anno_cols <- c("Compounds","HMDB","KEGG_A_class","KEGG_B_class","Pathway","C_ID")
  if ("PubChemID" %in% names(anno_basic)) anno_cols <- c(anno_cols, "PubChemID")
  anno_cols <- intersect(anno_cols, names(dat))
  anno <- dat[, anno_cols, drop=FALSE]
  
  bg <- read.csv(file.path(diff_outdir, "background_metabolites_union.csv"))
  id_cols_for_bg <- intersect(c("Compounds","HMDB","C_ID","PubChemID"), names(anno))
  bg_id <- merge(bg, anno[, id_cols_for_bg, drop=FALSE], by="Compounds", all.x=TRUE)
  write.csv(bg_id, file.path(diff_outdir, "background_with_IDs.csv"), row.names=FALSE)
  
  dem_with_ids_list <- list()
  for (tag in pair_order){
    f_sig <- file.path(diff_outdir, sprintf("DEMs_%s.csv", tag))
    if (!file.exists(f_sig)) next
    dem <- read.csv(f_sig); if (nrow(dem)==0) next
    dem2 <- merge(dem, anno, by="Compounds", all.x=TRUE)
    dem2$Pair <- tag
    dem2 <- dem2[, c("Pair", setdiff(names(dem2),"Pair"))]
    write.csv(dem2, file.path(diff_outdir, paste0("DEMs_with_IDs_", tag, ".csv")), row.names=FALSE)
    
    if ("HMDB" %in% names(dem2)){
      up   <- subset(dem2, Direction=="Up")
      down <- subset(dem2, Direction=="Down")
      write.table(na.omit(unique(up$HMDB)),   file.path(diff_outdir, paste0("HMDB_up_",   tag, ".txt")),
                  row.names=FALSE, col.names=FALSE, quote=FALSE)
      write.table(na.omit(unique(down$HMDB)), file.path(diff_outdir, paste0("HMDB_down_", tag, ".txt")),
                  row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
    dem_with_ids_list[[tag]] <- dem2
  }
  dems_with_ids_all <- do.call(rbind, dem_with_ids_list[pair_order[pair_order %in% names(dem_with_ids_list)]])
  if (!is.null(dems_with_ids_all) && nrow(dems_with_ids_all)>0){
    write.csv(dems_with_ids_all, file.path(diff_outdir, "DEMs_with_IDs_all_pairs.csv"), row.names=FALSE)
  }
  
  
  
  message(sprintf("===== 完成：%s（MIX） =====", infile))
}

# ─────────────────────────────────────────
# PCA：字体注册 + 缩放函数 + PCA主函数（支持 xlsx/csv）
# ─────────────────────────────────────────

# 加载 showtext 包并注册 Times New Roman 字体
# ─────────────────────────────────────────
# 字体注册（一次性放在这里）
# ─────────────────────────────────────────
suppressPackageStartupMessages({
  if (!requireNamespace("sysfonts", quietly = TRUE)) install.packages("sysfonts")
  if (!requireNamespace("showtext",  quietly = TRUE)) install.packages("showtext")
  library(sysfonts)
  library(showtext)
})

# 在 Windows 系统字体目录里找 Times New Roman
find_font <- function(patterns) {
  cand <- unlist(lapply(patterns, function(p)
    list.files("C:/Windows/Fonts", pattern = p, ignore.case = TRUE, full.names = TRUE)))
  cand[which(file.exists(cand))][1]
}

ttf_regular    <- find_font(c("^times\\.ttf$", "^times new roman.*regular.*\\.ttf$", "^times new roman\\.ttf$"))
ttf_bold       <- find_font(c("^timesbd\\.ttf$", "^times new roman.*bold.*\\.ttf$"))
ttf_italic     <- find_font(c("^timesi\\.ttf$", "^times new roman.*italic.*\\.ttf$"))
ttf_bolditalic <- find_font(c("^timesbi\\.ttf$", "^times new roman.*bold italic.*\\.ttf$", "^times new roman.*italic bold.*\\.ttf$"))

if (!is.na(ttf_regular) && nzchar(ttf_regular)) {
  # 关键：用 sysfonts::font_add（不是 showtext::font_add）
  sysfonts::font_add(
    family = "Times New Roman",
    regular = ttf_regular,
    bold = if (!is.na(ttf_bold)) ttf_bold else ttf_regular,
    italic = if (!is.na(ttf_italic)) ttf_italic else ttf_regular,
    bolditalic = if (!is.na(ttf_bolditalic)) ttf_bolditalic else ttf_regular
  )
  showtext::showtext_auto()  # 启用 showtext 字形渲染
} else {
  warning("⚠️ 未能找到 Times New Roman 字体")
}



# 缩放方法（与上游常量兼容）
pareto_scale <- function(X) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  sds <- apply(Xc, 2, sd, na.rm = TRUE)
  sds[sds == 0 | is.na(sds)] <- 1
  Xc / sqrt(sds)
}
apply_scaling <- function(X, method = c("autoscale","pareto","none")) {
  method <- match.arg(method)
  if (method == "autoscale") {
    scale(X, center = TRUE, scale = TRUE)
  } else if (method == "pareto") {
    pareto_scale(X)
  } else {
    X
  }
}

# 读取表格：自动识别 xlsx/csv
.read_table <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path, sheet = sheet, .name_repair = "minimal") |> as.data.frame()
  } else if (ext %in% c("csv", "txt")) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else stop("Unsupported file type: ", ext)
}

# PCA 主函数
run_pca_from_file <- function(
    path,
    sheet = NULL,                      # Excel 可选
    sample_regex = "^(AH|PD|PL|QC)",
    valid_pat    = "^(AH|PD|PL)(27|31|34)-\\d+$|^QC\\d+$",
    use_log2     = TRUE,
    pseudocount  = 0,
    scaling      = c("autoscale","pareto","none"),
    title        = "PCA of Metabolomics"
){
  scaling <- match.arg(scaling)
  
  df0 <- .read_table(path, sheet = sheet)
  stopifnot(ncol(df0) >= 2)
  
  # 第一列为代谢物名称则设为行名
  first_col <- df0[[1]]
  num_rate <- mean(suppressWarnings(!is.na(as.numeric(first_col))))
  if (num_rate < 0.5) {
    rownames(df0) <- as.character(first_col)
    df <- df0[, -1, drop = FALSE]
  } else {
    df <- df0
  }
  
  cn <- trimws(colnames(df)); colnames(df) <- cn
  sample_cols <- cn[stringr::str_detect(cn, sample_regex)]
  if (length(sample_cols) == 0) stop("No sample columns found matching regex: ", sample_regex)
  
  valid_cols   <- sample_cols[grepl(valid_pat, sample_cols)]
  invalid_cols <- setdiff(sample_cols, valid_cols)
  if (length(invalid_cols)) {
    message("Dropped non-standard sample columns:"); print(invalid_cols)
  }
  sample_cols <- valid_cols
  if (length(sample_cols) == 0) stop("No valid sample columns under strict pattern.")
  
  mat <- df |>
    dplyr::select(dplyr::all_of(sample_cols)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ suppressWarnings(as.numeric(.)))) |>
    as.matrix()
  
  if (use_log2) mat <- log2(mat + pseudocount)
  
  X  <- t(mat)
  Xs <- apply_scaling(X, method = scaling)
  
  pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1_lab <- paste0("PC1 (", sprintf("%.1f", 100 * var_exp[1]), "%)")
  pc2_lab <- paste0("PC2 (", sprintf("%.1f", 100 * var_exp[2]), "%)")
  
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$Sample <- rownames(scores)
  scores$Coral <- stringr::str_extract(scores$Sample, "^(AH|PD|PL|QC)")
  scores$Temp  <- ifelse(scores$Coral == "QC", "QC",
                         stringr::str_extract(scores$Sample, "(?<=^(AH|PD|PL))\\d{2}"))
  scores$Coral <- factor(scores$Coral, levels = c("AH","PD","PL","QC"))
  scores$Temp  <- factor(scores$Temp,  levels = c("27","31","34","QC"))
  
  p <- ggplot(scores, aes(PC1, PC2, color = Coral, shape = Temp)) +
    geom_point(size = 3, alpha = 0.95) +
    stat_ellipse(
      data = subset(scores, Coral != "QC"),
      aes(group = interaction(Coral, Temp)),
      level = 0.95, linewidth = 0.6, alpha = 0.25, show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "AH" = "#B3D9FF",
      "PD" = "#1f8f86",
      "PL" = "#7C8CF8",
      "QC" = "grey80"
    )) +
    scale_shape_manual(
      values = c("27" = 16, "31" = 17, "34" = 15, "QC" = 3),
      breaks = c("27","31","34","QC"),
      labels = c("27℃","31℃","34℃","QC")
    ) +
    labs(title = if (is.null(title)) paste0("PCA of ", basename(path)) else title,
         x = pc1_lab, y = pc2_lab) +
    theme_classic(base_size = 13, base_family = "Times New Roman") +
    theme(
      text = element_text(family = "Times New Roman"),
      axis.text = element_text(family = "Times New Roman", size = 12, color = "black"),
      axis.title = element_text(family = "Times New Roman", size = 13, face = "bold", color = "black"),
      legend.title = element_text(family = "Times New Roman", size = 12, face = "bold", color = "black"),
      legend.text = element_text(family = "Times New Roman", size = 12, color = "black"),
      plot.title = element_text(family = "Times New Roman", face = "bold", hjust = 0.5, size = 14),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
      axis.line.y.left   = element_line(color = "black", linewidth = 0.6),
      legend.position = "right"
    )
  
  print(p)
  
  # 保存 PDF（嵌入字体）
  pdf_name <- paste0(
    tools::file_path_sans_ext(basename(path)),
    "_PCA_", format(Sys.time(), "%Y%m%d_%H%M"), ".pdf"
  )
  out_path <- file.path(dirname(path), pdf_name)
  ggsave(out_path, plot = p, width = 7, height = 6, device = cairo_pdf)
  message("✅ PCA PDF saved: ", out_path)
  
  invisible(list(scores = scores, pca = pca, var_exp = var_exp, pdf = out_path))
}

pca_res <- run_pca_from_file(
  path       = "D:/metabolites/metabolites_data.xlsx",
  sheet      = 1,                
  use_log2   = TRUE,
  pseudocount= 1,                
  scaling    = "autoscale",      
  title      = "PCA of Metabolomics (All samples)"
)
pca_res$pdf  # 会打印出保存的PDF路径

# ========= 运行（混合数据）=========
run_pipeline_mixed("D:/metabolites/metabolites_data.xlsx")
