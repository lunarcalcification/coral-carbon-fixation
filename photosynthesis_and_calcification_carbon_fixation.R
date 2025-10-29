# =========================
# 依赖
# =========================
# install.packages(c("tidyverse","readxl","janitor","rstatix","emmeans",
#                    "multcompView","broom","showtext","sysfonts","patchwork","ragg"))

library(tidyverse); library(readxl); library(janitor)
library(rstatix);    library(emmeans); library(multcompView); library(broom)
library(showtext);   library(sysfonts)

# =========================
# 0) 路径
# =========================
fpath   <- "F:/毕业论文及数据分析/data.xlsx"
out_dir <- dirname(fpath)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 稳妥保存 ----------
safe_ggsave <- function(filename, plot, width=10, height=5.5, dpi=300){
  ext <- tools::file_ext(filename)
  if (ext %in% c("png","tif","jpg","jpeg")) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggsave(filename, plot, width=width, height=height, dpi=dpi, device=ragg::agg_png)
    } else ggsave(filename, plot, width=width, height=height, dpi=dpi)
  } else if (ext == "pdf") {
    if (capabilities("cairo"))
      ggsave(filename, plot, width=width, height=height, device=grDevices::cairo_pdf)
    else
      ggsave(filename, plot, width=width, height=height, device=grDevices::pdf)
  } else ggsave(filename, plot, width=width, height=height, dpi=dpi)
  message("Saved: ", normalizePath(filename, winslash = "/"))
}

# =========================
# 0.5) 字体：Times New Roman（含斜体）+ 中文回退
# =========================
tnr_regular_paths <- c(
  "C:/Windows/Fonts/times.ttf",
  "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
  "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"
)
tnr_italic_paths <- c(
  "C:/Windows/Fonts/timesi.ttf",
  "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
  "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Italic.ttf"
)
tnr_reg <- tnr_regular_paths[file.exists(tnr_regular_paths)][1]
tnr_ita <- tnr_italic_paths[file.exists(tnr_italic_paths)][1]

if (!length(tnr_reg)) {
  showtext_auto()
  base_family <- theme_get()$text$family %||% ""
} else {
  font_add(family = "Times New Roman", regular = tnr_reg,
           italic  = if (length(tnr_ita)) tnr_ita else tnr_reg)
  showtext_auto()
  base_family <- "Times New Roman"
}

pick_cn_font <- function(){
  cand <- c(
    "C:/Windows/Fonts/msyh.ttc","C:/Windows/Fonts/simhei.ttf","C:/Windows/Fonts/simsun.ttc",
    "/System/Library/Fonts/STHeiti Light.ttc","/System/Library/Fonts/Songti.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/truetype/arphic/uming.ttc"
  )
  cand[file.exists(cand)][1]
}
cn_font <- pick_cn_font()
if (!is.null(cn_font)) try(font_add(family="CNfallback", regular=cn_font), silent=TRUE)

# =========================
# 1) 读取与列名“智能匹配”
# =========================
raw <- read_excel(fpath, sheet = 1) %>% clean_names(replace = c("\u00B5"="u", "\u03BC"="u"))

find_col <- function(nms, patterns, required=TRUE, label=""){
  for (pat in patterns) {
    idx <- which(stringr::str_detect(nms, stringr::regex(pat, ignore_case = TRUE)))
    if (length(idx) >= 1) return(nms[idx[1]])
  }
  if (required) stop(paste0("未能自动识别列：", label, "\n当前表头为：\n", paste(nms, collapse=", ")), call. = FALSE)
  NULL
}

nms <- names(raw)
col_species <- find_col(nms, c("^species$", "^sp$", "^taxon$"), TRUE,  "species")
col_treat   <- find_col(nms, c("^treat$", "^treatment$", "^light_dark$", "^phase$", "^condition$"),
                        TRUE,  "treatment")
col_temp    <- find_col(nms, c("^incu.*temp", "^incubation.*temp", "^temp$", "^temperature$", "temp_c", "temperature_c"),
                        TRUE,  "incubator temperature")
col_calrate <- find_col(nms, c("^cal[_]?rate", "^ncc", "calcification"), TRUE, "calcification rate")
col_dic     <- find_col(nms, c("^remove[_]?tdic$", "remove.*tdic", "^dic(_rate)?$", "^ncp_dic$", "^cr_dic$"),
                        required = FALSE, label = "DIC-based rate (NCP/CR)")

raw_std <- raw %>%
  rename(
    species   = !!sym(col_species),
    treatment = !!sym(col_treat),
    incu_temperature = !!sym(col_temp)
  ) %>%
  mutate(
    treatment = tolower(as.character(treatment)),
    treatment = recode(treatment, "day"="light", "night"="dark", "l"="light", "d"="dark"),
    species   = factor(species),
    temp_num  = readr::parse_number(as.character(incu_temperature)),
    incu_temp = factor(temp_num, levels = sort(unique(temp_num)))
  )

keep_species <- c("A.hyacinthus","P.damicornis","P.lutea")

# =========================
# 2) 数据拆分
# =========================
rate_col <- col_calrate
dat_cal <- raw_std %>% filter(species %in% keep_species)

has_dic <- !is.null(col_dic)
dat_prod <- NULL
if (has_dic) {
  dat_prod <- raw_std %>%
    mutate(remove_tdic = .data[[col_dic]]) %>%
    filter(treatment %in% c("light","dark"), species %in% keep_species)
}

# =========================
# 3) 24h TCCF（由 cal24 表示）
# =========================
id_cal <- case_when("id" %in% names(dat_cal) ~ "id",
                    "ID" %in% names(dat_cal) ~ "ID",
                    TRUE ~ NA_character_)
if (is.na(id_cal)) { dat_cal <- dat_cal %>% mutate(tmp_id=row_number()); id_cal <- "tmp_id" }

hours_light <- 12; hours_dark <- 12

daily_pair_cal <- dat_cal %>%
  filter(treatment %in% c("light","dark")) %>%
  select(!!rlang::sym(id_cal), species, incu_temp, treatment,
         rate = .data[[rate_col]]) %>%
  mutate(treatment=as.character(treatment)) %>%
  tidyr::pivot_wider(names_from=treatment, values_from=rate) %>%
  mutate(
    cal24 = light*hours_light + dark*hours_dark,
    cal_h = cal24/(hours_light+hours_dark)
  )

daily_summary_cal <- daily_pair_cal %>%
  group_by(species, incu_temp) %>%
  summarise(
    n_pair  = n(),
    mean_day = mean(cal24, na.rm=TRUE),
    se_day   = sd(cal24,   na.rm=TRUE) / sqrt(n_pair), # 使用SE
    .groups="drop"
  )

# =========================
# 4) NCP / CR / TPCF（24h 积分量）
# =========================
p_NCP24_line <- p_CR24_line <- p_TPCF24_line <- NULL
daily_pair_prod <- NULL; daily_summary_prod <- NULL

if (!is.null(dat_prod) && has_dic){
  id_prod <- case_when("id" %in% names(dat_prod) ~ "id",
                       "ID" %in% names(dat_prod) ~ "ID",
                       TRUE ~ NA_character_)
  if (is.na(id_prod)) { dat_prod <- dat_prod %>% mutate(tmp_id=row_number()); id_prod <- "tmp_id" }
  
  daily_pair_prod <- dat_prod %>%
    select(!!rlang::sym(id_prod), species, incu_temp, treatment,
           rate = remove_tdic) %>%
    mutate(treatment=as.character(treatment)) %>%
    tidyr::pivot_wider(names_from=treatment, values_from=rate) %>%
    mutate(
      NCP   = light,
      CR    = dark,
      GPP   = NCP + abs(CR),
      NCP24 = NCP * hours_light,
      CR24  = abs(CR) * hours_dark,
      GPP24 = NCP24 + CR24
    )
  
  daily_summary_prod <- daily_pair_prod %>%
    group_by(species, incu_temp) %>%
    summarise(
      n_pair     = n(),
      mean_NCP24 = mean(NCP24, na.rm=TRUE),  se_NCP24 = sd(NCP24, na.rm=TRUE) / sqrt(n_pair),  # 使用SE
      mean_CR24  = mean(CR24,  na.rm=TRUE),  se_CR24  = sd(CR24,  na.rm=TRUE) / sqrt(n_pair),
      mean_GPP24 = mean(GPP24, na.rm=TRUE),  se_GPP24 = sd(GPP24, na.rm=TRUE) / sqrt(n_pair),
      .groups="drop"
    )
}

# =========================
# 5) 统计分组字母：KW + Dunn（两组用 Wilcoxon）
# =========================
get_letters_by_species <- function(df, y_col, alpha = 0.05) {
  # 需要 df 至少包含：species, incu_temp, y_col
  dat <- df %>%
    dplyr::select(species, incu_temp, y = !!rlang::sym(y_col)) %>%
    dplyr::filter(!is.na(species), !is.na(incu_temp), !is.na(y)) %>%
    dplyr::mutate(incu_temp = as.factor(incu_temp))
  
  sp <- as.character(unique(dat$species))[1]
  lv <- levels(dat$incu_temp); k <- length(lv)
  
  # 只有一个温度水平
  if (k <= 1L) {
    letters_tbl <- tibble::tibble(
      species = sp,
      incu_temp = factor(lv, levels = lv),
      .group = "a"
    )
    log_tbl <- tibble::tibble(
      species = sp, k = k,
      route = "single level",
      global_p = NA_real_,
      posthoc = NA_character_
    )
    return(list(letters = letters_tbl, log = log_tbl))
  }
  
  # 全局：Kruskal–Wallis
  kw <- tryCatch(rstatix::kruskal_test(dat, y ~ incu_temp), error = function(e) NULL)
  global_p <- if (!is.null(kw)) kw$p else NA_real_
  
  # 两组：Wilcoxon 秩和检验
  if (k == 2L) {
    w <- tryCatch(rstatix::wilcox_test(dat, y ~ incu_temp, exact = FALSE), error = function(e) NULL)
    labs <- c("a","a")
    if (!is.null(w) && !is.na(w$p) && w$p < alpha) labs <- c("a","b")
    
    letters_tbl <- tibble::tibble(
      species   = sp,
      incu_temp = factor(lv, levels = lv),
      .group    = labs
    )
    log_tbl <- tibble::tibble(
      species  = sp, k = k,
      route    = "Wilcoxon (KW for 2 groups)",
      global_p = global_p,
      posthoc  = "none"
    )
    return(list(letters = letters_tbl, log = log_tbl))
  }
  
  # ≥3 组：Dunn 事后检验（Holm 校正）
  dres <- rstatix::dunn_test(dat, y ~ incu_temp, p.adjust.method = "holm") %>%
    dplyr::mutate(p.adj = ifelse(is.na(p.adj), 1, p.adj))
  
  pvec <- dres$p.adj
  names(pvec) <- paste(dres$group1, dres$group2, sep = "-")
  
  let <- multcompView::multcompLetters(pvec, threshold = alpha)$Letters
  
  letters_tbl <- tibble::tibble(
    species   = sp,
    incu_temp = factor(names(let), levels = lv),
    .group    = unname(let)
  )
  
  log_tbl <- tibble::tibble(
    species  = sp, k = k,
    route    = "Kruskal–Wallis",
    global_p = global_p,
    posthoc  = "Dunn (Holm)"
  )
  
  return(list(letters = letters_tbl, log = log_tbl))
}

# =========================
# 5.1 通用绘图参数
# =========================
dodge <- position_dodge(width = 0.35)
jit   <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.35)

# =========================
# 5.15 物种颜色映射（浅色系 + Times风格）
# =========================
col_map_species <- c(
  "A.hyacinthus" = "#B3D9FF",  # 浅蓝
  "P.damicornis" = "#1f8f86",  # 墨绿青
  "P.lutea"      = "#7C8CF8"   # 紫蓝
)

# =========================
# 5.2 画图函数（KW+Dunn 自动字母）
# =========================
make_line_with_letters <- function(raw_df, sum_df, y_raw, y_mean, y_se, ylab,
                                   letters_size = 7, letters_vjust = 0){
  # ① 为每个物种计算 KW + Dunn 字母
  letters_list <- raw_df %>%
    group_by(species) %>%
    group_map(~ {
      tmp <- .x %>%
        mutate(species = .y$species[[1]]) %>%           # 把分组键加回数据
        dplyr::select(species, incu_temp, !!rlang::sym(y_raw)) %>%
        dplyr::rename(y = !!rlang::sym(y_raw))
      get_letters_by_species(tmp, "y")                  # 返回 list(letters=..., log=...)
    })
  
  letters_df <- dplyr::bind_rows(lapply(letters_list, `[[`, "letters"))
  route_log  <- dplyr::bind_rows(lapply(letters_list, `[[`, "log")) %>%
    dplyr::mutate(global_p = signif(global_p, 3))       # 这里只保留 global_p
  print(route_log)                                      # 可注释掉
  
  # —— ② 字母放置高度（补 NA + 自动抬高上界）——
  label_df <- sum_df %>%
    select(species, incu_temp, mean = !!sym(y_mean), se = !!sym(y_se)) %>%  # 使用 SE 代替 SD
    left_join(letters_df, by = c("species","incu_temp")) %>%
    mutate(
      .group = tidyr::replace_na(.group, "a"),
      se0    = tidyr::replace_na(se, 0)
    ) %>%
    group_by(species) %>%
    mutate(y_pos = mean + se0 + max(se0, na.rm = TRUE) * 0.25) %>%
    ungroup()
  
  # 为了不被裁掉：把上界抬到能容纳字母
  y_top <- max(c(label_df$y_pos, sum_df[[y_mean]] + sum_df[[y_se]]), na.rm = TRUE) * 1.05
  
  p <- ggplot() +
    geom_point(
      data = raw_df %>% mutate(temp_fac = factor(as.numeric(as.character(incu_temp)), levels = c(27,31,34))),
      aes(x = temp_fac, y = .data[[y_raw]], color = species),
      position = jit, alpha = 0.40, size = 1.8
    ) +
    geom_errorbar(
      data = sum_df %>% mutate(temp_fac = factor(as.numeric(as.character(incu_temp)), levels = c(27,31,34))),
      aes(x = temp_fac, ymin = .data[[y_mean]] - .data[[y_se]], ymax = .data[[y_mean]] + .data[[y_se]], color = species),  # 使用 SE
      width = 0.28, linewidth = 0.9, position = dodge
    ) +
    geom_line(
      data = sum_df %>% mutate(temp_fac = factor(as.numeric(as.character(incu_temp)), levels = c(27,31,34))),
      aes(x = temp_fac, y = .data[[y_mean]], color = species, group = species),
      linewidth = 1.2, alpha = 0.9, position = dodge
    ) +
    geom_point(
      data = sum_df %>% mutate(temp_fac = factor(as.numeric(as.character(incu_temp)), levels = c(27,31,34))),
      aes(x = temp_fac, y = .data[[y_mean]], color = species),
      size = 3, position = dodge
    ) +
    geom_text(
      data = label_df %>% mutate(temp_fac = factor(as.numeric(as.character(incu_temp)), levels = c(27,31,34))),
      aes(x = temp_fac, y = y_pos, label = .group, color = species),
      position = dodge, vjust = 0, size = letters_size, show.legend = FALSE
    ) +
    scale_color_manual(values = col_map_species) +
    scale_x_discrete(labels = c("27"="27\u2103","31"="31\u2103","34"="34\u2103")) +
    scale_y_continuous(name = ylab, expand = expansion(mult = c(0.05, 0.18))) +
    coord_cartesian(ylim = c(NA, y_top)) +
    labs(x = NULL) +
    theme_bw(base_size = 20) +
    theme(
      text = element_text(family = base_family),
      axis.text.x = element_text(size = 20, vjust = 0.5),
      axis.text.y = element_text(size = 20),
      legend.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(family = base_family, face = "italic", size = 18),
      panel.grid.minor = element_blank()
    )
  attr(p, "route_log") <- route_log
  p
}

# =========================
# 5.3 每日光/暗钙化（LC / DC）
# =========================
raw_light_cal24 <- dat_cal %>%
  filter(treatment == "light", species %in% keep_species) %>%
  mutate(LC = .data[[rate_col]] * hours_light)

raw_dark_cal24 <- dat_cal %>%
  filter(treatment == "dark", species %in% keep_species) %>%
  mutate(DC = .data[[rate_col]] * hours_dark)

daily_light_cal24 <- raw_light_cal24 %>%
  group_by(species, incu_temp) %>%
  summarise(
    mean_LC = mean(LC, na.rm = TRUE),
    se_LC   = sd(LC, na.rm = TRUE) / sqrt(n()),  # 计算标准误差
    n_LC    = n(),
    .groups = "drop"
  )

daily_dark_cal24 <- raw_dark_cal24 %>%
  group_by(species, incu_temp) %>%
  summarise(
    mean_DC = mean(DC, na.rm = TRUE),
    se_DC   = sd(DC, na.rm = TRUE) / sqrt(n()),  # 计算标准误差
    n_DC    = n(),
    .groups = "drop"
  )

# 作图（带 KW+Dunn 字母）
p_LC <- make_line_with_letters(
  raw_df = raw_light_cal24,
  sum_df = daily_light_cal24,
  y_raw  = "LC",
  y_mean = "mean_LC",
  y_se   = "se_LC",  # 使用 SE
  ylab   = expression(paste("LC (", mu, "mol C ", cm^{-2}, " ", day^{-1}, ")"))
)

p_DC <- make_line_with_letters(
  raw_df = raw_dark_cal24,
  sum_df = daily_dark_cal24,
  y_raw  = "DC",
  y_mean = "mean_DC",
  y_se   = "se_DC",  # 使用 SE
  ylab   = expression(paste("DC (", mu, "mol C ", cm^{-2}, " ", day^{-1}, ")"))
)

# =========================
# 6) 生成“折线 + 字母”图（TCCF / 以及可选的 NCP/CR/TPCF）
# =========================
p_TCCF24_line <- make_line_with_letters(
  raw_df = daily_pair_cal,
  sum_df = daily_summary_cal,
  y_raw = "cal24",
  y_mean = "mean_day",
  y_se = "se_day",  # 使用 SE
  ylab = expression(paste("TCCF (", mu, "mol C ", cm^{-2}, " ", day^{-1}, ")"))
)

if (!is.null(daily_pair_prod)) {
  p_NCP24_line <- make_line_with_letters(
    raw_df = daily_pair_prod, sum_df = daily_summary_prod,
    y_raw = "NCP24", y_mean = "mean_NCP24", y_se = "se_NCP24",  # 使用 SE
    ylab = expression(paste("NCP (", mu, "mol C ", cm^{-2}, " ", day^{-1}, ")"))
  )
  p_CR24_line <- make_line_with_letters(
    raw_df = daily_pair_prod, sum_df = daily_summary_prod,
    y_raw = "CR24", y_mean = "mean_CR24", y_se = "se_CR24",  # 使用 SE
    ylab = expression(paste("CR (", mu, "mol C ", cm^{-2}, " ", day^{-1}, ")"))
  )
  p_TPCF24_line <- make_line_with_letters(
    raw_df = daily_pair_prod, sum_df = daily_summary_prod,
    y_raw = "GPP24", y_mean = "mean_GPP24", y_se = "se_GPP24",  # 使用 SE
    ylab = expression(paste("TPCF (", mu, "mol C ", cm^{-2}, " ", day^{-1}, ")"))
  )
}

# =========================
# 7) 导出
# =========================
message("Output dir: ", normalizePath(out_dir, winslash = "/"))

safe_ggsave(file.path(out_dir, "TCCF_24h_line_mean_se_letters.pdf"),  p_TCCF24_line)
safe_ggsave(file.path(out_dir, "Light_Calcification_rate_line_mean_se_letters.pdf"),  p_LC)
safe_ggsave(file.path(out_dir, "Night_Calcification_rate_line_mean_se_letters.pdf"),  p_DC)
safe_ggsave(file.path(out_dir, "NCP_24h_line_mean_se_letters.pdf"),  p_NCP24_line)
safe_ggsave(file.path(out_dir, "CR_24h_line_mean_se_letters.pdf"),   p_CR24_line)
safe_ggsave(file.path(out_dir, "TPCF_24h_line_mean_se_letters.pdf"),  p_TPCF24_line)


# =========================
# 8) 预览（可选）
# =========================
print(p_CR24_line)
print(p_NCP24_line)
print(p_DC)
print(p_LC)
print(p_TPCF24_line)
print(p_TCCF24_line)
