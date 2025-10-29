# =====================================================
# 📦 安装/加载
# =====================================================
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("ropensci/rnaturalearthhires")
}

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(showtext)

# 字体（Times New Roman 用于地图角标）
font_add(
  family = "Times New Roman",
  regular = "C:/Windows/Fonts/times.ttf",
  bold = "C:/Windows/Fonts/timesbd.ttf",
  italic = "C:/Windows/Fonts/timesi.ttf",
  bolditalic = "C:/Windows/Fonts/timesbi.ttf"
)
showtext_auto()


# Define the file path
file_path <- "C:/Users/Administrator/Desktop/Data_1.xlsx"

# Read the second sheet (Sheet 2) from the Excel file
data_sheet2 <- read_excel(file_path, sheet = 2)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =====================================================
# ① Sheet1：温度时间序列（折线图）
# =====================================================
dat_temp <- read_excel(file.path(out_dir, "hobo.xlsx"), sheet = "Sheet1")
colnames(dat_temp) <- c("Time", "T27", "T31", "T34")

dat_temp$Time <- as.POSIXct(
  dat_temp$Time, format = "%m/%d/%y %p%I时%M分%S秒", tz = "Asia/Shanghai"
)

df_long <- dat_temp %>%
  pivot_longer(cols = c(T27, T31, T34),
               names_to = "Tank", values_to = "Temperature")

colors <- c("T27"="#A8DADC","T31"="#F9C784","T34"="#F4978E")

p_temp <- ggplot(df_long, aes(Time, Temperature, color = Tank)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors,
                     labels = c("27°C 缸","31°C 缸","34°C 缸"),
                     name = "实验缸") +
  scale_y_continuous(limits = c(26,35), breaks = 26:35,
                     labels = function(x) paste0(x,"°C")) +
  scale_x_datetime(
    limits = as.POSIXct(c("2024-09-24 00:00:00","2024-10-08 00:00:00"),
                        tz = "Asia/Shanghai"),
    date_breaks = "1 day", date_labels = "%m-%d"
  ) +
  labs(x=NULL, y=NULL, title=NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.border  = element_blank(),
    axis.line.x   = element_line(color = "grey40", linewidth = 0.6),
    axis.line.y   = element_line(color = "grey40", linewidth = 0.6),
    axis.ticks.x  = element_line(color = "grey40", linewidth = 0.5),
    axis.ticks.y  = element_line(color = "grey40", linewidth = 0.5),
    axis.ticks.length = unit(3, "pt"),
    axis.text.x = element_text(vjust = 0.5)
  )

ggsave(file.path(out_dir, "hobo温度变化图.pdf"),
       plot = p_temp, width = 10, height = 6, units = "in")

# =====================================================
# ② Sheet2：Fv/Fm 柱形图（加入 27℃恒温；三处理×两日期成对比较）
# =====================================================
raw <- read_excel("C:/Users/Administrator/Desktop/Data_1.xlsx", sheet = 3)

df <- raw %>%
  rename(
    species   = !!names(raw)[1],
    treatment = !!names(raw)[2],
    date_raw  = !!names(raw)[3],
    FvFm      = !!names(raw)[4]
  ) %>%
  mutate(
    date = as.Date(as.character(date_raw)),
    # 标准化处理名称：27℃恒温 / 31℃ 升温 / 34℃ 升温
    treatment = case_when(
      str_detect(treatment, "27") ~ "27℃ 恒温",
      str_detect(treatment, "31") ~ "31℃ 升温",
      str_detect(treatment, "34") ~ "34℃ 升温",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(treatment),
         date %in% as.Date(c("2024-09-28","2024-10-05"))) %>%
  mutate(
    date_lab = format(date, "%m-%d"),
    # 设定 x 轴顺序：三处理依次，每个处理内是 09-28, 10-05
    xgroup = factor(paste(treatment, date),
                    levels = c(
                      paste("27℃ 恒温", as.Date("2024-09-28")),
                      paste("27℃ 恒温", as.Date("2024-10-05")),
                      paste("31℃ 升温", as.Date("2024-09-28")),
                      paste("31℃ 升温", as.Date("2024-10-05")),
                      paste("34℃ 升温", as.Date("2024-09-28")),
                      paste("34℃ 升温", as.Date("2024-10-05"))
                    )),
    treatment = factor(treatment, levels = c("27℃ 恒温","31℃ 升温","34℃ 升温"))
  )

# 汇总（均值、SE）
sumdf <- df %>%
  group_by(species, treatment, xgroup, date_lab) %>%
  summarise(
    n = dplyr::n(),
    mean = mean(FvFm, na.rm = TRUE),
    sd   = sd(FvFm,   na.rm = TRUE),
    se   = sd/sqrt(n),
    .groups = "drop"
  )

# 星号函数
p_stars <- function(p){
  ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# 在 每个物种 × 处理 内：比较 09-28 vs 10-05
test_by_grp <- df %>%
  group_by(species, treatment) %>%
  summarise(
    p_value = {
      g1 <- FvFm[date == as.Date("2024-09-28")]
      g2 <- FvFm[date == as.Date("2024-10-05")]
      use_t <- FALSE
      if(length(g1) >= 5 && length(g2) >= 5){
        s1 <- tryCatch(shapiro.test(g1)$p.value, error = function(e) NA_real_)
        s2 <- tryCatch(shapiro.test(g2)$p.value, error = function(e) NA_real_)
        if(!is.na(s1) && !is.na(s2) && s1 > 0.05 && s2 > 0.05) use_t <- TRUE
      }
      if(use_t) t.test(g1, g2, var.equal = FALSE)$p.value
      else      wilcox.test(g1, g2, exact = FALSE)$p.value
    },
    .groups = "drop"
  ) %>%
  mutate(
    label = p_stars(p_value),
    # 三个处理对应三组括号的 x 位置：27(1,2), 31(3,4), 34(5,6)
    x1 = case_when(
      treatment == "27℃ 恒温" ~ 1,
      treatment == "31℃ 升温" ~ 3,
      treatment == "34℃ 升温" ~ 5
    ),
    x2 = x1 + 1,
    xm = (x1 + x2)/2
  )

# —— 统一三面板的括号高度（更靠近柱子）——
bracket_gap <- 0.008   # 括号距离柱顶(含误差棒)
bracket_h   <- 0.0035  # 括号两端竖线长度
star_off    <- 0.0055  # 星号相对括号上移

y_global <- max(sumdf$mean + sumdf$se, na.rm = TRUE) + bracket_gap
bracket_y <- tibble::tibble(species = unique(sumdf$species), y = y_global)
test_by_grp <- left_join(test_by_grp, bracket_y, by = "species")
y_top <- y_global + star_off + 0.01

# 颜色：浅→深（27 恒温最浅，34 升温最深）
cols_bar <- c("27℃ 恒温"="#FDE0DF", "31℃ 升温"="#F8B9B6", "34℃ 升温"="#E67C73")

p_bar <- ggplot(sumdf, aes(x = xgroup, y = mean, fill = treatment)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.18, linewidth = 0.4) +
  facet_wrap(~ species, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    values = cols_bar,
    name   = NULL,
    labels = c("Control", "Moderate heat stress", "Severe heat stress")
  ) +
  guides(fill = guide_legend(title.position = "top")) +
  theme(
    legend.title = element_text(family = "Times New Roman"),
    legend.text  = element_text(family = "Times New Roman")
  ) +
  scale_x_discrete(labels = c(
    "27℃ 恒温 2024-09-28"="T0","27℃ 恒温 2024-10-05"="T1",
    "31℃ 升温 2024-09-28"="T0","31℃ 升温 2024-10-05"="T1",
    "34℃ 升温 2024-09-28"="T0","34℃ 升温 2024-10-05"="T1"))+
  scale_y_continuous(breaks = seq(0.45, y_top, by = 0.05),
                     minor_breaks = NULL, expand = expansion(mult = c(0,0))) +
  coord_cartesian(ylim = c(0.45, y_top), clip = "on") +
  labs(x=NULL, y="Fv/Fm", title=NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.text = element_text(family = "Times New Roman"),
    strip.text = element_text(family = "Times New Roman", face = "italic"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line.x = element_line(color="grey40", linewidth=0.6),
    axis.line.y = element_line(color="grey40", linewidth=0.6),
    axis.ticks.x = element_line(color="grey40", linewidth=0.5),
    axis.ticks.y = element_line(color="grey40", linewidth=0.5),
    axis.ticks.length = unit(3, "pt")
  ) +
  # 统一高度的括号 + 更靠近柱形图
  geom_segment(data = test_by_grp,
               aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = test_by_grp,
               aes(x = x1, xend = x1, y = y, yend = y - bracket_h),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = test_by_grp,
               aes(x = x2, xend = x2, y = y, yend = y - bracket_h),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = test_by_grp,
            aes(x = xm, y = y + star_off, label = label),
            size = 5, inherit.aes = FALSE)

ggsave(file.path(out_dir, "FvFm_三处理成对柱状图_统一括号.pdf"),
       plot = p_bar, width = 12, height = 5, units = "in")




# =====================================================
# ③ 海南采样点地图（10m，五角星 + 两处英文标注）
# =====================================================
coast <- ne_coastline(scale = 10, returnclass = "sf")
china_prov <- rnaturalearth::ne_download(
  scale = 10, type = "states", category = "cultural", returnclass = "sf"
)

name_cols <- intersect(c("name","name_en","gn_name","name_local"), names(china_prov))
key <- do.call(paste, c(china_prov[name_cols], list(sep = " ")))
hainan <- china_prov[grepl("Hainan|海南", key, ignore.case = TRUE), ]

lon <- 109 + 28/60   # 109°28′E
lat <- 18 + 12/60    # 18°12′N

xlim <- c(108, 111.3)
ylim <- c(17, 20.2)
x_mid <- mean(xlim)
y_3q  <- ylim[1] + 0.75 * diff(ylim)

p_map <- ggplot() +
  geom_sf(data = coast,  color = "grey60", linewidth = 0.3) +
  geom_sf(data = hainan, fill = "#E9ECEF", color = "grey30", linewidth = 0.6) +
  annotate("text", x = lon, y = lat, label = "\u2606",
           color = "red", size = 7, fontface = "bold") +
  annotate("text", x = lon, y = lat - 0.10, label = "Coral sampling site",
           color = "red", size = 4, family = "Times New Roman", hjust = 0.5, vjust = 1) +
  annotate("text", x = 110.3, y = 17.25, label = "South China Sea",
           color = "grey25", size = 6, family = "Times New Roman") +
  annotate("text", x = x_mid, y = y_3q, label = "Hai Nan Island",
           color = "grey25", size = 6, family = "Times New Roman") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3, line_width = 0.4) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(x=NULL, y=NULL, title=NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid   = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text    = element_text(color = "black"),
    axis.title   = element_blank(),
    plot.title   = element_blank(),
    axis.ticks   = element_line(color = "grey40", linewidth = 0.4)
  )

ggsave(file.path("C:/Users/Administrator/Desktop/", "采样点位置图_五角星精美版_标注.pdf"),
       plot = p_map, width = 7, height = 6, units = "in")

# ===== 完成 =====
