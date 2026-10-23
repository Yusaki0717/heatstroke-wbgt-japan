# =============================================================================
# 04_generate_grouped_figures.R -- Figures 1-7 and the TOC graphic from the
# locked summer-grouped models. Written to output/figures/.
#
# Figure 1: pooled exposure-response curves (lag 0-5), three periods
# Figure 2: Tokyo lag-response at the 95th percentile
# Figure 3: national daily time series, broken between summers
# Figure 4: post-pandemic BLUP forest plot
# Figure 5: Table 4 sensitivity estimates (incl. severe / non-severe strata)
# Figure 6: choropleth maps (a-c periods, d % change post vs pre)
# Figure 7: national attributable burden + post-pandemic prefecture AF
# TOC:      annualized attributable EAD with pooled RR at P95
#
# Prerequisites: 01_core_reanalysis.R and 02_attributable_reanalysis.R.
# Extra packages: ggplot2, tidyr, tibble, patchwork, sf, rnaturalearth
# (+ rnaturalearthhires for ne_states), viridis, cowplot.
# =============================================================================

source("R/00_setup.R")

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(sf)
  library(rnaturalearth)
  library(viridis)
  library(cowplot)
  library(grid)
})

start_log("04_generate_grouped_figures")

obj <- readRDS(file.path(OUT_DIR, "grouped_stage1_main_severity.rds"))
MAIN <- obj$MAIN
t2 <- read_out("table2_grouped_exact.csv")
t4 <- read_out("table4_grouped_exact.csv")
pref_rr <- read_out("prefecture_blup_p95_grouped.csv")
rr_change <- read_out("prefecture_change_grouped.csv")
af_pref <- read_out("attributable_prefecture_all_grouped.csv")
af_sum <- read_out("attributable_summary_grouped.csv")

period_labels <- c(
  pre_pandemic = "Pre-pandemic\n(2015-2019)",
  pandemic = "Pandemic\n(2020-2021)",
  post_pandemic = "Post-pandemic\n(2022-2023)"
)
period_cols <- c("#2166AC", "#B2182B", "#4DAF4A")

# Figure 1: pooled exposure-response curves.
grid_heat <- sort(unique(c(seq(wbgt_cen, WBGT_RANGE[2], length.out = 120),
                           Q_REPORT[["P95"]])))
B <- ns(grid_heat, knots = knots_var, Boundary.knots = WBGT_RANGE)
B0 <- as.numeric(ns(wbgt_cen, knots = knots_var, Boundary.knots = WBGT_RANGE))
Bc <- sweep(B, 2, B0, "-")
curve_list <- list()
for (p in periods) {
  po <- pool_grouped(MAIN, p)
  lr <- as.numeric(Bc %*% coef(po$fit))
  se <- sqrt(diag(Bc %*% vcov(po$fit) %*% t(Bc)))
  curve_list[[p]] <- data.frame(wbgt = grid_heat, rr = exp(lr),
                                 lo = exp(lr - 1.96 * se), hi = exp(lr + 1.96 * se),
                                 period = p)
}
curve_df <- bind_rows(curve_list) %>%
  mutate(period = factor(period, levels = periods, labels = period_labels[periods]))
fig1 <- ggplot(curve_df, aes(wbgt, rr, colour = period, fill = period)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = Q_REPORT[["P95"]], linetype = "dotted", colour = "grey35") +
  scale_colour_manual(values = period_cols) + scale_fill_manual(values = period_cols) +
  coord_cartesian(xlim = c(wbgt_cen, WBGT_RANGE[2])) +
  labs(x = expression("Daily maximum WBGT ("*degree*"C)"),
       y = "Cumulative Relative Risk",
       title = "Pooled exposure-response: WBGT and heatstroke ambulance dispatches",
       subtitle = expression("Reference: WBGT = 27.5"*degree*"C (median); 47 prefectures"),
       colour = "Period", fill = "Period") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(FIG_OUT, "Figure1_grouped.png"), fig1, width = 12, height = 7, dpi = 300)
ggsave(file.path(FIG_OUT, "Figure1_grouped.pdf"), fig1, width = 12, height = 7)

# The curve must pass exactly through the Table 2 value at the 95th percentile.
for (p in periods) {
  z <- curve_list[[p]][which.min(abs(curve_list[[p]]$wbgt - Q_REPORT[["P95"]])), ]
  expected <- t2$rr[t2$percentile == "P95" & t2$period == p]
  stopifnot(abs(z$rr - expected) < 1e-8)
}
cat("FIG1|PASS\n")

# Figure 2: Tokyo lag-response with summer grouping.
lag_list <- list(); tok_cum <- numeric()
for (p in periods) {
  d <- dat %>% filter(pref_code == 13, period == p, !is.na(wbgt_max)) %>% arrange(date)
  cb <- crossbasis(d$wbgt_max, lag = 5,
                   argvar = list(fun = "ns", knots = knots_var, Boundary.knots = WBGT_RANGE),
                   arglag = list(fun = "ns", knots = knots_lag),
                   group = factor(d$year))
  mo <- gnm(total ~ cb + factor(dow) + holiday, eliminate = d$stratum,
            family = quasipoisson(), data = d)
  cp <- crosspred(cb, mo, at = Q_REPORT[["P95"]], cen = wbgt_cen)
  lag_list[[p]] <- data.frame(lag = 0:5, rr = cp$matRRfit[1, ],
                              lo = cp$matRRlow[1, ], hi = cp$matRRhigh[1, ], period = p)
  tok_cum[p] <- prod(cp$matRRfit[1, ])
}
lag_df <- bind_rows(lag_list) %>%
  mutate(period = factor(period, levels = periods, labels = period_labels[periods]))
fig2 <- ggplot(lag_df, aes(lag, rr, colour = period, fill = period)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_x_continuous(breaks = 0:5) +
  scale_colour_manual(values = period_cols) + scale_fill_manual(values = period_cols) +
  labs(x = "Lag (days)", y = "Relative Risk",
       title = expression("Lag-response at WBGT 95th percentile (31.9"*degree*"C) - Tokyo"),
       subtitle = expression("Reference: WBGT = 27.5"*degree*"C"),
       colour = "Period", fill = "Period") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(FIG_OUT, "Figure2_grouped.png"), fig2, width = 8, height = 6, dpi = 300)
ggsave(file.path(FIG_OUT, "Figure2_grouped.pdf"), fig2, width = 8, height = 6)
tok_check <- read_out("tokyo_cumulative_p95_grouped.csv")
stopifnot(max(abs(tok_cum[periods] - tok_check$rr[match(periods, tok_check$period)])) < 1e-6)
cat("FIG2|PASS|", paste(sprintf("%.6f", tok_cum[periods]), collapse = ","), "\n", sep = "")

# Figure 3: descriptive national time series, explicitly broken between summers.
daily <- dat %>% group_by(date, period, year) %>%
  summarise(total = sum(total), wbgt_mean = mean(wbgt_max, na.rm = TRUE), .groups = "drop") %>%
  mutate(series = interaction(period, year, drop = TRUE),
         period_f = factor(period, levels = periods, labels = period_labels[periods]))
fig3a <- ggplot(daily, aes(date, total, colour = period_f, group = series)) +
  geom_line(alpha = 0.8, linewidth = 0.35) +
  scale_colour_manual(values = period_cols) +
  labs(x = NULL, y = "Daily national dispatches",
       title = "Heatstroke ambulance dispatches, Japan 2015-2023") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")
fig3b <- ggplot(daily, aes(date, wbgt_mean, colour = period_f, group = series)) +
  geom_line(alpha = 0.8, linewidth = 0.35) +
  scale_colour_manual(values = period_cols) +
  labs(x = NULL, y = expression("Mean daily WBGTmax ("*degree*"C)"), colour = "Period") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
fig3 <- fig3a / fig3b + plot_layout(heights = c(2, 1))
ggsave(file.path(FIG_OUT, "Figure3_year_gaps.png"), fig3, width = 12, height = 7, dpi = 300)
ggsave(file.path(FIG_OUT, "Figure3_year_gaps.pdf"), fig3, width = 12, height = 7)
cat("FIG3|PASS|year_groups=", length(unique(daily$series)), "\n", sep = "")

# Figure 4: post-pandemic BLUP forest plot.
fp <- pref_rr[pref_rr$period == "post_pandemic", ]
fp <- fp[order(fp$rr), ]; fp$pref_name <- factor(fp$pref_name, levels = fp$pref_name)
fig4 <- ggplot(fp, aes(rr, pref_name)) +
  geom_point(size = 1.5) + geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(x = expression("Cumulative RR at WBGT P95 (31.9"*degree*"C) vs median"), y = NULL,
       title = "Post-pandemic (2022-2023)") +
  theme_minimal(base_size = 9) + theme(axis.text.y = element_text(size = 7))
ggsave(file.path(FIG_OUT, "Figure4_grouped.png"), fig4, width = 6, height = 10, dpi = 300)
ggsave(file.path(FIG_OUT, "Figure4_grouped.pdf"), fig4, width = 6, height = 10)
cat(sprintf("FIG4|PASS|min=%.6f|max=%.6f\n", min(fp$rr), max(fp$rr)))

# Figure 5: sensitivity estimates at 300 dpi.
f5 <- t4
f5$analysis <- factor(f5$analysis,
  levels = c("Main analysis (lag 0-5)", "Lag 0-7", "Excluding 2020",
             "Elderly only (>=65)", "Severe cases", "Non-severe cases"),
  labels = c("Main analysis", "Lag 0-7", "Exclude 2020", "Elderly only",
             "Severe cases", "Non-severe cases"))
f5$period <- factor(f5$period, levels = periods,
                    labels = c("Pre-pandemic", "Pandemic", "Post-pandemic"))
fig5 <- ggplot(f5, aes(analysis, rr, colour = period, group = period)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6),
                width = 0.25, linewidth = 0.9) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  scale_colour_manual(values = period_cols, name = "Period") +
  scale_x_discrete(limits = rev) + coord_flip() +
  labs(y = expression("Cumulative RR at WBGT P95 (31.9"*degree*"C)"), x = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.major.x = element_line(colour = "grey90"),
        panel.grid.minor = element_blank())
ggsave(file.path(FIG_OUT, "Figure5_grouped.png"), fig5, width = 8, height = 5, dpi = 300)
ggsave(file.path(FIG_OUT, "Figure5_grouped.pdf"), fig5, width = 8, height = 5)
stopifnot(nrow(f5) == 18L); cat("FIG5|PASS\n")

# Shared Japan geometry for Figures 6 and 7 (needs rnaturalearthhires).
japan_sf <- ne_states(country = "japan", returnclass = "sf")
japan_sf$pref_code <- as.integer(sub("^JP-", "", japan_sf$iso_3166_2))
stopifnot(nrow(japan_sf) == 47, setequal(japan_sf$pref_code, 1:47))

# Figure 6: compact 2x2 layout with two legends; all 47 prefectures estimable,
# so no prefecture is greyed out.
rr_lim <- c(floor(min(pref_rr$rr) * 2) / 2, ceiling(max(pref_rr$rr) * 2) / 2)
ch_lim <- c(floor(min(rr_change$pct_change) / 5) * 5,
            ceiling(max(rr_change$pct_change) / 5) * 5)
make_rr_map <- function(p, title) {
  z <- merge(japan_sf, pref_rr[pref_rr$period == p, ], by = "pref_code", all.x = TRUE)
  ggplot(z) + geom_sf(aes(fill = rr), colour = "white", linewidth = 0.15) +
    scale_fill_viridis(option = "inferno", direction = -1, name = "RR",
                       limits = rr_lim, na.value = "grey85") +
    labs(title = title) + theme_void(base_size = 9) +
    theme(plot.title = element_text(hjust = .5, face = "bold", size = 10))
}
p6a <- make_rr_map("pre_pandemic", "(a) Pre-pandemic (2015-2019)")
p6b <- make_rr_map("pandemic", "(b) Pandemic (2020-2021)")
p6c <- make_rr_map("post_pandemic", "(c) Post-pandemic (2022-2023)")
map_change <- merge(japan_sf, rr_change, by = "pref_code", all.x = TRUE)
p6d <- ggplot(map_change) +
  geom_sf(aes(fill = pct_change), colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = mean(ch_lim), name = "Change (%)", limits = ch_lim,
                       na.value = "grey85") +
  labs(title = "(d) RR change (post vs pre)") + theme_void(base_size = 9) +
  theme(plot.title = element_text(hjust = .5, face = "bold", size = 10))
leg_rr <- cowplot::get_legend(p6c + theme(legend.position = "bottom") +
                               guides(fill = guide_colorbar(title.position = "top")))
leg_ch <- cowplot::get_legend(p6d + theme(legend.position = "bottom") +
                               guides(fill = guide_colorbar(title.position = "top")))
maps <- cowplot::plot_grid(p6a + theme(legend.position = "none"),
                           p6b + theme(legend.position = "none"),
                           p6c + theme(legend.position = "none"),
                           p6d + theme(legend.position = "none"), ncol = 2)
legs <- cowplot::plot_grid(leg_rr, leg_ch, ncol = 2)
fig6 <- cowplot::plot_grid(maps, legs, ncol = 1, rel_heights = c(1, 0.10))
ggsave(file.path(FIG_OUT, "Figure6_grouped_2x2.png"), fig6, width = 14, height = 12, dpi = 300)
ggsave(file.path(FIG_OUT, "Figure6_grouped_2x2.pdf"), fig6, width = 14, height = 12)
cat(sprintf("FIG6|PASS|change=%.6f,%.6f\n", min(rr_change$pct_change), max(rr_change$pct_change)))

# Figure 7: national burden plus post-pandemic prefecture AF.
ab <- af_sum %>% filter(outcome == "all_severity", reference == "P50") %>%
  mutate(non_attr = total_ead - an,
         period_lab = factor(period_labels[period], levels = period_labels[periods]),
         annual_an_round = round(an_annual))
long_ab <- ab %>% select(period_lab, Attributable = an, `Non-attributable` = non_attr) %>%
  pivot_longer(-period_lab, names_to = "type", values_to = "count")
p7a <- ggplot(long_ab, aes(period_lab, count / 1000, fill = type)) +
  geom_col(width = 0.65) +
  geom_text(data = ab, aes(period_lab, total_ead / 1000 + 8,
                           label = sprintf("AF = %.1f%%", af), fill = NULL),
            size = 4.2, fontface = "bold") +
  geom_text(data = ab, aes(period_lab, (an / 2) / 1000,
                           label = sprintf("~%s/yr", format(annual_an_round, big.mark = ",")),
                           fill = NULL), size = 3.6, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("Attributable" = "#D73027", "Non-attributable" = "#4575B4"),
                    name = NULL) +
  labs(x = NULL, y = "Total EAD (thousands)", tag = "(a)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 15),
        panel.grid.minor = element_blank())
af_pref$af <- 100 * af_pref$an_p50 / af_pref$total_ead
map_af <- merge(japan_sf, af_pref[af_pref$period == "post_pandemic", ],
                by = "pref_code", all.x = TRUE)
p7b <- ggplot(map_af) + geom_sf(aes(fill = af), colour = "white", linewidth = 0.2) +
  scale_fill_viridis(option = "plasma", direction = -1, name = "AF (%)", na.value = "grey80") +
  labs(tag = "(b)", title = "Post-pandemic (2022-2023)") +
  theme_minimal(base_size = 13) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        plot.tag = element_text(face = "bold", size = 15),
        plot.title = element_text(size = 11, hjust = 0.5), legend.position = "right")
fig7 <- p7a + p7b + plot_layout(widths = c(1.3, 1))
ggsave(file.path(FIG_OUT, "Figure7_grouped.png"), fig7, width = 14, height = 6.5, dpi = 300)
ggsave(file.path(FIG_OUT, "Figure7_grouped.pdf"), fig7, width = 14, height = 6.5)
cat("FIG7|PASS\n")

# TOC graphic at the required physical size and 300 dpi.
t2p <- t2[t2$percentile == "P95", ]
toc_df <- ab %>%
  mutate(period_short = factor(c("2015-2019", "2020-2021", "2022-2023"),
                               levels = c("2015-2019", "2020-2021", "2022-2023")),
         rr = t2p$rr[match(period, t2p$period)])
sc <- max(toc_df$an_annual / 1000) / max(toc_df$rr)
toc <- ggplot(toc_df, aes(period_short)) +
  geom_col(aes(y = an_annual / 1000), width = .55, fill = "grey78") +
  geom_line(aes(y = rr * sc, group = 1), colour = "#B2182B", linewidth = 0.9) +
  geom_point(aes(y = rr * sc), colour = "#B2182B", size = 2.2) +
  scale_y_continuous(name = "Attributable EAD\n(thousands/year)",
                     sec.axis = sec_axis(~ . / sc, name = "Cumulative RR")) +
  labs(x = NULL, title = "Falling heat risk, rising heat burden") +
  theme_minimal(base_size = 7.5) +
  theme(plot.title = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank(), legend.position = "none",
        axis.title.y.right = element_text(colour = "#B2182B"),
        axis.text.y.right = element_text(colour = "#B2182B"),
        plot.margin = margin(3, 5, 3, 5))
ggsave(file.path(FIG_OUT, "TOC_grouped.png"), toc, width = 7, height = 3.77, dpi = 300)
ggsave(file.path(FIG_OUT, "TOC_grouped.tiff"), toc, width = 7, height = 3.77,
       dpi = 300, compression = "lzw")
cat("TOC|PASS|2100x1131@300dpi\n")

expected_png <- c("Figure1_grouped.png", "Figure2_grouped.png", "Figure3_year_gaps.png",
                  "Figure4_grouped.png", "Figure5_grouped.png", "Figure6_grouped_2x2.png",
                  "Figure7_grouped.png", "TOC_grouped.png")
stopifnot(all(file.exists(file.path(FIG_OUT, expected_png))))
end_log("04_generate_grouped_figures")
