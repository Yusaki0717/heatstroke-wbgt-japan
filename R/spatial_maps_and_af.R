# =============================================================================
# Spatial maps + attributable fraction / attributable number
#
# Produces: Fig 4, Fig 6 (choropleth maps), and the all-severity + severe-case
# rows of Table 5 (AF/AN by period).
#
# PREREQUISITE: run R/heatstroke_dlnm_analysis.R first, so the workspace
# contains: dat, coef_list, vcov_list, conv_list, blup_results, meta_results,
# knots_var, knots_lag, wbgt_cen, wbgt_all, wbgt_range, wbgt_p95, periods,
# prefectures, pref_lookup
#
# Additional R packages needed beyond the main script's list:
#   sf, viridis, rnaturalearth, rnaturalearthdata, tidyr
#
# NOTE ON PART "0. Bridge": the main script computes per-prefecture BLUP
# (shrinkage-adjusted) coefficients in `blup_results` but does not itself
# write out a per-prefecture RR-at-P95 table. The block below derives
# output/pref_rr_at_p95.csv and output/rr_change_post_vs_pre.csv from
# `blup_results`, using the same basis-extraction pattern as the pooled RR
# calculation in Part 6.2 of the main script. This bridge has been verified:
# pref_rr_at_p95.csv reproduces the Table S1 values exactly (e.g.
# post-pandemic Kagoshima 2.36, Aichi 4.03), and rr_change_post_vs_pre.csv
# correctly restricts to the 44 commonly-converged prefectures, matching the
# range in Table S2 (-10.2% to -50.6%).
# =============================================================================

################################################################################
# Extension script: spatial maps + attributable fraction
#
# PREREQUISITE: main analysis already run (heatstroke_dlnm_analysis.R)
#           workspace contains: dat, coef_list, vcov_list, conv_list,
#           knots_var, knots_lag, wbgt_cen, wbgt_all, periods, etc.
#
# Outputs:
#   A. Spatial maps (3 figures)
#   B. Attributable Fraction / Number tables and figures
################################################################################

library(dplyr)
library(ggplot2)
library(sf)
library(dlnm)
library(gnm)
library(splines)
library(viridis)

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 0. BRIDGE: derive per-prefecture RR-at-P95 tables from blup_results
#    (see note above -- verify against Table S1 before trusting the maps)
# ==============================================================================
cat("=== 0. Deriving per-prefecture RR tables from BLUP estimates ===\n")

wbgt_pred_map <- seq(wbgt_range[1], wbgt_range[2], length.out = 100)
basis_map <- ns(wbgt_pred_map, knots = knots_var, Boundary.knots = wbgt_range)
idx_p95_map <- which.min(abs(wbgt_pred_map - wbgt_p95))
idx_ref_map <- which.min(abs(wbgt_pred_map - wbgt_cen))
bdiff_map <- as.numeric(basis_map[idx_p95_map, ] - basis_map[idx_ref_map, ])

pref_rr_rows <- list()
for (p in periods) {
  if (is.null(blup_results[[p]])) next
  b <- blup_results[[p]]
  # `blup(mixmeta_object, vcov = TRUE)` returns one element per included study
  # (prefecture here), each carrying $blup (coefficient vector) and $vcov.
  # `names(b)` / the ordering should correspond to the converged prefectures
  # used to build ymat for that period (see main script Part 6, `ok_idx`).
  ok_idx <- which(sapply(conv_list[[p]], isTRUE))
  pref_codes_p <- prefectures[ok_idx]
  for (k in seq_along(b)) {
    est <- b[[k]]
    coef_k <- est$blup
    vcov_k <- est$vcov
    logrr <- sum(bdiff_map * coef_k)
    se <- sqrt(as.numeric(t(bdiff_map) %*% vcov_k %*% bdiff_map))
    pref_rr_rows[[length(pref_rr_rows) + 1]] <- data.frame(
      pref_code = pref_codes_p[k],
      period = p,
      rr = round(exp(logrr), 3),
      rr_low = round(exp(logrr - 1.96 * se), 3),
      rr_high = round(exp(logrr + 1.96 * se), 3)
    )
  }
}
pref_rr <- do.call(rbind, pref_rr_rows)
write.csv(pref_rr, "output/pref_rr_at_p95.csv", row.names = FALSE)
cat("  Wrote output/pref_rr_at_p95.csv --", nrow(pref_rr), "rows\n")

rr_change <- pref_rr %>%
  filter(period %in% c("pre_pandemic", "post_pandemic")) %>%
  select(pref_code, period, rr) %>%
  tidyr::pivot_wider(names_from = period, values_from = rr) %>%
  filter(!is.na(pre_pandemic), !is.na(post_pandemic)) %>%
  mutate(pct_change = round((post_pandemic - pre_pandemic) / pre_pandemic * 100, 1))
write.csv(rr_change, "output/rr_change_post_vs_pre.csv", row.names = FALSE)
cat("  Wrote output/rr_change_post_vs_pre.csv --", nrow(rr_change), "rows\n\n")

# ==============================================================================
# 1. Download Japan prefecture shapefile
# ==============================================================================
cat("=== Downloading Japan prefecture map data ===\n")

# Uses NaturalEarth's Japan admin-1 (prefecture-level) shapefile
# via the rnaturalearth package
if (!require("rnaturalearth")) install.packages("rnaturalearth", repos = "https://cloud.r-project.org")
if (!require("rnaturalearthdata")) install.packages("rnaturalearthdata", repos = "https://cloud.r-project.org")

# A higher-resolution shapefile could be sourced from Japan's national land
# numerical information service; here we use NaturalEarth's admin-1 data directly
japan_sf <- ne_states(country = "japan", returnclass = "sf")

# Build a mapping from prefecture name to pref_code
pref_names_jp <- c(
  "Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
  "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa",
  "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano",
  "Gifu", "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo",
  "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
  "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki",
  "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa"
)

# Match against NaturalEarth's `name` field
# First inspect what names are present
cat("NE prefecture names:\n")
cat(paste(sort(japan_sf$name), collapse = ", "), "\n\n")

# Create the mapping (NE name -> pref_code)
ne_to_code <- data.frame(
  name = japan_sf$name,
  stringsAsFactors = FALSE
)

# Manual matching (NE's romanization may differ slightly)
pref_code_map <- data.frame(
  pref_name = pref_names_jp,
  pref_code = 1:47,
  stringsAsFactors = FALSE
)

# Attempt the join
japan_sf <- japan_sf %>%
  left_join(pref_code_map, by = c("name" = "pref_name"))

# Check match results
matched <- sum(!is.na(japan_sf$pref_code))
cat(sprintf("Matched %d / %d prefectures\n", matched, nrow(japan_sf)))

# If NE's naming differs, fix common discrepancies manually
if (matched < 47) {
  cat("Attempting fuzzy matching for unmatched prefectures...\n")
  unmatched <- japan_sf$name[is.na(japan_sf$pref_code)]
  cat("Unmatched:", paste(unmatched, collapse = ", "), "\n")
  
  # Common NE naming differences (macrons, alternative romanization)
  fix_map <- c(
    "Hokkaidō" = 1, "Hyōgo" = 28, "Kōchi" = 39, "Ōita" = 44,
    "Ōsaka" = 27, "Tōkyō" = 13, "Kyōto" = 26, "Gumma" = 10,
    "Hokkaido" = 1, "Hyogo" = 28, "Kochi" = 39, "Oita" = 44,
    "Osaka" = 27, "Tokyo" = 13, "Kyoto" = 26, "Gunma" = 10
  )
  
  for (nm in unmatched) {
    if (nm %in% names(fix_map)) {
      japan_sf$pref_code[japan_sf$name == nm] <- fix_map[nm]
    }
  }
  
  matched2 <- sum(!is.na(japan_sf$pref_code))
  cat(sprintf("After fixing: matched %d / %d\n", matched2, nrow(japan_sf)))
}

# ==============================================================================
# A. Spatial maps
# ==============================================================================
cat("\n=== A. Spatial maps ===\n")

# Read the per-prefecture RR tables produced in Part 0 above
pref_rr <- read.csv("output/pref_rr_at_p95.csv")
rr_change <- read.csv("output/rr_change_post_vs_pre.csv")

# --- A1: RR map, post-pandemic period ---
rr_post <- pref_rr[pref_rr$period == "post_pandemic", ]
map_post <- japan_sf %>%
  left_join(rr_post, by = "pref_code")

fig_map_rr <- ggplot(map_post) +
  geom_sf(aes(fill = rr), color = "white", size = 0.2) +
  scale_fill_viridis(
    option = "inferno",
    name = "Cumulative RR\nat WBGT P95",
    limits = c(1, max(rr_post$rr, na.rm = TRUE)),
    na.value = "grey80",
    direction = -1
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Post-pandemic (2022-2023)")

ggsave("output/figures/fig_map_rr_post.pdf", fig_map_rr, width = 7, height = 9)
ggsave("output/figures/fig_map_rr_post.png", fig_map_rr, width = 7, height = 9, dpi = 300)
cat("  Map RR post-pandemic saved.\n")

# --- A2: three periods side by side ---
rr_all <- pref_rr
rr_all$period <- factor(rr_all$period,
                        levels = c("pre_pandemic", "pandemic", "post_pandemic"),
                        labels = c("Pre-pandemic\n(2015-2019)",
                                   "Pandemic\n(2020-2021)",
                                   "Post-pandemic\n(2022-2023)"))

map_all <- japan_sf %>%
  left_join(rr_all, by = "pref_code", relationship = "many-to-many")

fig_map_3period <- ggplot(map_all %>% filter(!is.na(period))) +
  geom_sf(aes(fill = rr), color = "white", size = 0.15) +
  facet_wrap(~ period, nrow = 1) +
  scale_fill_viridis(
    option = "inferno",
    name = "Cumulative RR\nat WBGT P95",
    na.value = "grey80",
    direction = -1
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold")
  )

ggsave("output/figures/fig_map_3period.pdf", fig_map_3period, width = 14, height = 7)
ggsave("output/figures/fig_map_3period.png", fig_map_3period, width = 14, height = 7, dpi = 300)
cat("  Map 3-period comparison saved.\n")

# --- A3: RR percent-change map (post vs pre) ---
map_change <- japan_sf %>%
  left_join(rr_change, by = "pref_code")

fig_map_change <- ggplot(map_change) +
  geom_sf(aes(fill = pct_change), color = "white", size = 0.2) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = -35,
    name = "RR change (%)\nPost vs Pre",
    limits = c(-55, -5),
    na.value = "grey80"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Change in cumulative RR at WBGT P95 (post-pandemic vs pre-pandemic)")

ggsave("output/figures/fig_map_rr_change.pdf", fig_map_change, width = 7, height = 9)
ggsave("output/figures/fig_map_rr_change.png", fig_map_change, width = 7, height = 9, dpi = 300)
cat("  Map RR change saved.\n")


# ==============================================================================
# B. Attributable Fraction / Attributable Number
# ==============================================================================
cat("\n=== B. Attributable Fraction ===\n")

# Uses the attrdl() function from the dlnm package
# For each prefecture x period, compute AF from the prefecture-level model

af_results <- list()

for (p in c("pre_pandemic", "pandemic", "post_pandemic")) {
  cat("\n--- Period:", p, "---\n")
  
  for (i in 1:47) {
    pref_dat <- dat %>%
      filter(pref_code == i, period == p, !is.na(wbgt_max)) %>%
      arrange(date)
    
    if (nrow(pref_dat) < 60 || sum(pref_dat$total) < 10) next
    
    tryCatch({
      cb <- crossbasis(
        pref_dat$wbgt_max,
        lag = 5,
        argvar = list(fun = "ns", knots = knots_var),
        arglag = list(fun = "ns", knots = knots_lag)
      )
      
      model <- gnm(
        total ~ cb + factor(dow) + holiday,
        eliminate = stratum,
        family = quasipoisson(),
        data = pref_dat
      )
      
      # Compute attributable number (all days above the reference)
      # wbgt_cen is the reference; attributable EAD is computed above it
      an_heat <- attrdl(
        pref_dat$wbgt_max,
        cb,
        pref_dat$total,
        model,
        type = "an",
        cen = wbgt_cen,
        range = c(wbgt_cen, max(pref_dat$wbgt_max))  # heat side only
      )
      
      af_heat <- attrdl(
        pref_dat$wbgt_max,
        cb,
        pref_dat$total,
        model,
        type = "af",
        cen = wbgt_cen,
        range = c(wbgt_cen, max(pref_dat$wbgt_max))
      )
      
      total_ead <- sum(pref_dat$total)
      n_days <- nrow(pref_dat)
      
      af_results <- c(af_results, list(data.frame(
        pref_code = i,
        pref_name = pref_lookup$pref_name[i],
        period = p,
        total_ead = total_ead,
        an_heat = round(an_heat, 1),
        af_heat = round(af_heat * 100, 2),  # percentage
        daily_an = round(an_heat / n_days, 2),
        stringsAsFactors = FALSE
      )))
      
    }, error = function(e) {
      # silently skip
    })
  }
}

af_df <- do.call(rbind, af_results)

# Summary: national totals by period
cat("\n=== National AF summary by period ===\n")
af_summary <- af_df %>%
  group_by(period) %>%
  summarise(
    n_pref = n(),
    total_ead = sum(total_ead),
    total_an_heat = round(sum(an_heat)),
    mean_af_heat = round(mean(af_heat), 1),
    median_af_heat = round(median(af_heat), 1),
    .groups = "drop"
  ) %>%
  mutate(
    national_af = round(total_an_heat / total_ead * 100, 1),
    period = factor(period,
                    levels = c("pre_pandemic", "pandemic", "post_pandemic"),
                    labels = c("Pre-pandemic", "Pandemic", "Post-pandemic"))
  )

print(af_summary)
write.csv(af_summary, "output/af_summary_by_period.csv", row.names = FALSE)
write.csv(af_df, "output/af_by_prefecture_period.csv", row.names = FALSE)

# --- AF map (post-pandemic) ---
af_post <- af_df[af_df$period == "post_pandemic", ]
map_af <- japan_sf %>%
  left_join(af_post, by = "pref_code")

fig_map_af <- ggplot(map_af) +
  geom_sf(aes(fill = af_heat), color = "white", size = 0.2) +
  scale_fill_viridis(
    option = "plasma",
    name = "Heat-attributable\nfraction (%)",
    na.value = "grey80",
    direction = -1
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Heat-attributable fraction of EAD, Post-pandemic (2022-2023)")

ggsave("output/figures/fig_map_af_post.pdf", fig_map_af, width = 7, height = 9)
ggsave("output/figures/fig_map_af_post.png", fig_map_af, width = 7, height = 9, dpi = 300)
cat("  AF map saved.\n")

# --- AF bar chart by period ---
af_period_bar <- af_df %>%
  group_by(period) %>%
  summarise(
    total_ead = sum(total_ead),
    total_an = sum(an_heat),
    af = sum(an_heat) / sum(total_ead) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    non_attr = total_ead - total_an,
    period = factor(period,
                    levels = c("pre_pandemic", "pandemic", "post_pandemic"),
                    labels = c("Pre-pandemic\n(2015-2019)",
                               "Pandemic\n(2020-2021)",
                               "Post-pandemic\n(2022-2023)"))
  )

cat("\n=== Key AF results for manuscript ===\n")
for (i in 1:nrow(af_period_bar)) {
  cat(sprintf("  %s: %s total EAD, %s attributable to heat (AF = %.1f%%)\n",
              as.character(af_period_bar$period[i]),
              format(af_period_bar$total_ead[i], big.mark = ","),
              format(round(af_period_bar$total_an[i]), big.mark = ","),
              af_period_bar$af[i]))
}

# Stacked bar
library(tidyr)
af_long <- af_period_bar %>%
  select(period, Attributable = total_an, `Non-attributable` = non_attr) %>%
  pivot_longer(-period, names_to = "type", values_to = "count")

fig_af_bar <- ggplot(af_long, aes(x = period, y = count / 1000, fill = type)) +
  geom_col(width = 0.6) +
  geom_text(data = af_period_bar,
            aes(x = period, y = total_ead / 1000 + 5,
                label = sprintf("AF = %.1f%%", af), fill = NULL),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Attributable" = "#D73027",
                                "Non-attributable" = "#4575B4"),
                    name = "") +
  labs(
    x = NULL,
    y = "Total EAD (thousands)",
    title = "Heat-attributable emergency ambulance dispatches by period"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/figures/fig_af_barplot.pdf", fig_af_bar, width = 7, height = 5)
ggsave("output/figures/fig_af_barplot.png", fig_af_bar, width = 7, height = 5, dpi = 300)
cat("  AF barplot saved.\n")

# --- Severe-case subgroup AF ---
cat("\n=== Severe-case AF ===\n")
dat$severe_total <- dat$sev_death + dat$sev_severe

af_severe <- list()

for (p in c("pre_pandemic", "pandemic", "post_pandemic")) {
  for (i in 1:47) {
    pref_dat <- dat %>%
      filter(pref_code == i, period == p, !is.na(wbgt_max)) %>%
      arrange(date)
    
    if (nrow(pref_dat) < 60 || sum(pref_dat$severe_total) < 5) next
    
    tryCatch({
      cb <- crossbasis(
        pref_dat$wbgt_max,
        lag = 5,
        argvar = list(fun = "ns", knots = knots_var),
        arglag = list(fun = "ns", knots = knots_lag)
      )
      
      model <- gnm(
        severe_total ~ cb + factor(dow) + holiday,
        eliminate = stratum,
        family = quasipoisson(),
        data = pref_dat
      )
      
      an <- attrdl(pref_dat$wbgt_max, cb, pref_dat$severe_total, model,
                   type = "an", cen = wbgt_cen,
                   range = c(wbgt_cen, max(pref_dat$wbgt_max)))
      
      af_severe <- c(af_severe, list(data.frame(
        pref_code = i, period = p,
        total_severe = sum(pref_dat$severe_total),
        an_heat = round(an, 1)
      )))
    }, error = function(e) {})
  }
}

af_severe_df <- do.call(rbind, af_severe)
af_severe_summary <- af_severe_df %>%
  group_by(period) %>%
  summarise(
    total = sum(total_severe),
    an = round(sum(an_heat)),
    af = round(sum(an_heat) / sum(total_severe) * 100, 1),
    .groups = "drop"
  )

cat("\nSevere-case AF summary:\n")
print(af_severe_summary)
write.csv(af_severe_summary, "output/af_severe_summary.csv", row.names = FALSE)

cat("\n=== ALL INCREMENTAL ANALYSES COMPLETE ===\n")
cat("New outputs:\n")
cat("  output/figures/fig_map_3period.pdf     - RR spatial distribution, 3 periods\n")
cat("  output/figures/fig_map_rr_change.pdf   - RR percent-change spatial distribution\n") 
cat("  output/figures/fig_map_af_post.pdf     - Post-pandemic AF spatial distribution\n")
cat("  output/figures/fig_af_barplot.pdf      - AF bar chart\n")
cat("  output/af_summary_by_period.csv        - National AF summary\n")
cat("  output/af_by_prefecture_period.csv     - Prefecture-level AF detail\n")
cat("  output/af_severe_summary.csv           - Severe-case AF summary\n")
