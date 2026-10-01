# =============================================================================
# 11  COMBINED HEALTH + DISEASE SUMMARIES AND STAGE 3 FIGURES
# Female - male per stratum from both cohorts (2 healthy datasets + 4 disease strata), random-effects
# meta-analysis with disease state as a moderator, for the six programs, X-escape, the mural/pericyte-to-capillary
# ratio and X-Y combined dosage in capillaries. Then figures for stage 3.
# Needs: paper 1 aim1.rds, aim4.rds, pb.rds; this folder p2_primary.rds, p2_pericytes_xy.rds, stage3_*.rds
# Makes: tables/S3E_combined.csv, figures/S3_*.pdf/.tiff, combined.rds
# =============================================================================
A1 <- readRDS(file.path(PAPER1_DIR, "aim1.rds")); A2 <- readRDS("p2_primary.rds")
per <- bind_rows(A1$per_ds |> transmute(program, stratum = dataset_id, est, se, n_w, n_m, state = "healthy"),
                 A2$per_st |> transmute(program, stratum, est, se, n_w, n_m, state = "disease")) |>
  filter(program %in% c(PRIMARY, "X_escape"))

# mural (disease) / pericyte (healthy) per capillary EC, per stratum
C <- load_cohorts()
# note: healthy = cells labeled pericyte; disease = mural cells (pericyte + smooth muscle labels), see Deviation 1
ratio <- bind_rows(core_cols(C$healthy$samples), core_cols(C$disease$samples) |> filter(n_EC >= MIN_EC)) |>
  mutate(y = log((n_ratio_cells + 0.5) / (n_capillary + 0.5)))
per <- bind_rows(per, bind_rows(lapply(split(ratio, ratio$stratum), function(d) {
  r <- sex_fit(d, "y"); if (is.null(r)) NULL else mutate(r, program = "mural_per_capillary", stratum = d$stratum[1], state = d$cohort[1]) })))

comb <- per |> group_by(program) |> group_modify(~ {
  d <- .x |> filter(is.finite(est), is.finite(se)); if (nrow(d) < 2) return(tibble())
  m  <- tryCatch(rma(yi = d$est, sei = d$se, method = "REML"), error = function(e) rma(yi = d$est, sei = d$se, method = "DL"))
  mm <- if (n_distinct(d$state) == 2) tryCatch(rma(yi = d$est, sei = d$se, mods = ~ state, data = d, method = "REML"), error = function(e) NULL) else NULL
  tibble(k = nrow(d), n_w = sum(d$n_w), n_m = sum(d$n_m), est = m$b[1], lo = m$ci.lb, hi = m$ci.ub, p = m$pval, I2 = m$I2,
         healthy_vs_disease_p = if (is.null(mm)) NA else mm$pval[2])
}) |> ungroup() |>
  mutate(healthy_vs_disease_p = ifelse(program == "mural_per_capillary", NA, healthy_vs_disease_p),  # definitions differ
         p_holm = ifelse(program %in% PRIMARY, p.adjust(ifelse(program %in% PRIMARY, p, NA), "holm"), NA))
cat("\n== Female - male, both cohorts combined (SD units; 2 healthy + 4 disease strata) ==\n")
comb |> r3() |> print(width = Inf)

# X-Y combined dosage in capillaries: pooled healthy and disease meta estimates (inverse variance)
xh <- readRDS(file.path(PAPER1_DIR, "aim4.rds"))$by_group |> filter(group == "capillary") |> transmute(pair, measure, est, lo, hi, cohort = "healthy")
xd <- readRDS("p2_pericytes_xy.rds")$xy_res |> transmute(pair, measure, est, lo, hi, cohort = "disease")
xy <- bind_rows(xh, xd) |> mutate(se = (hi - lo) / (2 * qnorm(0.975)))
xy_comb <- xy |> group_by(pair, measure) |> filter(n() == 2) |>
  summarise(est_comb = sum(est / se^2) / sum(1 / se^2), se_comb = sqrt(1 / sum(1 / se^2)),
            het_p = pchisq(sum(((est - sum(est / se^2) / sum(1 / se^2)) / se)^2), 1, lower.tail = FALSE), .groups = "drop") |>
  mutate(lo = est_comb - 1.96 * se_comb, hi = est_comb + 1.96 * se_comb, p = 2 * pnorm(-abs(est_comb / se_comb)))
cat("\n== X-Y paralogs in capillary ECs, female - male (log2), healthy and disease combined ==\n")
xy_comb |> r3() |> arrange(measure, pair) |> print(n = Inf, width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(bind_rows(comb, xy_comb |> rename(est = est_comb)), "tables/S3E_combined.csv", row.names = FALSE)
saveRDS(list(per = per, combined = comb, xy = xy, xy_comb = xy_comb), "combined.rds")

# ---- figures ----
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank()))
save_fig <- function(p, name, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = 170, height = h, units = "mm")
  ggsave(file.path("figures", paste0(name, ".tiff")), p, width = 170, height = h, units = "mm", dpi = 300, compression = "lzw")
}
lab <- c(NO_eNOS = "NO / eNOS", Endothelin_ACE = "Endothelin / ACE", Prostacyclin = "Prostacyclin", Barrier = "Barrier / junction",
         FA_transport = "Fatty-acid transport", Angiogenic_tip = "Angiogenic tip", mural_per_capillary = "Pericytes (healthy) / mural cells (disease) per capillary EC",
         X_escape = "X-escape (positive control)")
pd <- per |> filter(program %in% names(lab)) |> mutate(label = factor(lab[program], levels = rev(lab)))
cd <- comb |> filter(program %in% names(lab)) |> mutate(label = factor(lab[program], levels = rev(lab)))
f1 <- ggplot() + annotate("rect", xmin = -MARGIN, xmax = MARGIN, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(data = pd, aes(est, label, colour = state, size = n_w + n_m), shape = 1,
             position = position_jitter(width = 0, height = 0.18, seed = 1)) +
  geom_errorbar(data = cd, aes(xmin = lo, xmax = hi, y = label), width = 0, orientation = "y", linewidth = 0.6) +
  geom_point(data = cd, aes(est, label), shape = 18, size = 3.2) +
  scale_colour_manual(values = c(healthy = "grey40", disease = "#B2182B"), name = "Stratum") +
  scale_size_area(max_size = 3, name = "Donors") +
  labs(x = "Female minus male (donor-level SD units); diamonds = all 138 donors combined", y = NULL)
save_fig(f1, "S3_Fig_combined_health_disease", 95)

S3 <- readRDS("stage3_de.rds")$summary
f2 <- S3 |> mutate(cell_type = reorder(cell_type, rho_autosomal)) |>
  ggplot(aes(rho_autosomal, cell_type)) + geom_vline(xintercept = 0, linewidth = 0.3) + geom_col(fill = "grey55", width = 0.6) +
  geom_text(aes(label = paste0(replicated_either, " genes")), hjust = -0.15, size = 2.6) +
  labs(x = "Agreement of autosomal sex effects between healthy and diseased cohorts (Spearman rho)\nlabel: autosomal genes replicated across cohorts", y = NULL) +
  expand_limits(x = max(S3$rho_autosomal, na.rm = TRUE) * 1.35)
save_fig(f2, "S3_Fig_where_sex_differences", 70)

Hm <- readRDS("stage3_hormone.rds")$expr |> group_by(gene, cell_type) |> summarise(v = median(log2cpm), .groups = "drop")
f3 <- ggplot(Hm, aes(cell_type, gene, fill = v)) + geom_tile() + geom_text(aes(label = sprintf("%.1f", v)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Median\nlog2 CPM") + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_fig(f3, "S3_Fig_hormone_receptors", 60)

St <- readRDS("stage3_states_var.rds")$states |> filter(cohort == "both") |>
  mutate(outcome = sub("^hi_", "", outcome))
f4 <- ggplot(St, aes(est, reorder(outcome, est))) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y") + geom_point(shape = 18, size = 3) +
  labs(x = "Female minus male (SD units), capillary cell states, both cohorts combined", y = NULL)
save_fig(f4, "S3_Fig_capillary_states", 55)
message("Stage 3 figures in ", file.path(getwd(), "figures"))
