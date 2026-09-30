# =============================================================================
# 07  FIGURES (PDF + 300 dpi TIFF, 170 mm wide)
# Figure 2  primary: female - male in disease, pooled and per stratum, with equivalence bounds
# Figure 3  sex x DCM (same consortium): female - male in health and in DCM, with interaction
# Figure 4  pericytes and X-Y dosage in disease vs health
# =============================================================================
A <- readRDS("p2_primary.rds"); B <- readRDS("p2_interaction.rds"); C <- readRDS("p2_pericytes_xy.rds")
dir.create("figures", showWarnings = FALSE)
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank()))
save_fig <- function(p, name, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = 170, height = h, units = "mm")
  ggsave(file.path("figures", paste0(name, ".tiff")), p, width = 170, height = h, units = "mm", dpi = 300, compression = "lzw")
}
lab <- c(NO_eNOS = "NO / eNOS", Endothelin_ACE = "Endothelin / ACE", Prostacyclin = "Prostacyclin",
         Barrier = "Barrier / junction", FA_transport = "Fatty-acid transport", Angiogenic_tip = "Angiogenic tip",
         X_escape = "X-escape (positive control)")
keep <- intersect(names(lab), A$main$program)

m <- A$main |> filter(program %in% keep) |> mutate(label = factor(lab[program], levels = rev(lab[keep])))
d <- A$per_st |> filter(program %in% keep) |> mutate(label = factor(lab[program], levels = levels(m$label)))
f2 <- ggplot() +
  annotate("rect", xmin = -MARGIN, xmax = MARGIN, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(data = d, aes(est, label, size = n_w + n_m, colour = stratum), shape = 1,
             position = position_jitter(width = 0, height = 0.15, seed = 1)) +
  geom_errorbar(data = m, aes(xmin = lo, xmax = hi, y = label), width = 0, orientation = "y", linewidth = 0.6) +
  geom_point(data = m, aes(est, label), shape = 18, size = 3.2) +
  scale_size_area(max_size = 3, name = "Donors") +
  labs(x = "Female minus male in diseased hearts (donor-level SD units)", y = NULL, colour = "Stratum")
save_fig(f2, "P2_Figure2_primary", 100)

i <- B$inter |> filter(program %in% keep) |>
  select(program, healthy = sex_in_healthy, DCM = sex_in_DCM) |>
  pivot_longer(c(healthy, DCM), names_to = "state", values_to = "est") |>
  mutate(label = factor(lab[program], levels = rev(lab[keep])), state = factor(state, levels = c("healthy", "DCM")))
f3 <- ggplot(i, aes(est, label, colour = state)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_line(aes(group = label), colour = "grey70") + geom_point(size = 2.5) +
  scale_colour_manual(values = c(healthy = "grey45", DCM = "#B2182B"), name = NULL) +
  labs(x = "Female minus male (SD units), same-consortium model", y = NULL)
save_fig(f3, "P2_Figure3_sex_by_DCM", 80)

xy <- C$xy_res |> filter(pair %in% c("KDM6A/UTY", "KDM5C/KDM5D", "USP9X/USP9Y", "DDX3X/DDX3Y"))
if ("healthy" %in% names(xy)) {
  xy2 <- bind_rows(xy |> transmute(pair, measure, state = "disease", est),
                   xy |> transmute(pair, measure, state = "healthy (paper 1)", est = healthy))
  f4 <- ggplot(xy2, aes(est, pair, colour = state)) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(size = 2.5, position = position_dodge(width = 0.4)) + facet_wrap(~ measure) +
    scale_colour_manual(values = c(disease = "#B2182B", "healthy (paper 1)" = "grey45"), name = NULL) +
    labs(x = "Female minus male in capillary ECs (log2)", y = NULL)
  save_fig(f4, "P2_Figure4_XY_dosage", 70)
}
writeLines(capture.output(sessionInfo()), "sessionInfo_paper2.txt")
message("Figures in ", file.path(getwd(), "figures"))
