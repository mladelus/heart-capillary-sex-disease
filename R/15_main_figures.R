# =============================================================================
# 15  MAIN FIGURES of the combined paper - multi-panel, colour (PDF + TIFF + PNG, 180 mm wide)
# Figure 2  Capillary function programs are alike in females and males (A forest, B heatmap, C donors, D control)
# Figure 3  X-Y paralog dosage (A X + Y by cell group, B X alone, C capillaries in two cohorts, D Y share in males)
# Figure 4  Where sex differences lie (A cell types, B cardiomyocytes, C capillaries, D GTEx, E top genes in 3 cohorts)
# Figure 5  GTEx bulk left ventricle (A donors, B programs by age, C endothelin-1 across datasets, D X + Y by age)
# Figure S  Hormone receptors, estrogen response, immune programs, capillary cell states
# Every panel fills its own cell (no shared margins), so there is no empty space between panels.
# Colours: sex = magenta / blue; female-minus-male effects = red (higher in females) to blue (higher in males)
# through white; datasets and diseases = colour-blind-safe palette (Okabe-Ito).
# Needs: combined.rds, p2_primary.rds, stage3_*.rds; paper 1 aim1.rds, aim4.rds, inv_datasets.rds
# =============================================================================
for (p in c("patchwork", "scales")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(patchwork)
FEM <- "#D81B60"; MAL <- "#1E88E5"; HI <- "#B2182B"; LO <- "#2166AC"
STRATA_COL <- c("Healthy 1" = "#4D4D4D", "Healthy 2" = "#A6A6A6", "DCM" = "#0072B2", "ACM" = "#CC79A7",
                "Infarction" = "#E69F00", "Myocarditis" = "#009E73")   # Healthy 1 = harmonized cohort, 2 = Heart Cell Atlas
CELL_COL <- c(capillary = "#D55E00", arterial = "#E69F00", venous = "#CC79A7", mural = "#009E73", pericyte = "#009E73",
              fibroblast = "#56B4E9", cardiomyocyte = "#0072B2", myeloid = "#8C6BB1", lymphoid = "#666666", "smooth muscle" = "#A6761D")
theme_set(theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        axis.line = element_line(colour = "grey30", linewidth = 0.3), axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(colour = "grey35", size = 8),
        plot.title.position = "plot", plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.01, 0.99),
        strip.text = element_text(face = "bold", size = 8), legend.key.size = unit(3.5, "mm"), legend.title = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0), plot.margin = margin(4, 6, 2, 12)))
# each panel is drawn in full inside its own cell, with its own letter
P <- function(p, tag) wrap_elements(full = p + labs(tag = tag))
save_fig <- function(p, name, h, w = 180) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path("figures", paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 300, compression = "lzw")
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 200)
}
div_fill <- function(name, lim) scale_fill_gradient2(low = LO, mid = "white", high = HI, midpoint = 0, limits = c(-lim, lim),
                                                     oob = scales::squish, name = name)
titles1 <- readRDS(file.path(PAPER1_DIR, "inv_datasets.rds")) |> select(dataset_id, dataset_title)
nice_stratum <- function(x) {
  x <- as.character(x); t1 <- titles1$dataset_title[match(x, titles1$dataset_id)]; x <- ifelse(is.na(t1), x, t1)
  case_when(grepl("^dilated", x) ~ "DCM", grepl("^arrhythmogenic", x) ~ "ACM",
            grepl("^myocardial infarction", x) ~ "Infarction", grepl("^myocarditis", x) ~ "Myocarditis",
            x == "Heart" ~ "Healthy 1", TRUE ~ "Healthy 2")
}
lab <- c(NO_eNOS = "Nitric oxide / eNOS", Endothelin_ACE = "Endothelin / ACE", Barrier = "Barrier / junction",
         FA_transport = "Fatty-acid transport", Angiogenic_tip = "Angiogenic tip",
         mural_per_capillary = "Pericytes or mural cells\nper capillary cell", X_escape = "X-escape\n(positive control)")

# =============================== FIGURE 2 ===============================
K  <- readRDS("combined.rds")
pd <- K$per |> filter(program %in% names(lab)) |>
  mutate(label = factor(lab[program], levels = rev(lab)), Stratum = factor(nice_stratum(stratum), levels = names(STRATA_COL)))
cd <- K$combined |> filter(program %in% names(lab)) |> mutate(label = factor(lab[program], levels = rev(lab)))
f2a <- ggplot() + annotate("rect", xmin = -MARGIN, xmax = MARGIN, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#7B3294") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(data = pd, aes(est, label, colour = Stratum, size = n_w + n_m), alpha = 0.85,
             position = position_jitter(width = 0, height = 0.22, seed = 3)) +
  geom_errorbar(data = cd, aes(xmin = lo, xmax = hi, y = label), width = 0, orientation = "y", linewidth = 0.8) +
  geom_point(data = cd, aes(est, label), shape = 23, size = 3, fill = "white", stroke = 0.9) +
  scale_colour_manual(values = STRATA_COL, name = NULL) + scale_size_area(max_size = 4, guide = "none") +
  labs(title = "All 138 hearts combined", subtitle = "Diamond = pooled estimate (95% CI); band = equivalence bounds",
       x = "Female minus male (SD units)", y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5), nrow = 1)) + theme(legend.position = "bottom")
f2b <- ggplot(pd, aes(Stratum, label, fill = est)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.1f", est)), size = 2.4) + div_fill("Female minus male (SD)", 2) +
  labs(title = "Every dataset and disease", subtitle = "Healthy 1 = harmonized cohort; 2 = Heart Cell Atlas", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", legend.key.width = unit(9, "mm"), legend.title = element_text(vjust = 0.9))

A1 <- readRDS(file.path(PAPER1_DIR, "aim1.rds")); A2 <- readRDS("p2_primary.rds")
progs <- intersect(names(lab), intersect(names(A1$scores), names(A2$scores)))
don <- bind_rows(A1$scores |> transmute(stratum = as.character(dataset_id), state = "Healthy", sex = as.character(sex), across(all_of(progs))),
                 A2$scores |> transmute(stratum = as.character(stratum), state = "Diseased", sex = as.character(sex), across(all_of(progs)))) |>
  pivot_longer(all_of(progs), names_to = "program", values_to = "score") |> filter(is.finite(score)) |>
  group_by(program, stratum) |> mutate(z = as.numeric(scale(score))) |> ungroup() |>
  mutate(label = factor(lab[program], levels = lab), state = factor(state, levels = c("Healthy", "Diseased")),
         Sex = factor(sex, levels = c("female", "male"), labels = c("Female", "Male")))
sex_violin <- function(d) ggplot(d, aes(state, z, fill = Sex, colour = Sex)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_violin(alpha = 0.25, linewidth = 0.3, position = position_dodge(width = 0.85), scale = "width") +
  geom_point(size = 0.6, alpha = 0.7, position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.85, seed = 1)) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5, linewidth = 0.3,
               position = position_dodge(width = 0.85), colour = "black") +
  scale_fill_manual(values = c(Female = FEM, Male = MAL), name = NULL) + scale_colour_manual(values = c(Female = FEM, Male = MAL), name = NULL) +
  facet_wrap(~ label, nrow = 1, labeller = label_wrap_gen(14)) + labs(x = NULL, y = "Program score (SD units)") +
  theme(legend.position = "bottom")
f2c <- sex_violin(don |> filter(program != "X_escape")) +
  labs(title = "Each dot is one donor: females and males overlap", subtitle = "Black bar = median; scores standardized within each dataset or disease")
f2d <- sex_violin(don |> filter(program == "X_escape")) + labs(title = "Positive control", subtitle = "A true sex difference") +
  theme(legend.position = "none")
fig2 <- (P(f2a, "A") | P(f2b, "B")) / ((P(f2c, "C") | P(f2d, "D")) + plot_layout(widths = c(3.6, 1))) + plot_layout(heights = c(1.1, 1))
save_fig(fig2, "Figure2_capillary_programs", 190)

# =============================== FIGURE 3 ===============================
A4 <- readRDS(file.path(PAPER1_DIR, "aim4.rds"))
hm <- function(meas, ttl) { d <- A4$by_group |> filter(measure == meas) |> mutate(star = ifelse(FDR < 0.05, "*", ""))
  ggplot(d, aes(group, pair, fill = est)) + geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = paste0(sprintf("%+.1f", est), star)), size = 2.3) + div_fill("Female minus male (log2)", 1.6) +
    labs(title = ttl, subtitle = "Healthy hearts, by cell type; * FDR < 0.05", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(), axis.line = element_blank(),
          legend.position = "bottom", legend.key.width = unit(9, "mm"), legend.title = element_text(vjust = 0.9)) }
f3a <- hm("X_plus_Y", "X + Y combined: less in females")
f3b <- hm("X_only", "X copy alone: more in females")
xy <- K$xy |> mutate(measure = recode(measure, X_only = "X copy alone", X_plus_Y = "X + Y combined"),
                     Cohort = factor(cohort, levels = c("healthy", "disease"), labels = c("Healthy (45 donors)", "Diseased (93 donors)")))
f3c <- ggplot(xy, aes(est, pair, colour = Cohort)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.6)) + facet_wrap(~ measure, scales = "free_x") +
  scale_colour_manual(values = c("Healthy (45 donors)" = "#4D4D4D", "Diseased (93 donors)" = "#E66101"), name = NULL) +
  labs(title = "Capillary endothelial cells: the same in two independent cohorts", x = "Female minus male (log2), 95% CI", y = NULL) +
  theme(legend.position = "bottom")
ys <- A4$pair_long |> filter(sex == "male", pair %in% c("KDM6A/UTY", "KDM5C/KDM5D", "USP9X/USP9Y", "DDX3X/DDX3Y"), is.finite(Y_share)) |>
  mutate(group = reorder(group, Y_share, FUN = median), endo = ifelse(group %in% c("capillary", "arterial"), "Endothelial", "Other heart cells"))
f3d <- ggplot(ys, aes(group, 100 * Y_share, fill = endo)) + geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.85) +
  facet_wrap(~ pair, nrow = 1) + scale_fill_manual(values = c(Endothelial = "#D55E00", "Other heart cells" = "#56B4E9"), name = NULL) +
  labs(title = "In males, endothelium draws less of each pair from the Y copy", x = NULL, y = "Y share of X + Y transcripts (%)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")
fig3 <- (P(f3a, "A") | P(f3b, "B")) / P(f3c, "C") / P(f3d, "D") + plot_layout(heights = c(1.15, 1, 1.05))
save_fig(fig3, "Figure3_XY_dosage", 225)

# =============================== FIGURE 4 ===============================
S3 <- readRDS("stage3_de.rds"); GX <- readRDS("stage3_gtex.rds")
Fr <- S3$fair |> mutate(cell_type = reorder(cell_type, replicated_equal_donors))
f4a <- ggplot(Fr, aes(replicated_equal_donors, cell_type, fill = cell_type)) + geom_col(width = 0.65, show.legend = FALSE) +
  geom_point(aes(null_replicated_max, cell_type), shape = 124, size = 6, colour = "black", show.legend = FALSE) +
  geom_text(aes(x = pmax(replicated_equal_donors, null_replicated_max), label = replicated_equal_donors), hjust = -0.9, size = 2.8) +
  scale_fill_manual(values = CELL_COL) + scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Only muscle cells stand out", subtitle = "Black bar = most seen with\nshuffled sex labels",
       x = "Non-X/Y genes differing by sex\nin both cohorts (same donors)", y = NULL)
zz <- function(ct, ttl) { d <- S3$de[[ct]] |> filter(chr_class == "autosome") |>
    mutate(rep = rep_h2d | rep_d2h, dir = factor(ifelse(z_comb > 0, "Higher in females", "Higher in males"),
                                                  levels = c("Higher in females", "Higher in males")))
  ggplot(d, aes(z_h, z_d)) + geom_hline(yintercept = 0, linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.25) +
    geom_point(data = d |> filter(!rep), colour = "grey80", size = 0.25, alpha = 0.5) +
    geom_point(data = d |> filter(rep), aes(colour = dir), size = 1.4) +
    geom_text(data = d |> filter(rep) |> arrange(desc(abs(z_comb))) |> head(8), aes(label = gene), size = 2.1, vjust = -0.7, check_overlap = TRUE) +
    scale_colour_manual(values = c("Higher in females" = FEM, "Higher in males" = MAL), name = NULL, drop = FALSE) +
    coord_cartesian(xlim = c(-9, 9), ylim = c(-9, 9)) +
    labs(title = ttl, subtitle = paste0(sum(d$rep), " of ", nrow(d), " genes replicated"),
         x = "Sex effect, healthy cohort (z)", y = "Sex effect, diseased cohort (z)") + theme(legend.position = "bottom") }
f4b <- zz("cardiomyocyte", "Cardiomyocytes")
f4c <- zz("capillary", "Capillary endothelium")
G <- GX$cm_genes |> filter(!is.na(logFC_gtex)) |>
  mutate(dir = ifelse(z_comb > 0, "Higher in females", "Higher in males"), Sig = ifelse(validated, "p < 0.05 in GTEx", "Not significant"))
f4d <- ggplot(G, aes(z_comb, logFC_gtex, colour = dir)) + geom_hline(yintercept = 0, linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_point(aes(shape = Sig), size = 1.8) +
  geom_text(data = G |> arrange(p_gtex) |> head(12), aes(label = gene), size = 2.1, vjust = -0.8, colour = "black", check_overlap = TRUE) +
  scale_colour_manual(values = c("Higher in females" = FEM, "Higher in males" = MAL), guide = "none") +
  scale_shape_manual(values = c("p < 0.05 in GTEx" = 16, "Not significant" = 1), name = NULL) +
  labs(title = "Confirmed in a third cohort (GTEx, 432 hearts)",
       subtitle = sprintf("%d of %d genes in the same direction", sum(G$same_direction), nrow(G)),
       x = "Single-cell cardiomyocytes (combined z)", y = "GTEx bulk left ventricle (log2)") + theme(legend.position = "bottom")
top <- G |> arrange(desc(abs(z_comb))) |> head(24) |>
  select(gene, z_comb, "Healthy\n(single-cell)" = logFC_healthy, "Diseased\n(single-cell)" = logFC_disease, "GTEx\n(bulk)" = logFC_gtex) |>
  pivot_longer(-c(gene, z_comb), names_to = "cohort", values_to = "logFC") |>
  mutate(gene = reorder(gene, z_comb), cohort = factor(cohort, levels = c("Healthy\n(single-cell)", "Diseased\n(single-cell)", "GTEx\n(bulk)")))
f4e <- ggplot(top, aes(cohort, gene, fill = logFC)) + geom_tile(colour = "white", linewidth = 0.5) + div_fill("Female minus male (log2)", 1.5) +
  labs(title = "Top genes in three cohorts", x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.text.y = element_text(size = 6.5),
        legend.position = "bottom", legend.key.width = unit(8, "mm"), legend.title = element_text(vjust = 0.9))
fig4 <- (P(f4a, "A") | P(f4b, "B") | P(f4c, "C")) / ((P(f4d, "D") | P(f4e, "E")) + plot_layout(widths = c(1.25, 1))) + plot_layout(heights = c(1, 1.2))
save_fig(fig4, "Figure4_where_sex_differences", 210)

# =============================== FIGURE 5 ===============================
f5a <- GX$meta |> count(AGE, sex) |> mutate(Sex = factor(sex, levels = c("female", "male"), labels = c("Female", "Male"))) |>
  ggplot(aes(AGE, n, fill = Sex)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.3) +
  scale_fill_manual(values = c(Female = FEM, Male = MAL), name = NULL) + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "GTEx has the younger women", subtitle = "Left ventricle, 432 donors", x = "Age (years)", y = "Donors") +
  theme(legend.position = "bottom")
gl <- c(cap_content = "Capillary content", NO_eNOS = "Nitric oxide / eNOS", Endothelin_ACE = "Endothelin / ACE", Prostacyclin = "Prostacyclin",
        Barrier = "Barrier / junction", FA_transport = "Fatty-acid transport", Angiogenic_tip = "Angiogenic tip")
ages <- c("50-79 years", "20-49 years"); AGE_COL <- c("20-49 years" = "#E66101", "50-79 years" = "#5E3C99")
band <- function(d, labs_vec) d |> filter(outcome %in% names(labs_vec), !grepl("^difference", contrast)) |>
  mutate(label = factor(labs_vec[outcome], levels = rev(labs_vec)), Age = factor(ifelse(grepl("20-49", contrast), ages[2], ages[1]), levels = ages))
band_plot <- function(d) ggplot(d, aes(est, label, colour = Age)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = AGE_COL, name = NULL, guide = guide_legend(reverse = TRUE)) + labs(y = NULL) + theme(legend.position = "bottom")
f5b <- band_plot(band(GX$age_bands, gl)) +
  labs(title = "Vessel-cell programs: higher in women only before 50?", subtitle = "Suggestive; not significant after correction",
       x = "Female minus male (SD units), adjusted for capillary content")
ed <- bind_rows(
  bind_rows(lapply(c("capillary", "venous", "arterial"), function(ct) { r <- S3$de[[ct]] |> filter(gene == "EDN1"); if (!nrow(r)) return(NULL)
    bind_rows(tibble(what = paste0(ct, ", healthy"), est = r$logFC_h, se = r$se_h, kind = "Single-cell, healthy"),
              tibble(what = paste0(ct, ", diseased"), est = r$logFC_d, se = r$se_d, kind = "Single-cell, diseased")) })),
  GX$endothelin |> filter(gene == "EDN1") |> transmute(what = "GTEx whole tissue", est = logFC_gtex, se = abs(logFC_gtex / t_gtex), kind = "GTEx bulk (432 hearts)")) |>
  mutate(what = factor(what, levels = rev(what)))
f5c <- ggplot(ed, aes(est, what, colour = kind)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), width = 0, orientation = "y", linewidth = 0.6) + geom_point(size = 2.2) +
  scale_colour_manual(values = c("Single-cell, healthy" = "#4D4D4D", "Single-cell, diseased" = "#E66101", "GTEx bulk (432 hearts)" = "#1B9E77"), name = NULL) +
  labs(title = "Endothelin-1 (EDN1): not confirmed in bulk tissue", x = "Female minus male (log2), 95% CI", y = NULL) +
  theme(legend.position = "bottom")
xl <- c(XplusY_KDM5C_KDM5D = "KDM5C + KDM5D", XplusY_KDM6A_UTY = "KDM6A + UTY", XplusY_USP9X_USP9Y = "USP9X + USP9Y",
        XplusY_DDX3X_DDX3Y = "DDX3X + DDX3Y", X_escape = "X-escape program")
f5d <- band_plot(band(GX$age_bands, xl)) + labs(title = "Sex-chromosome dosage: the same at every age", x = "Female minus male (SD units)")
fig5 <- ((P(f5a, "A") | P(f5b, "B")) + plot_layout(widths = c(1, 1.3))) / (P(f5c, "C") | P(f5d, "D"))
save_fig(fig5, "Figure5_GTEx", 170)

# =============================== FIGURE S (stage 3 nulls) ===============================
Hm <- readRDS("stage3_hormone.rds")$expr |> group_by(gene, cell_type) |> summarise(v = median(log2cpm), .groups = "drop")
fsa <- ggplot(Hm, aes(cell_type, gene, fill = v)) + geom_tile(colour = "white", linewidth = 0.6) + geom_text(aes(label = sprintf("%.1f", v)), size = 2.4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", name = "Median log2 CPM") +
  labs(title = "Sex-hormone receptors by cell type", subtitle = "ESR2 is likely inflated by reads from a neighbouring gene", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", legend.key.width = unit(8, "mm"), legend.title = element_text(vjust = 0.9))
Im <- readRDS("stage3_immune.rds")
er <- Im$estrogen_oxphos |> filter(outcome == "Estrogen_early", adjusted == "none") |> left_join(Im$esr1, by = "cell_type")
fsb <- ggplot(er, aes(ESR1_log2cpm, est, colour = cell_type)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5) + geom_point(size = 2.4) +
  scale_colour_manual(values = CELL_COL, name = NULL) +
  labs(title = "Estrogen response does not track receptor level", x = "ESR1 expression (median log2 CPM)", y = "Estrogen response,\nfemale minus male (SD)") +
  guides(colour = guide_legend(nrow = 2)) + theme(legend.position = "bottom")
im <- Im$immune |> filter(cohort == "both", adjusted == "none", outcome != "X_escape") |> mutate(what = paste0(cell_type, ": ", gsub("_", " ", outcome)))
fsc <- ggplot(im, aes(est, reorder(what, est), colour = cell_type)) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.6) + geom_point(size = 2) +
  scale_colour_manual(values = CELL_COL, guide = "none") + labs(title = "Immune programs: no sex difference", x = "Female minus male (SD units), 138 donors", y = NULL)
St <- readRDS("stage3_states_var.rds")$states |> filter(cohort == "both", adjusted == "none") |> mutate(outcome = sub("^hi_", "", outcome))
fsd <- ggplot(St, aes(est, reorder(outcome, est))) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.6, colour = "#D55E00") + geom_point(size = 2, colour = "#D55E00") +
  labs(title = "Capillary cell states: no sex difference", x = "Female minus male (SD units), 138 donors", y = NULL)
figS <- (P(fsa, "A") | P(fsb, "B")) / ((P(fsc, "C") | P(fsd, "D")) + plot_layout(widths = c(1.25, 1)))
save_fig(figS, "FigureS_hormones_immune_states", 175)

# ---- protein network of the cardiomyocyte genes, without added partners ----
E <- readRDS("stage3_cm_pathways.rds")$string_edges
if (!is.null(E) && nrow(E) && requireNamespace("igraph", quietly = TRUE)) {
  g <- igraph::graph_from_data_frame(E[, c("preferredName_A", "preferredName_B", "score")], directed = FALSE)
  z <- setNames(GX$cm_genes$z_comb, GX$cm_genes$gene)[igraph::V(g)$name]
  pdf("figures/FigureS_cm_network.pdf", width = 5, height = 5); par(mar = c(0, 0, 0, 0)); set.seed(1)
  plot(g, vertex.color = ifelse(z > 0, FEM, MAL), vertex.size = 26, vertex.label.cex = 0.6, vertex.label.color = "white",
       vertex.frame.color = NA, edge.width = 1 + 3 * igraph::E(g)$score, edge.color = "grey60")
  legend("bottomleft", c("higher in females", "higher in males"), pt.bg = c(FEM, MAL), pch = 21, bty = "n", cex = 0.7)
  dev.off()
}
message("Figures saved in ", file.path(getwd(), "figures"), " (PDF, TIFF and PNG)")
