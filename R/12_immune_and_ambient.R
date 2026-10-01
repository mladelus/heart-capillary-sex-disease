# =============================================================================
# 12  STAGE 3E (exploratory) - immune cells, estrogen response by cell type, and the ambient-RNA check
# See ANALYSIS_PLAN.md, Stage 3 E. Same donor-level methods as the primary analysis.
# Makes: tables/S3E_immune_programs.csv, S3E_composition.csv, S3E_immune_x_genes.csv, S3E_estrogen_oxphos.csv,
#        stage3_immune.rds, figures/S3_Fig_estrogen_by_celltype.*
# =============================================================================
C <- load_cohorts(); H <- hallmark_sets()
IFN <- c("IFIT1", "IFIT3", "ISG15", "MX1", "IFI6", "IFI44L", "OAS1", "STAT1", "IRF7", "XAF1")
IMMUNE <- list(
  myeloid = list(Interferon_I = IFN,
    MHC_II = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD74", "CIITA"),
    Inflammatory = c("IL1B", "TNF", "NLRP3", "CCL3", "CCL4", "CXCL8", "NFKBIA", "PTGS2"),
    Resident_macrophage = c("LYVE1", "F13A1", "MRC1", "CD163", "FOLR2", "STAB1"),
    Recruited_monocyte = c("CCR2", "PLAC8", "S100A8", "S100A9", "VCAN", "FCN1"),
    Complement_efferocytosis = c("C1QA", "C1QB", "C1QC", "MERTK", "TREM2", "APOE")),
  lymphoid = list(Cytotoxic = c("GZMB", "GZMK", "GZMA", "PRF1", "NKG7", "GNLY", "CCL5"),
    Interferon_I = IFN, B_cell = c("MS4A1", "CD79A", "CD79B", "BANK1")))

# donor scores for one cell type in both cohorts (TMM log-CPM and mean-z scores within each stratum)
score_ct <- function(ct, sets) bind_rows(lapply(1:2, function(i) {
  co <- c("healthy", "disease")[i]; m <- C[[co]]$mats[[CELLTYPES[[ct]][i]]]; if (is.null(m)) return(NULL)
  S <- core_cols(C[[co]]$samples) |> filter(sample %in% colnames(m))
  bind_rows(lapply(split(S, S$stratum), function(d) {
    if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX) return(NULL)
    y <- DGEList(m[, d$sample, drop = FALSE]); y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE])
    lc <- cpm(y, log = TRUE, prior.count = 1)
    bind_cols(d, as_tibble(sapply(sets, function(g) prog_score(lc, g))))
  }))
}))
meta_ct <- function(sc, vars, extra = NULL) bind_rows(lapply(vars, function(v) {
  per <- bind_rows(lapply(split(sc, sc$stratum), function(d) {
    if (sum(is.finite(d[[v]])) < 2 * MIN_PER_SEX) return(NULL)
    r <- sex_fit(d, v, extra); if (is.null(r)) NULL else mutate(r, stratum = d$stratum[1], cohort = d$cohort[1]) }))
  if (!nrow(per)) return(NULL)
  bind_rows(per |> group_by(cohort) |> group_modify(~ meta_one(.x)) |> ungroup(),
            meta_one(per) |> mutate(cohort = "both")) |> mutate(outcome = v)
}))

# ---- 1. immune programs ----
imm <- bind_rows(lapply(names(IMMUNE), function(ct) {
  sets <- c(IMMUNE[[ct]], pos_control, list(CM_ambient = CM_AMBIENT))
  sc <- score_ct(ct, sets)
  bind_rows(meta_ct(sc, names(c(IMMUNE[[ct]], pos_control))) |> mutate(adjusted = "none"),
            meta_ct(sc, names(IMMUNE[[ct]]), "CM_ambient") |> mutate(adjusted = "+ ambient heart-muscle RNA")) |>
    mutate(cell_type = ct)
})) |> group_by(cell_type, cohort, adjusted) |>
  mutate(p_holm = ifelse(outcome == "X_escape", NA, p.adjust(ifelse(outcome == "X_escape", NA, p), "holm"))) |> ungroup()
cat("\n== 1. Immune programs: female - male (SD units) ==\n")
imm |> select(cell_type, adjusted, outcome, cohort, k, n_w, n_m, est, lo, hi, p, p_holm, I2) |> r3() |>
  arrange(cell_type, adjusted, outcome, cohort) |> print(n = Inf, width = Inf)

# ---- 2. immune composition (cell counts per class are stored with each donor's download) ----
dirs <- c(healthy = PAPER1_DIR, disease = getwd())
ncl <- bind_rows(lapply(names(dirs), function(co) bind_rows(lapply(list.files(file.path(dirs[[co]], "dl"), full.names = TRUE), function(f) {
  x <- readRDS(f); n <- x$n_class
  tibble(sample = paste(substr(x$dataset_id, 1, 8), x$donor_id, sep = "|"), total = sum(n),
         myeloid = if ("myeloid" %in% names(n)) n[["myeloid"]] else 0L, lymphoid = if ("lymphoid" %in% names(n)) n[["lymphoid"]] else 0L)
}))))
comp_d <- bind_rows(core_cols(C$healthy$samples), core_cols(C$disease$samples)) |> inner_join(ncl, by = "sample") |>
  mutate(myeloid_share = qlogis((myeloid + 0.5) / (total + 1)), lymphoid_share = qlogis((lymphoid + 0.5) / (total + 1)))
comp <- meta_ct(comp_d, c("myeloid_share", "lymphoid_share"))
cat("\n== 2. Immune cell share of all cells: female - male (SD units of logit share) ==\n")
comp |> select(outcome, cohort, k, n_w, n_m, est, lo, hi, p, I2) |> r3() |> arrange(outcome, cohort) |> print(n = Inf, width = Inf)
cat("\nMedian % of all cells, by cohort and sex:\n")
comp_d |> group_by(cohort, sex) |> summarise(myeloid_pct = median(100 * myeloid / total), lymphoid_pct = median(100 * lymphoid / total),
                                             .groups = "drop") |> r3() |> print()

# ---- 3. immune X-linked genes reported to escape X inactivation (gene-level results of script 08) ----
XI <- c("TLR7", "TLR8", "CXCR3", "CD40LG", "BTK", "IL2RG", "CYBB", "KDM6A", "DDX3X")
de <- readRDS("stage3_de.rds")$de
xg <- bind_rows(lapply(c("myeloid", "lymphoid"), function(ct) if (is.null(de[[ct]])) NULL else
  de[[ct]] |> filter(gene %in% XI) |> transmute(cell_type = ct, gene, logFC_healthy = logFC_h, p_healthy = p_h,
                                               logFC_disease = logFC_d, p_disease = p_d, z_combined = z_comb,
                                               p_combined = 2 * pnorm(-abs(z_comb)))))
cat("\n== 3. Immune X-linked genes, female - male (log2); genes missing were too lowly expressed ==\n")
xg |> r3() |> arrange(cell_type, p_combined) |> print(n = Inf, width = Inf)

# ---- 4. estrogen response and oxidative phosphorylation, donor level, every cell type; ambient check ----
hs <- list(Estrogen_early = H$HALLMARK_ESTROGEN_RESPONSE_EARLY, OxPhos = H$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
           MYC_targets = H$HALLMARK_MYC_TARGETS_V1, CM_ambient = CM_AMBIENT, ESR1_gene = "ESR1")
eo <- bind_rows(lapply(names(CELLTYPES), function(ct) {
  sc <- score_ct(ct, hs[c("Estrogen_early", "OxPhos", "MYC_targets", "CM_ambient")]); if (is.null(sc) || !nrow(sc)) return(NULL)
  v <- c("Estrogen_early", "OxPhos", "MYC_targets")
  bind_rows(meta_ct(sc, c(v, "CM_ambient")) |> mutate(adjusted = "none"),
            if (ct != "cardiomyocyte") meta_ct(sc, v, "CM_ambient") |> mutate(adjusted = "+ ambient heart-muscle RNA"),
            meta_ct(sc |> filter(age < 50), "Estrogen_early") |> mutate(adjusted = "age < 50"),
            meta_ct(sc |> filter(age >= 50), "Estrogen_early") |> mutate(adjusted = "age >= 50")) |> mutate(cell_type = ct)
})) |> filter(cohort == "both")
esr1 <- readRDS("stage3_hormone.rds")$expr |> filter(gene == "ESR1") |> group_by(cell_type) |> summarise(ESR1_log2cpm = median(log2cpm))
cat("\n== 4. Donor-level scores by cell type, female - male (SD units), both cohorts combined ==\n")
eo |> left_join(esr1, by = "cell_type") |> select(outcome, adjusted, cell_type, ESR1_log2cpm, k, n_w, n_m, est, lo, hi, p) |>
  r3() |> arrange(outcome, adjusted, desc(ESR1_log2cpm)) |> print(n = Inf, width = Inf)
er <- eo |> filter(outcome == "Estrogen_early", adjusted == "none") |> left_join(esr1, by = "cell_type")
ct_cor <- cor.test(er$ESR1_log2cpm, er$est, method = "spearman", exact = FALSE)
cat(sprintf("\nAcross %d cell types: Spearman correlation between ESR1 expression and the female - male estrogen-response difference = %.2f (p = %.3f)\n",
            nrow(er), ct_cor$estimate, ct_cor$p.value))

dir.create("tables", showWarnings = FALSE)
write.csv(imm, "tables/S3E_immune_programs.csv", row.names = FALSE)
write.csv(comp, "tables/S3E_composition.csv", row.names = FALSE)
write.csv(xg, "tables/S3E_immune_x_genes.csv", row.names = FALSE)
write.csv(eo |> left_join(esr1, by = "cell_type"), "tables/S3E_estrogen_oxphos.csv", row.names = FALSE)
saveRDS(list(immune = imm, composition = comp, x_genes = xg, estrogen_oxphos = eo, esr1 = esr1), "stage3_immune.rds")

f <- ggplot(er, aes(ESR1_log2cpm, est, label = cell_type)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = "grey50") + geom_point(size = 2.2) +
  geom_text(vjust = -0.9, size = 2.7) + theme_bw(base_size = 9) +
  labs(x = "ESR1 expression in the cell type (median log2 CPM)", y = "Estrogen-response score, female minus male (SD units)")
ggsave("figures/S3_Fig_estrogen_by_celltype.pdf", f, width = 110, height = 85, units = "mm")
ggsave("figures/S3_Fig_estrogen_by_celltype.tiff", f, width = 110, height = 85, units = "mm", dpi = 300, compression = "lzw")
