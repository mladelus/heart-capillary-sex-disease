# =============================================================================
# 08  STAGE 3A (exploratory) - where in the heart do female and male cells differ?
# Discovery in one cohort, replication in the other, for each cell type (see ANALYSIS_PLAN.md, Stage 3).
# Per cohort and cell type: limma-voom per dataset/disease, inverse-variance meta of sex effects.
# Autosomal sex signal per cell type: (i) Spearman correlation of gene z-scores between cohorts;
# (ii) genes at FDR < 0.05 in one cohort replicated (p < 0.05, same direction) in the other.
# Pathways: Hallmark sets on combined (Stouffer) autosomal z-scores, rank-based gene-set test.
# Makes: stage3_de.rds, tables/S3A_celltype_summary.csv, tables/S3A_replicated_genes.csv, tables/S3A_hallmark.csv
# =============================================================================
C  <- load_cohorts()
ch <- chrom_map()
H  <- hallmark_sets()

res <- list(); summ <- list(); repl <- list(); hall <- list()
for (ct in names(CELLTYPES)) {
  mh <- C$healthy$mats[[CELLTYPES[[ct]][1]]]; md <- C$disease$mats[[CELLTYPES[[ct]][2]]]
  if (is.null(mh) || is.null(md)) { message(ct, ": missing in one cohort, skipped"); next }
  message("cell type: ", ct)
  dh <- de_meta(mh, C$healthy$samples); dd <- de_meta(md, C$disease$samples)
  if (!nrow(dh) || !nrow(dd)) { message(ct, ": too few donors, skipped"); next }
  nh <- attr(dh, "donors"); nd <- attr(dd, "donors")              # donors that entered each model
  j <- inner_join(dh, dd, by = "gene", suffix = c("_h", "_d")) |> left_join(ch, by = "gene") |>
    mutate(chr_class = ifelse(is.na(chr_class), "unknown", chr_class),
           z_comb = (z_h * sqrt(nrow(nh)) + z_d * sqrt(nrow(nd))) / sqrt(nrow(nh) + nrow(nd)),
           rep_h2d = FDR_h < 0.05 & p_d < 0.05 & sign(z_h) == sign(z_d),
           rep_d2h = FDR_d < 0.05 & p_h < 0.05 & sign(z_h) == sign(z_d))
  res[[ct]] <- j
  a <- j |> filter(chr_class == "autosome")
  ct_test <- cor.test(a$z_h, a$z_d, method = "spearman", exact = FALSE)
  summ[[ct]] <- tibble(cell_type = ct,
    healthy_F = sum(nh$sex == "female"), healthy_M = sum(nh$sex == "male"),
    disease_F = sum(nd$sex == "female"), disease_M = sum(nd$sex == "male"),
    autosomal_genes = nrow(a), rho_autosomal = unname(ct_test$estimate), rho_p = ct_test$p.value,
    FDR05_healthy = sum(a$FDR_h < 0.05), replicated_in_disease = sum(a$rep_h2d),
    FDR05_disease = sum(a$FDR_d < 0.05), replicated_in_healthy = sum(a$rep_d2h),
    replicated_either = sum(a$rep_h2d | a$rep_d2h),
    pct_rep_h2d = if (sum(a$FDR_h < 0.05) >= 10) 100 * sum(a$rep_h2d) / sum(a$FDR_h < 0.05) else NA,
    pct_rep_d2h = if (sum(a$FDR_d < 0.05) >= 10) 100 * sum(a$rep_d2h) / sum(a$FDR_d < 0.05) else NA,
    rho_logFC_disc_genes = if (sum(a$FDR_h < 0.05 | a$FDR_d < 0.05) >= 10)
      cor(a$logFC_h[a$FDR_h < 0.05 | a$FDR_d < 0.05], a$logFC_d[a$FDR_h < 0.05 | a$FDR_d < 0.05], method = "spearman") else NA,
    X_FDR05_both = sum(j$chr_class == "X" & j$FDR_h < 0.05 & j$FDR_d < 0.05 & sign(j$z_h) == sign(j$z_d)),
    Y_FDR05_both = sum(j$chr_class == "Y" & j$FDR_h < 0.05 & j$FDR_d < 0.05 & sign(j$z_h) == sign(j$z_d)))
  repl[[ct]] <- a |> filter(rep_h2d | rep_d2h) |> mutate(cell_type = ct) |>
    select(cell_type, gene, chr, logFC_h, FDR_h, logFC_d, FDR_d, z_comb) |> arrange(desc(abs(z_comb)))
  stat <- setNames(a$z_comb, a$gene)
  hall[[ct]] <- bind_rows(lapply(names(H), function(s) {
    idx <- which(names(stat) %in% H[[s]]); if (length(idx) < 10) return(NULL)
    tibble(cell_type = ct, set = sub("^HALLMARK_", "", s), genes = length(idx), mean_z = mean(stat[idx]),
           p_up = limma::geneSetTest(idx, stat, alternative = "up", ranks.only = TRUE),
           p_down = limma::geneSetTest(idx, stat, alternative = "down", ranks.only = TRUE))
  }))
  if (!nrow(hall[[ct]])) next
  hall[[ct]] <- hall[[ct]] |> mutate(direction = ifelse(mean_z > 0, "higher in females", "higher in males"),
                p = pmin(1, 2 * pmin(p_up, p_down)), FDR = p.adjust(p, "BH")) |> arrange(p)
}
summ <- bind_rows(summ) |> arrange(desc(rho_autosomal))

# ---- fairness check: same donors for every core cell type, and a sex-label permutation null for rho ----
# More donors = more power, so cell types are also compared on the donors that have all six core pseudobulks.
core <- c("capillary", "arterial", "venous", "mural", "fibroblast", "cardiomyocyte")
common <- lapply(c(healthy = 1, disease = 2), function(i) {
  co <- c("healthy", "disease")[i]
  Reduce(intersect, lapply(core, function(ct) colnames(C[[co]]$mats[[CELLTYPES[[ct]][i]]])))
})
perm_sex <- function(S) S |> group_by(stratum) |> mutate(sex = sample(sex)) |> ungroup()
NPERM <- 20
fair <- bind_rows(lapply(core, function(ct) {
  message("equal-donor check: ", ct)
  mh <- C$healthy$mats[[CELLTYPES[[ct]][1]]]; md <- C$disease$mats[[CELLTYPES[[ct]][2]]]
  Sh <- C$healthy$samples |> filter(sample %in% common$healthy); Sd <- C$disease$samples |> filter(sample %in% common$disease)
  rho_of <- function(Sh, Sd) {
    j <- inner_join(de_meta(mh, Sh), de_meta(md, Sd), by = "gene", suffix = c("_h", "_d")) |>
      inner_join(ch |> filter(chr_class == "autosome"), by = "gene")
    c(rho = cor(j$z_h, j$z_d, method = "spearman"), rep = sum((j$FDR_h < 0.05 & j$p_d < 0.05 | j$FDR_d < 0.05 & j$p_h < 0.05) & sign(j$z_h) == sign(j$z_d)))
  }
  obs <- rho_of(Sh, Sd)
  null <- sapply(seq_len(NPERM), function(i) rho_of(perm_sex(Sh), perm_sex(Sd)))
  tibble(cell_type = ct, donors_healthy = nrow(Sh), donors_disease = nrow(Sd),
         rho_equal_donors = obs[["rho"]], null_rho_max = max(abs(null["rho", ])),
         replicated_equal_donors = obs[["rep"]], null_replicated_max = max(null["rep", ]))
}))
cat("\n== Same donors for every core cell type (", length(common$healthy), " healthy, ", length(common$disease),
    " diseased), with ", NPERM, " sex-label permutations as a null ==\n", sep = "")
fair |> r3() |> print(width = Inf)
cat("\n== Where is the autosomal sex signal? (ranked by between-cohort agreement) ==\n")
summ |> r3() |> print(n = Inf, width = Inf)
cat("\nPositive controls, X and Y genes significant in both cohorts, by cell type: see X_FDR05_both / Y_FDR05_both above.\n")

repl <- bind_rows(repl)
cat("\n== Replicated autosomal genes (top 10 per cell type, by combined z) ==\n")
repl |> group_by(cell_type) |> slice_head(n = 10) |> ungroup() |> r3() |> print(n = Inf, width = Inf)

hall <- bind_rows(hall)
cat("\n== Hallmark pathways with FDR < 0.05 (combined autosomal sex z, per cell type) ==\n")
hall |> filter(FDR < 0.05) |> select(cell_type, set, genes, mean_z, direction, p, FDR) |> r3() |> print(n = Inf, width = Inf)
cat("\n== Hormone-response Hallmarks in every cell type ==\n")
hall |> filter(grepl("ESTROGEN|ANDROGEN", set)) |> select(cell_type, set, mean_z, direction, p, FDR) |>
  r3() |> print(n = Inf, width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(summ, "tables/S3A_celltype_summary.csv", row.names = FALSE)
write.csv(fair, "tables/S3A_equal_donor_check.csv", row.names = FALSE)
write.csv(repl, "tables/S3A_replicated_genes.csv", row.names = FALSE)
write.csv(hall, "tables/S3A_hallmark.csv", row.names = FALSE)
saveRDS(list(summary = summ, fair = fair, de = res, replicated = repl, hallmark = hall), "stage3_de.rds")
