# =============================================================================
# 09  STAGE 3B (exploratory) - sex-hormone receptor map
# ESR1, ESR2, GPER1, AR, PGR, CYP19A1: log2 CPM by cell type, sex and cohort; in capillaries, female - male by
# age band; Hallmark estrogen- and androgen-response scores in capillaries, female - male, both cohorts combined.
# Makes: tables/S3B_receptor_expression.csv, tables/S3B_capillary_hormone.csv, stage3_hormone.rds
# =============================================================================
C <- load_cohorts()
H <- hallmark_sets()
REC <- c("ESR1", "ESR2", "GPER1", "AR", "PGR", "CYP19A1")

lc_all <- function(mat) {                        # TMM log2 CPM over all donors of one cohort and cell type
  y <- calcNormFactors(DGEList(mat)); cpm(y, log = TRUE, prior.count = 1)
}
expr <- bind_rows(lapply(names(CELLTYPES), function(ct) bind_rows(lapply(c("healthy", "disease"), function(co) {
  m <- C[[co]]$mats[[CELLTYPES[[ct]][if (co == "healthy") 1 else 2]]]; if (is.null(m)) return(NULL)
  S <- C[[co]]$samples |> filter(sample %in% colnames(m)); m <- m[, S$sample, drop = FALSE]
  lc <- lc_all(m); cp <- cpm(calcNormFactors(DGEList(m)))
  bind_rows(lapply(intersect(REC, rownames(lc)), function(g) tibble(cell_type = ct, cohort = co, gene = g,
    sex = S$sex, log2cpm = lc[g, ], detected = cp[g, ] > 1)))
}))))
tab <- expr |> group_by(cell_type, cohort, gene, sex) |>
  summarise(donors = n(), median_log2cpm = median(log2cpm), pct_donors_cpm_gt1 = 100 * mean(detected), .groups = "drop")
cat("\n== Sex-hormone receptors: median log2 CPM by cell type (both sexes, both cohorts) ==\n")
expr |> group_by(gene, cell_type) |> summarise(median = median(log2cpm), .groups = "drop") |>
  pivot_wider(names_from = cell_type, values_from = median) |> r3() |> print(width = Inf)
cat("\n== Percent of donor pseudobulks with CPM > 1 ==\n")
expr |> group_by(gene, cell_type) |> summarise(pct = round(100 * mean(detected)), .groups = "drop") |>
  pivot_wider(names_from = cell_type, values_from = pct) |> print(width = Inf)

# capillaries: female - male per receptor and hormone-response score, per stratum, by age band; meta across all strata
cap_scores <- bind_rows(lapply(c("healthy", "disease"), function(co) {
  m <- C[[co]]$mats$capillary; S <- core_cols(C[[co]]$samples) |> filter(sample %in% colnames(m), n_capillary >= MIN_CAP)
  bind_rows(lapply(split(S, S$stratum), function(d) {
    y <- DGEList(m[, d$sample, drop = FALSE]); y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE])
    lc <- cpm(y, log = TRUE, prior.count = 1); lca <- lc_all(m[, d$sample, drop = FALSE])
    out <- d |> select(sample, stratum, sex, age, assay) |> mutate(cohort = co)
    for (g in REC) out[[g]] <- if (g %in% rownames(lca)) lca[g, d$sample] else NA_real_
    for (s in c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_ANDROGEN_RESPONSE"))
      out[[sub("HALLMARK_", "", s)]] <- prog_score(lc, H[[s]])
    out
  }))
}))
outs <- c(REC, "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "ANDROGEN_RESPONSE")
fit_band <- function(dat, label) bind_rows(lapply(outs, function(v) {
  per <- bind_rows(lapply(split(dat, dat$stratum), function(d) {
    if (sum(is.finite(d[[v]])) < 2 * MIN_PER_SEX || isTRUE(sd(d[[v]], na.rm = TRUE) == 0)) return(NULL)
    r <- sex_fit(d, v); if (is.null(r)) NULL else mutate(r, stratum = d$stratum[1], cohort = d$cohort[1])
  }))
  if (!nrow(per)) return(NULL)
  bind_rows(per |> group_by(cohort) |> group_modify(~ meta_one(.x)) |> ungroup(),
            meta_one(per) |> mutate(cohort = "both")) |> mutate(outcome = v, ages = label)
}))
hb <- bind_rows(fit_band(cap_scores, "all ages"),
                fit_band(cap_scores |> filter(age < 50), "< 50 years"),
                fit_band(cap_scores |> filter(age >= 50), ">= 50 years"))
cat("\n== Capillary ECs: female - male (SD units), receptors and hormone-response scores ==\n")
hb |> select(ages, cohort, outcome, k, n_w, n_m, est, lo, hi, p) |> r3() |> arrange(ages, outcome, cohort) |>
  print(n = Inf, width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(tab, "tables/S3B_receptor_expression.csv", row.names = FALSE)
write.csv(hb, "tables/S3B_capillary_hormone.csv", row.names = FALSE)
saveRDS(list(expr = expr, table = tab, capillary = hb, scores = cap_scores), "stage3_hormone.rds")
