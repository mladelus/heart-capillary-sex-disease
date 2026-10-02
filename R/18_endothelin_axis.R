# =============================================================================
# 18  STAGE 3J (exploratory candidate-gene analysis) - the endothelin axis
# Endothelin-1 (EDN1) is made by endothelium; it acts on EDNRA (mural cells, cardiomyocytes, fibroblasts) and EDNRB (endothelium).
# (a) gene-level female - male (log2) for ligand, enzyme and receptors, per dataset/disease and pooled
# (b) percentage of capillary ECs with any EDN1 transcript per donor (cell level)
# (c) capillary-to-mural signalling score (EDN1 + EDNRA)
# (d) what tracks capillary EDN1 across donors, and whether it explains the sex difference
# (e) EDN1 sex difference below and above age 50
# Makes: tables/S3J_*.csv, stage3_endothelin.rds, figures/S3_Fig_endothelin.*
# =============================================================================
C <- load_cohorts(); H <- hallmark_sets()
ct_index <- function(ct, co) CELLTYPES[[ct]][if (co == "healthy") 1 else 2]
# log2 CPM of chosen genes for every donor of one cell type, both cohorts (TMM within dataset/disease, no gene filter)
gene_expr <- function(ct, genes) bind_rows(lapply(c("healthy", "disease"), function(co) {
  m <- C[[co]]$mats[[ct_index(ct, co)]]; if (is.null(m)) return(NULL)
  S <- core_cols(C[[co]]$samples) |> filter(sample %in% colnames(m))
  bind_rows(lapply(split(S, S$stratum), function(d) {
    if (nrow(d) < 2) return(NULL)
    lc <- cpm(calcNormFactors(DGEList(m[, d$sample, drop = FALSE])), log = TRUE, prior.count = 1)
    for (g in genes) { d[[g]] <- if (g %in% rownames(lc)) lc[g, d$sample] else NA_real_
      d[[paste0(g, "__cnt")]] <- if (g %in% rownames(m)) m[g, d$sample] else 0 }          # raw counts, to drop genes that are not expressed
    d
  }))
}))
lfit <- function(d, v, extra = NULL) {             # female - male in the units of v (not standardized)
  d <- droplevels(d[is.finite(d[[v]]) & !is.na(d$age), ])
  if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX || sd(d[[v]]) == 0) return(NULL)
  cn <- paste0(v, "__cnt"); if (cn %in% names(d) && mean(d[[cn]] > 0) < 0.5) return(NULL)    # detected in fewer than half the donors
  rhs <- c("sex", "age", extra, if (n_distinct(d$assay) > 1) "assay")
  co <- tryCatch(summary(lm(as.formula(paste0("`", v, "` ~ ", paste(rhs, collapse = " + "))), data = d))$coefficients, error = function(e) NULL)
  if (is.null(co) || !"sexfemale" %in% rownames(co)) return(NULL)
  tibble(est = co["sexfemale", 1], se = co["sexfemale", 2], p = co["sexfemale", 4], n_w = sum(d$sex == "female"), n_m = sum(d$sex == "male"))
}
pool <- function(dat, v, extra = NULL, fit = lfit) {
  per <- bind_rows(lapply(split(dat, dat$stratum), function(d) { r <- fit(d, v, extra)
    if (is.null(r)) NULL else mutate(r, stratum = d$stratum[1], cohort = d$cohort[1]) }))
  if (!nrow(per)) return(list(per = tibble(), meta = tibble()))
  list(per = per, meta = bind_rows(per |> group_by(cohort) |> group_modify(~ meta_one(.x)) |> ungroup(), meta_one(per) |> mutate(cohort = "both")))
}
short <- function(x) case_when(grepl("^dilated", x) ~ "DCM", grepl("^arrhythmogenic", x) ~ "ACM", grepl("^myocardial infarction", x) ~ "Infarction",
                               grepl("^myocarditis", x) ~ "Myocarditis", grepl("^Heart$", x) ~ "Healthy 1", TRUE ~ "Healthy 2")

# ---- (a) ligand, enzyme, receptors ----
EC_GENES <- c("EDN1", "EDN2", "EDN3", "ECE1", "ECE2", "EDNRB", "EDNRA")     # the three endothelins, their activating enzymes, both receptors
plan_a <- list(capillary = EC_GENES, arterial = EC_GENES, venous = EC_GENES,
               mural = c("EDNRA", "EDNRB"), cardiomyocyte = c("EDNRA", "EDNRB", "EDN1"), fibroblast = c("EDNRA", "EDNRB", "EDN1"))
A <- list(); Aper <- list()
for (ct in names(plan_a)) { ge <- gene_expr(ct, plan_a[[ct]])
  for (g in plan_a[[ct]]) { r <- pool(ge, g); if (!nrow(r$meta)) next
    A[[paste(ct, g)]] <- r$meta |> mutate(cell_type = ct, gene = g); Aper[[paste(ct, g)]] <- r$per |> mutate(cell_type = ct, gene = g) } }
A <- bind_rows(A); Aper <- bind_rows(Aper) |> mutate(group = short(stratum))
cat("\n== (a) Endothelin axis genes: female - male (log2 CPM), pooled ==\n")
A |> select(cell_type, gene, cohort, k, n_w, n_m, est, lo, hi, p, I2) |> r3() |> arrange(gene, cell_type, cohort) |> print(n = Inf, width = Inf)
cat("\n== (a) EDN1 in endothelium, each dataset and disease separately ==\n")
Aper |> filter(gene == "EDN1") |> select(cell_type, group, n_w, n_m, est, se, p) |> r3() |> arrange(cell_type, group) |> print(n = Inf, width = Inf)

# ---- (b) percentage of capillary ECs with any EDN1 transcript ----
dirs <- c(healthy = PAPER1_DIR, disease = getwd())
cellfrac <- bind_rows(lapply(names(dirs), function(co) {
  ec <- readRDS(file.path(dirs[[co]], "ec_subtypes.rds")); S <- C[[co]]$samples
  bind_rows(lapply(list.files(file.path(dirs[[co]], "dl"), full.names = TRUE), function(f) {
    x <- readRDS(f); s <- paste(substr(x$dataset_id, 1, 8), x$donor_id, sep = "|"); if (!s %in% S$sample) return(NULL)
    a <- ec |> filter(dataset_id == x$dataset_id, donor_id == x$donor_id, subtype == "capillary"); if (nrow(a) < MIN_CAP) return(NULL)
    m <- x$ec_counts; umi <- Matrix::colSums(m)[a$soma_joinid]
    e <- if ("EDN1" %in% rownames(m)) m["EDN1", a$soma_joinid] else rep(0, nrow(a))
    tibble(sample = s, n_cap = nrow(a), n_pos = sum(e > 0), pct_EDN1_pos = 100 * mean(e > 0), log_umi = median(log(umi)))
  }))
})) |> left_join(bind_rows(core_cols(C$healthy$samples), core_cols(C$disease$samples)) |> select(sample, stratum, cohort, sex, age, assay), by = "sample") |>
  mutate(logit_pos = log((n_pos + 0.5) / (n_cap - n_pos + 0.5)))
cat("\n== (b) Capillary ECs with any EDN1 transcript: median % per donor ==\n")
cellfrac |> group_by(cohort, sex) |> summarise(donors = n(), median_pct = median(pct_EDN1_pos), .groups = "drop") |> r3() |> print()
B <- bind_rows(pool(cellfrac, "logit_pos", fit = sex_fit)$meta |> mutate(adjusted = "none"),
               pool(cellfrac, "logit_pos", "log_umi", fit = sex_fit)$meta |> mutate(adjusted = "+ sequencing depth per cell"))
cat("\nFemale - male in the share of EDN1-positive capillary cells (SD units of logit share):\n")
B |> select(adjusted, cohort, k, n_w, n_m, est, lo, hi, p) |> r3() |> print(n = Inf, width = Inf)

# ---- (c) capillary-to-mural signalling: EDN1 (capillary) + EDNRA (mural), z-scored within dataset/disease ----
cap <- gene_expr("capillary", "EDN1") |> select(sample, stratum, cohort, sex, age, assay, EDN1)
mur <- gene_expr("mural", "EDNRA") |> select(sample, EDNRA)
pair <- inner_join(cap, mur, by = "sample") |> group_by(stratum) |>
  mutate(pair_score = as.numeric(scale(EDN1)) + as.numeric(scale(EDNRA))) |> ungroup()
Cc <- pool(pair, "pair_score", fit = sex_fit)$meta
cat("\n== (c) EDN1 (capillary) -> EDNRA (mural) signalling score: female - male (SD units) ==\n")
Cc |> select(cohort, k, n_w, n_m, est, lo, hi, p) |> r3() |> print()

# ---- (d) what tracks capillary EDN1, and does it explain the sex difference? ----
feat_sets <- list(Stress = STRESS, CM_ambient = CM_AMBIENT,
                  Hypoxia = c("VEGFA", "ADM", "SLC2A1", "PDK1", "BNIP3", "EGLN3", "ANKRD37"),
                  Shear_KLF2 = c("KLF2", "KLF4", "NOS3", "THBD"),
                  Inflammation = c("VCAM1", "ICAM1", "SELE", "CXCL2", "CCL2", "IL6"),
                  NO_eNOS = programs$NO_eNOS, Prostacyclin = programs$Prostacyclin, ETB_receptor = c("EDNRB", "EDNRB"),
                  Estrogen_response = setdiff(H$HALLMARK_ESTROGEN_RESPONSE_EARLY, "EDN1"))
Dd <- bind_rows(lapply(c("healthy", "disease"), function(co) {
  m <- C[[co]]$mats$capillary; S <- core_cols(C[[co]]$samples) |> filter(sample %in% colnames(m), n_capillary >= MIN_CAP)
  bind_rows(lapply(split(S, S$stratum), function(d) { if (nrow(d) < 6) return(NULL)
    y <- DGEList(m[, d$sample, drop = FALSE]); if (!"EDN1" %in% rownames(y)) return(NULL)
    keep <- filterByExpr(y, group = d$sex) | rownames(y) %in% c("EDN1", "EDNRB")
    lc <- cpm(calcNormFactors(y[keep, , keep.lib.sizes = FALSE]), log = TRUE, prior.count = 1)
    d$EDN1 <- lc["EDN1", ]; d$EDN1_z <- as.numeric(scale(d$EDN1))
    for (f in names(feat_sets)) d[[f]] <- if (f == "ETB_receptor") { if ("EDNRB" %in% rownames(lc)) as.numeric(scale(lc["EDNRB", ])) else NA_real_ } else
      as.numeric(scale(prog_score(lc, feat_sets[[f]])))
    d$age_z <- as.numeric(scale(d$age)); d
  }))
}))
feats <- c(names(feat_sets), "age_z")
track <- bind_rows(lapply(feats, function(f) { d <- Dd[is.finite(Dd[[f]]) & is.finite(Dd$EDN1_z), ]
  if (nrow(d) < 10) return(NULL)
  fm <- if (n_distinct(d$stratum) > 1) paste("EDN1_z ~", f, "+ stratum") else paste("EDN1_z ~", f)
  co <- tryCatch(summary(lm(as.formula(fm), data = d))$coefficients, error = function(e) NULL)
  if (is.null(co) || !f %in% rownames(co)) return(NULL)
  tibble(feature = f, slope_SD_per_SD = co[f, 1], p = co[f, 4], donors = nrow(d), strata = n_distinct(d$stratum)) }))
cat("\n== (d) What tracks capillary EDN1 across donors (within dataset/disease)? ==\n")
track |> r3() |> arrange(p) |> print()
adj <- bind_rows(pool(Dd, "EDN1")$meta |> mutate(adjusted_for = "nothing extra"),
                 bind_rows(lapply(names(feat_sets), function(f) { r <- pool(Dd[is.finite(Dd[[f]]), ], "EDN1", f)$meta
                   if (!nrow(r)) NULL else mutate(r, adjusted_for = f) })),
                 pool(Dd, "EDN1", c("Stress", "CM_ambient", "Hypoxia", "Inflammation"))$meta |> mutate(adjusted_for = "stress + ambient + hypoxia + inflammation")) |>
  filter(cohort == "both")
cat("\n== (d) EDN1 in capillaries, female - male (log2), all 138 donors, after adjusting for each feature ==\n")
adj |> select(adjusted_for, k, n_w, n_m, est, lo, hi, p) |> r3() |> print(n = Inf, width = Inf)

# ---- (e) age ----
E1 <- bind_rows(pool(Dd |> filter(age < 50), "EDN1")$meta |> mutate(ages = "< 50 years"),
                pool(Dd |> filter(age >= 50), "EDN1")$meta |> mutate(ages = ">= 50 years"))
cat("\n== (e) EDN1 in capillaries, female - male (log2), by age ==\n")
E1 |> select(ages, cohort, k, n_w, n_m, est, lo, hi, p) |> r3() |> print(n = Inf, width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(A, "tables/S3J_endothelin_genes.csv", row.names = FALSE)
write.csv(Aper, "tables/S3J_endothelin_per_stratum.csv", row.names = FALSE)
write.csv(B, "tables/S3J_EDN1_positive_cells.csv", row.names = FALSE)
write.csv(bind_rows(track |> mutate(table = "tracks EDN1"), adj |> mutate(table = "adjusted sex difference")), "tables/S3J_EDN1_correlates.csv", row.names = FALSE)
saveRDS(list(genes = A, per_stratum = Aper, positive_cells = B, cell_fractions = cellfrac, pair = Cc, tracks = track, adjusted = adj, age = E1, donors = Dd),
        "stage3_endothelin.rds")

pl <- Aper |> filter(gene == "EDN1") |> mutate(group = factor(group, levels = c("Healthy 1", "Healthy 2", "DCM", "ACM", "Infarction", "Myocarditis")),
                                               cell_type = factor(cell_type, levels = c("capillary", "arterial", "venous")))
f <- ggplot(pl, aes(est, group, colour = grepl("Healthy", group))) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), width = 0, orientation = "y", linewidth = 0.6) +
  geom_point(aes(size = n_w + n_m)) + facet_wrap(~ cell_type) + scale_size_area(max_size = 3.5, name = "Donors") +
  scale_colour_manual(values = c("TRUE" = "#4D4D4D", "FALSE" = "#E66101"), labels = c("TRUE" = "Healthy", "FALSE" = "Diseased"), name = NULL) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  labs(title = "Endothelin-1 (EDN1) in endothelial cells, by dataset and disease", x = "Female minus male (log2), 95% CI", y = NULL)
ggsave("figures/S3_Fig_endothelin.pdf", f, width = 180, height = 90, units = "mm")
ggsave("figures/S3_Fig_endothelin.png", f, width = 180, height = 90, units = "mm", dpi = 200)
