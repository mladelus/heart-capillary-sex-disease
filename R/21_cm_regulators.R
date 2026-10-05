# =============================================================================
# 21  PLAN SECTION R (exploratory) - what regulates the cardiomyocyte sex-biased genes?
# R1: are targets of sex-hormone receptors or of X-Y regulators more sex-biased than targets of other transcription factors?
# R2: within each sex, does the female-versus-male pattern of the 55 genes follow the dose of a regulator?
# R3: is the sex difference of the 55 genes larger before age 50 (GTEx)?
# Run in a fresh R window:
#   source("~/Desktop/heart-capillary-sex-disease/R/00_setup.R"); source("~/Desktop/heart-capillary-sex-disease/R/21_cm_regulators.R")
# Needs stage3_de.rds, stage3_cm_pathways.rds, pb.rds (both projects), gtex/ (script 14) and internet access.
# Makes: tables/S3R_*.csv, stage3_regulators.rds
# =============================================================================
set.seed(2041)
dir.create("tables", showWarnings = FALSE); dir.create("genesets", showWarnings = FALSE)
while (sink.number() > 0) sink()
sink("tables/S3R_output.txt", split = TRUE)

cm <- readRDS("stage3_de.rds")$de$cardiomyocyte |> filter(chr_class == "autosome")
G55 <- readRDS("stage3_cm_pathways.rds")$genes |> filter(consistent)
up_f <- G55$gene[G55$direction == "higher in females"]; up_m <- G55$gene[G55$direction == "higher in males"]
cat("Cardiomyocyte genes tested (autosomal):", nrow(cm), "| consistent genes:", nrow(G55), "(", length(up_f), "higher in females,", length(up_m), "higher in males )\n")

# regulators fixed in the plan
CLASSES <- list(hormone_receptors = c("AR", "ESR1", "ESR2", "PGR"),
                XY_demethylases   = c("KDM5C", "KDM5D", "KDM6A", "UTY"),
                ZFX_ZFY           = c("ZFX", "ZFY"))
RELATED <- c("KDM5A", "KDM5B", "KDM6B")          # same families, not sex-linked; shown for comparison only

# ---- R1. transcription-factor target libraries ------------------------------------------------------------------------
lib <- list()
gtrd <- tryCatch({
  m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT:GTRD"),
                error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD"))
  m <- m[grepl("_TARGET_GENES$", m$gs_name), ]
  split(m$gene_symbol, sub("_TARGET_GENES$", "", m$gs_name))
}, error = function(e) { message("GTRD sets not available: ", conditionMessage(e)); NULL })
if (!is.null(gtrd)) lib$GTRD <- gtrd
read_enrichr <- function(name) {
  f <- file.path("genesets", paste0(name, ".gmt"))
  if (!file.exists(f) || file.size(f) < 1e4)
    tryCatch(download.file(paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", name), f, mode = "wb", quiet = TRUE),
             error = function(e) message("Could not download ", name))
  if (!file.exists(f) || file.size(f) < 1e4) return(NULL)
  L <- strsplit(readLines(f, warn = FALSE), "\t")
  term <- vapply(L, `[`, "", 1); genes <- lapply(L, function(x) unique(sub(",.*$", "", x[-(1:2)][nzchar(x[-(1:2)])])))
  human <- !grepl("mouse|rat\\b|mm9|mm10", term, ignore.case = TRUE)         # human experiments only
  tf <- toupper(sub("[ _].*$", "", term))
  tapply(genes[human], tf[human], function(z) unique(unlist(z)))
}
for (nm in c("ChEA_2022", "ENCODE_TF_ChIP-seq_2015")) { x <- read_enrichr(nm); if (!is.null(x)) lib[[nm]] <- as.list(x) }
cat("\nTarget libraries loaded:", paste(names(lib), vapply(lib, length, 1L), collapse = "; "), "\n")

bg <- cm$gene; stat_abs <- setNames(abs(cm$z_comb), cm$gene); stat_signed <- setNames(cm$z_comb, cm$gene)
in55 <- bg %in% G55$gene
C <- load_cohorts()
mexp <- function(m) { y <- calcNormFactors(DGEList(m)); rowMeans(edgeR::cpm(y, log = TRUE, prior.count = 1)) }
ex_h <- mexp(C$healthy$mats[["cardiomyocyte"]]); ex_d <- mexp(C$disease$mats[["cardiomyocyte"]])
mean_expr <- rowMeans(cbind(ex_h[bg], ex_d[bg]), na.rm = TRUE); mean_expr[!is.finite(mean_expr)] <- median(mean_expr, na.rm = TRUE)

one_set <- function(S) {
  S <- intersect(S, bg); if (length(S) < 50 || length(S) > 0.6 * length(bg)) return(NULL)
  idx <- which(bg %in% S)
  ca <- limma::cameraPR(stat_abs, list(s = idx), inter.gene.cor = 0.01)
  cs <- limma::cameraPR(stat_signed, list(s = idx), inter.gene.cor = 0.01)
  fit <- suppressWarnings(glm(in55 ~ I(bg %in% S) + scale(mean_expr), family = binomial)); co <- summary(fit)$coefficients
  tibble(targets = length(S), mean_abs_z = mean(stat_abs[idx]), mean_abs_z_rest = mean(stat_abs[-idx]),
         p_more_biased = ifelse(ca$Direction == "Up", ca$PValue / 2, 1 - ca$PValue / 2),
         signed_direction = ifelse(cs$Direction == "Up", "targets higher in females", "targets higher in males"), p_signed = cs$PValue,
         in55 = sum(in55[idx]), OR_55 = exp(co[2, 1]), p_55 = co[2, 4])
}
r1 <- bind_rows(lapply(names(lib), function(ln) {
  per <- bind_rows(lapply(names(lib[[ln]]), function(tf) { x <- one_set(lib[[ln]][[tf]]); if (is.null(x)) NULL else mutate(x, library = ln, regulator = tf, type = "single factor", members = NA_character_) }))
  cls <- bind_rows(lapply(names(CLASSES), function(cn) {
    S <- unique(unlist(lib[[ln]][intersect(CLASSES[[cn]], names(lib[[ln]]))])); if (!length(S)) return(NULL)
    x <- one_set(S); if (is.null(x)) NULL else mutate(x, library = ln, regulator = cn, type = "class (union)",
                                                      members = paste(intersect(CLASSES[[cn]], names(lib[[ln]])), collapse = ", "))
  }))
  if (!nrow(per)) return(NULL)
  rk <- function(v) vapply(v, function(q) mean(per$mean_abs_z >= q), 1)       # share of all factors at least as extreme
  bind_rows(per, cls) |> mutate(n_factors = nrow(per), top_fraction = rk(mean_abs_z))
}))
if (!is.null(r1) && nrow(r1)) {
  of_interest <- c(unlist(CLASSES), RELATED, names(CLASSES))
  cat("\n== R1. PRIMARY: classes in the GTRD library (Bonferroni 0.0167; also needs top 10% of all factors) ==\n")
  r1 |> filter(type == "class (union)") |> select(library, regulator, members, targets, mean_abs_z, mean_abs_z_rest, p_more_biased, top_fraction, n_factors, in55, OR_55, p_55, signed_direction, p_signed) |>
    r3() |> print(n = Inf, width = Inf)
  cat("\n== R1. Single regulators of interest, every library ==\n")
  r1 |> filter(type == "single factor", regulator %in% of_interest) |> arrange(library, regulator) |>
    select(library, regulator, targets, mean_abs_z, mean_abs_z_rest, p_more_biased, top_fraction, n_factors, in55, OR_55, p_55, signed_direction, p_signed) |> r3() |> print(n = Inf, width = Inf)
  cat("\n== R1. For context: the 15 factors whose targets are most sex-biased, per library ==\n")
  r1 |> filter(type == "single factor") |> group_by(library) |> arrange(desc(mean_abs_z), .by_group = TRUE) |> slice_head(n = 15) |> ungroup() |>
    select(library, regulator, targets, mean_abs_z, p_more_biased, in55, OR_55, p_55) |> r3() |> print(n = Inf, width = Inf)
  write.csv(r1, "tables/S3R_regulator_targets.csv", row.names = FALSE)
} else cat("\nR1 could not be run: no target library was available.\n")

# ---- R2. within-sex dose-response in cardiomyocyte pseudobulks ---------------------------------------------------------
pattern_score <- function(lc, f, m) {
  zs <- function(g) { g <- intersect(g, rownames(lc)); z <- t(scale(t(lc[g, , drop = FALSE]))); z[!is.finite(z)] <- NA; colMeans(z, na.rm = TRUE) }
  zs(f) - zs(m)                                           # high = female pattern
}
PRED <- list(AR = "AR", ESR1 = "ESR1", PGR = "PGR", KDM5C_KDM5D = c("KDM5C", "KDM5D"), KDM6A_UTY = c("KDM6A", "UTY"), ZFX_ZFY = c("ZFX", "ZFY"))
donor_tab <- function(coh) {
  mat <- C[[coh]]$mats[["cardiomyocyte"]]; S <- C[[coh]]$samples |> filter(sample %in% colnames(mat), !is.na(age))
  bind_rows(lapply(split(S, S$stratum), function(d) {
    if (nrow(d) < 6) return(NULL)
    m <- mat[, d$sample, drop = FALSE]; y <- DGEList(m); y <- calcNormFactors(y); lc <- edgeR::cpm(y, log = TRUE, prior.count = 1); cp <- edgeR::cpm(y)
    out <- d |> transmute(sample, stratum, cohort = coh, sex = as.character(sex), age, assay = as.character(assay), score = pattern_score(lc, up_f, up_m))
    for (p in names(PRED)) { g <- intersect(PRED[[p]], rownames(cp)); out[[p]] <- if (length(g)) log2(colSums(cp[g, , drop = FALSE]) + 1) else NA_real_ }
    expressed <- rownames(lc)[rowMeans(lc) > 1]; pool <- setdiff(expressed, c(G55$gene, unlist(PRED)))
    rnd <- vapply(1:300, function(i) pattern_score(lc, sample(pool, length(up_f)), sample(pool, length(up_m))), numeric(ncol(lc)))
    out$random <- rnd; out
  }))
}
DT <- bind_rows(donor_tab("healthy"), donor_tab("disease"))
slope <- function(d, yv, xv) {                               # standardized slope of y on regulator dose, within one sex and stratum
  keep <- is.finite(d[[xv]]) & is.finite(yv); d <- d[keep, ]; yv <- yv[keep]
  if (nrow(d) < 6 || sd(d[[xv]]) == 0) return(NULL)
  d$yy <- as.numeric(scale(yv)); d$xx <- as.numeric(scale(d[[xv]]))
  f <- if (n_distinct(d$assay) > 1) yy ~ xx + age + assay else yy ~ xx + age
  co <- tryCatch(summary(lm(f, d))$coefficients, error = function(e) NULL)
  if (is.null(co) || !"xx" %in% rownames(co) || !is.finite(co["xx", 2])) return(NULL)
  c(est = co["xx", 1], se = co["xx", 2])
}
pooled <- function(per) { per <- per[is.finite(per$est) & is.finite(per$se) & per$se > 0, ]; if (nrow(per) < 2) return(c(est = NA, se = NA, k = nrow(per)))
  w <- 1 / per$se^2; c(est = sum(w * per$est) / sum(w), se = sqrt(1 / sum(w)), k = nrow(per)) }
r2 <- bind_rows(lapply(c("female", "male"), function(sx) bind_rows(lapply(names(PRED), function(p) {
  groups <- split(DT[DT$sex == sx, ], DT$stratum[DT$sex == sx])
  real <- pooled(bind_rows(lapply(groups, function(d) { s <- slope(d, d$score, p); if (is.null(s)) NULL else tibble(est = s[1], se = s[2]) })))
  null <- vapply(1:300, function(i) pooled(bind_rows(lapply(groups, function(d) { s <- slope(d, d$random[, i], p); if (is.null(s)) NULL else tibble(est = s[1], se = s[2]) })))[1], 1)
  tibble(sex = sx, regulator = p, strata = real[["k"]], donors = sum(vapply(groups, nrow, 1L)), slope_SD = real[["est"]], se = real[["se"]],
         p_model = 2 * pnorm(-abs(real[["est"]] / real[["se"]])), p_vs_random_gene_sets = mean(abs(null) >= abs(real[["est"]]), na.rm = TRUE),
         null_sd = sd(null, na.rm = TRUE))
}))))
r2 <- r2 |> mutate(p_holm = p.adjust(p_vs_random_gene_sets, "holm"))
cat("\n== R2. Within each sex: slope of the female-pattern score of the 55 genes on regulator dose (SD per SD; age-adjusted; pooled across strata) ==\n")
cat("   Positive = more regulator goes with a more female pattern. Empirical p against 300 random gene sets of the same size.\n")
r2 |> r3() |> print(n = Inf, width = Inf)
write.csv(r2, "tables/S3R_within_sex_dose.csv", row.names = FALSE)

# ---- R3. GTEx: is the sex difference of the 55 genes larger before age 50? -----------------------------------------------
loc <- file.path("gtex", c("gene_reads_2017-06-05_v8_heart_left_ventricle.gct.gz", "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                           "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
r3res <- NULL
if (all(file.exists(loc))) {
  g <- data.table::fread(loc[1], skip = 2, data.table = FALSE); sample_cols <- grep("^GTEX-", names(g), value = TRUE)
  cnt <- rowsum(as.matrix(g[, sample_cols]), group = g$Description); rm(g)
  subj <- read.delim(loc[2]); samp <- read.delim(loc[3], quote = "")
  meta <- tibble(SAMPID = colnames(cnt)) |> mutate(SUBJID = sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID)) |>
    left_join(samp |> select(SAMPID, SMRIN, SMTSISCH), by = "SAMPID") |> left_join(subj, by = "SUBJID") |>
    mutate(sex = relevel(factor(ifelse(SEX == 2, "female", "male")), ref = "male"),
           band = factor(ifelse(AGE %in% c("20-29", "30-39", "40-49"), "20-49", "50-79"), levels = c("50-79", "20-49")),
           hardy = factor(ifelse(is.na(DTHHRDY), "unknown", DTHHRDY)), isch = SMTSISCH / 60) |> filter(!is.na(SMRIN), !is.na(isch))
  cnt <- cnt[, meta$SAMPID]; y <- DGEList(cnt); y <- calcNormFactors(y[filterByExpr(y, group = meta$sex), , keep.lib.sizes = FALSE]); lc <- cpm(y, log = TRUE, prior.count = 1)
  meta$score <- pattern_score(lc, up_f, up_m)
  pool <- setdiff(rownames(lc), c(G55$gene, unlist(PRED)))
  fitb <- function(v) { d <- meta; d$yy <- as.numeric(scale(v)); fit <- lm(yy ~ sex * band + SMRIN + isch + hardy, d); b <- coef(fit); V <- vcov(fit); df <- fit$df.residual
    one <- function(w, lab) { est <- sum(w * b[names(w)]); se <- sqrt(as.numeric(t(w) %*% V[names(w), names(w)] %*% w)); tibble(contrast = lab, est = est, lo = est - qt(.975, df) * se, hi = est + qt(.975, df) * se, p = 2 * pt(-abs(est / se), df)) }
    bind_rows(one(c(sexfemale = 1), "female - male, age 50-79"), one(c(sexfemale = 1, `sexfemale:band20-49` = 1), "female - male, age 20-49"),
              one(c(`sexfemale:band20-49` = 1), "difference (20-49 minus 50-79)")) }
  r3res <- fitb(meta$score)
  nullint <- vapply(1:300, function(i) fitb(pattern_score(lc, sample(pool, length(up_f)), sample(pool, length(up_m))))$est[3], 1)
  r3res$p_vs_random_gene_sets <- c(NA, NA, mean(abs(nullint) >= abs(r3res$est[3])))
  cat("\n== R3. GTEx left ventricle: female-pattern score of the 55 genes, female - male (SD units), by age band ==\n")
  cat("   n (20-49 F/M; 50-79 F/M):", sum(meta$sex == "female" & meta$band == "20-49"), sum(meta$sex == "male" & meta$band == "20-49"),
      sum(meta$sex == "female" & meta$band == "50-79"), sum(meta$sex == "male" & meta$band == "50-79"), "\n")
  r3res |> r3() |> print(width = Inf)
  cat("\n   Within females, score by age decade (mean):\n"); meta |> filter(sex == "female") |> group_by(AGE) |> summarise(n = n(), mean_score = mean(score)) |> r3() |> print()
  cat("   Within males, score by age decade (mean):\n"); meta |> filter(sex == "male") |> group_by(AGE) |> summarise(n = n(), mean_score = mean(score)) |> r3() |> print()
  write.csv(r3res, "tables/S3R_gtex_age.csv", row.names = FALSE)
} else cat("\nR3 skipped: GTEx files not found in ", normalizePath("gtex", mustWork = FALSE), " (run script 14 first).\n")

saveRDS(list(targets = r1, dose = r2, gtex_age = r3res, donors = DT |> select(-random)), "stage3_regulators.rds")
sink()
message("Finished. Output saved in tables/S3R_output.txt")
