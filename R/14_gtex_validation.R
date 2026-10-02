# =============================================================================
# 14  GTEx (third, independent cohort; exploratory) - open-access bulk RNA-seq of the left ventricle (GTEx v8)
# Bulk tissue cannot isolate cell types. It is used here for three checks:
#   A. EDN1 and EDNRB (endothelial genes): female - male, with and without adjustment for capillary content,
#      and by age band (20-49 vs 50-79 years)
#   B. cardiomyocyte sex-biased genes from script 13 (bulk ventricle is mostly cardiomyocyte RNA): same direction?
#   C. younger females: capillary content, the six programs adjusted for capillary content, X-Y paralogs,
#      female - male at 20-49 and 50-79 years
# Models adjust for age band, RNA integrity (RIN), ischemic time and Hardy death class.
# Data: three GTEx v8 open-access files (about 60 MB), downloaded once into <data folder>/gtex.
# Makes: tables/S3G_gtex_*.csv, stage3_gtex.rds
# =============================================================================
dir.create("gtex", showWarnings = FALSE)
base <- "https://storage.googleapis.com/adult-gtex"
files <- c(counts = "bulk-gex/v8/rna-seq/counts-by-tissue/gene_reads_2017-06-05_v8_heart_left_ventricle.gct.gz",
           subj   = "annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
           samp   = "annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
loc <- file.path("gtex", basename(files))
options(timeout = 1800)
for (i in seq_along(files)) if (!file.exists(loc[i])) {
  ok <- tryCatch({ download.file(file.path(base, files[i]), loc[i], mode = "wb"); TRUE }, error = function(e) FALSE)
  if (!ok) stop("Could not download ", basename(files[i]), ". Download it by hand from gtexportal.org -> Downloads -> ",
                "Adult GTEx -> Bulk tissue expression / Metadata (v8), and put it in ", normalizePath("gtex"))
}
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

g <- data.table::fread(loc[1], skip = 2, data.table = FALSE)
sample_cols <- grep("^GTEX-", names(g), value = TRUE)
cnt <- rowsum(as.matrix(g[, sample_cols]), group = g$Description)          # sum duplicated gene symbols
subj <- read.delim(loc[2]); samp <- read.delim(loc[3], quote = "")
meta <- tibble(SAMPID = colnames(cnt)) |>
  mutate(SUBJID = sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID)) |>
  left_join(samp |> select(SAMPID, SMRIN, SMTSISCH), by = "SAMPID") |>
  left_join(subj, by = "SUBJID") |>
  mutate(sex = relevel(factor(ifelse(SEX == 2, "female", "male")), ref = "male"),
         band = factor(ifelse(AGE %in% c("20-29", "30-39", "40-49"), "20-49", "50-79"), levels = c("50-79", "20-49")),
         hardy = factor(ifelse(is.na(DTHHRDY), "unknown", DTHHRDY)),
         isch = SMTSISCH / 60) |>
  filter(!is.na(SMRIN), !is.na(isch))
cnt <- cnt[, meta$SAMPID]
cat("\n== GTEx left ventricle samples by sex and age ==\n")
meta |> count(AGE, sex) |> pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print()

y  <- DGEList(cnt); y <- calcNormFactors(y[filterByExpr(y, group = meta$sex), , keep.lib.sizes = FALSE])
lc <- cpm(y, log = TRUE, prior.count = 1)
CAP_BULK <- c("CA4", "RGCC", "BTNL9", "GPIHBP1")
cat("Capillary genes found:", paste(intersect(CAP_BULK, rownames(lc)), collapse = ", "), "\n")

sc <- meta |> mutate(cap_content = prog_score(lc, CAP_BULK))
for (p in names(programs)) {
  sc[[p]] <- prog_score(lc, programs[[p]])
  sc[[paste0("capcov_", p)]] <- prog_score(lc, setdiff(CAP_BULK, programs[[p]]))
}
sc$X_escape <- prog_score(lc, pos_control$X_escape)
cp <- edgeR::cpm(y)
for (i in seq_len(nrow(xy_pairs))) {
  if (!xy_pairs$X[i] %in% rownames(cp)) next
  X <- cp[xy_pairs$X[i], ]
  sc[[paste0("Xonly_", xy_pairs$X[i])]] <- log2(X + 1)
  if (xy_pairs$Y[i] %in% rownames(cp))
    sc[[paste0("XplusY_", xy_pairs$X[i], "_", xy_pairs$Y[i])]] <- log2(X + cp[xy_pairs$Y[i], ] + 1)
}

band_fit <- function(v, adj = NULL) {
  d <- sc[is.finite(sc[[v]]), ]; d$yy <- as.numeric(scale(d[[v]]))
  rhs <- paste(c("sex * band", "SMRIN", "isch", "hardy", adj), collapse = " + ")
  fit <- lm(as.formula(paste("yy ~", rhs)), data = d); df <- fit$df.residual
  b <- coef(fit); V <- vcov(fit)
  one <- function(w, label) {                    # linear combination w of coefficients
    est <- sum(w * b[names(w)]); se <- sqrt(as.numeric(t(w) %*% V[names(w), names(w)] %*% w))
    tibble(outcome = v, contrast = label, est = est, lo = est - qt(0.975, df) * se, hi = est + qt(0.975, df) * se,
           p = 2 * pt(-abs(est / se), df), lo90 = est - qt(0.95, df) * se, hi90 = est + qt(0.95, df) * se)
  }
  bind_rows(one(c(sexfemale = 1), "female - male, age 50-79"),
            one(c(sexfemale = 1, `sexfemale:band20-49` = 1), "female - male, age 20-49"),
            one(c(`sexfemale:band20-49` = 1), "difference (20-49 minus 50-79)")) |>
    mutate(n_young_f = sum(d$sex == "female" & d$band == "20-49"), n_young_m = sum(d$sex == "male" & d$band == "20-49"),
           n_old_f = sum(d$sex == "female" & d$band == "50-79"), n_old_m = sum(d$sex == "male" & d$band == "50-79"))
}
res <- bind_rows(
  band_fit("cap_content"),
  bind_rows(lapply(names(programs), function(p) band_fit(p, paste0("capcov_", p)))),
  band_fit("X_escape"),
  bind_rows(lapply(grep("^Xonly_|^XplusY_", names(sc), value = TRUE), band_fit))) |>
  mutate(within_0.8 = ifelse(grepl("^difference", contrast), NA, lo90 > -MARGIN & hi90 < MARGIN))
cat("\n== C. GTEx left ventricle: female - male (SD units) by age band ==\n")
res |> select(outcome, contrast, est, lo, hi, p, within_0.8) |>
  mutate(across(where(is.double), ~ signif(.x, 3))) |> print(n = Inf, width = Inf)
cat("\nSample sizes (20-49 F/M; 50-79 F/M):", unlist(res[1, c("n_young_f", "n_young_m", "n_old_f", "n_old_m")]), "\n")

# ---- genome-wide sex effect in bulk left ventricle (limma-voom) ----
Xg  <- model.matrix(~ sex + band + SMRIN + isch + hardy, sc)
Xg  <- Xg[, colSums(abs(Xg)) > 0, drop = FALSE]
fit <- eBayes(lmFit(voom(y, Xg), Xg))
tt  <- topTable(fit, coef = "sexfemale", number = Inf, sort.by = "none") |> tibble::rownames_to_column("gene") |>
  transmute(gene, logFC_gtex = logFC, t_gtex = t, p_gtex = P.Value, FDR_gtex = adj.P.Val) |> as_tibble()
Xc  <- cbind(Xg, cap_content = sc$cap_content)
ttc <- topTable(eBayes(lmFit(voom(y, Xc), Xc)), coef = "sexfemale", number = Inf, sort.by = "none") |>
  tibble::rownames_to_column("gene") |> transmute(gene, logFC_adj_cap = logFC, p_adj_cap = P.Value) |> as_tibble()

# ---- A. endothelin genes ----
ETG <- intersect(c("EDN1", "EDNRB", "ECE1", "EDNRA", "ACE"), rownames(lc))
cat("\n== A. Endothelin genes in GTEx left ventricle, female - male (log2); n = ", nrow(sc), " (",
    sum(sc$sex == "female"), " F, ", sum(sc$sex == "male"), " M) ==\n", sep = "")
etg <- tt |> filter(gene %in% ETG) |> left_join(ttc, by = "gene")
etg |> mutate(across(where(is.double), ~ signif(.x, 3))) |> print(width = Inf)
for (gname in ETG) sc[[paste0("gene_", gname)]] <- lc[gname, ]
et_band <- bind_rows(lapply(paste0("gene_", ETG), function(v) band_fit(v, "cap_content")))
cat("\nBy age band (SD units, adjusted for capillary content):\n")
et_band |> select(outcome, contrast, est, lo, hi, p) |> mutate(across(where(is.double), ~ signif(.x, 3))) |> print(n = Inf, width = Inf)

# ---- B. cardiomyocyte sex-biased genes from the single-cell cohorts ----
cmg <- readRDS("stage3_cm_pathways.rds")$genes |> left_join(tt, by = "gene") |>
  mutate(in_gtex = !is.na(logFC_gtex), same_direction = sign(z_comb) == sign(logFC_gtex),
         validated = same_direction & p_gtex < 0.05)
cat("\n== B. Cardiomyocyte sex-biased genes (single-cell) in GTEx bulk left ventricle ==\n")
cmg |> select(gene, direction, z_comb, replicated, logFC_gtex, p_gtex, FDR_gtex, same_direction, validated) |>
  mutate(across(where(is.double), ~ signif(.x, 3))) |> print(n = Inf, width = Inf)
v <- cmg |> filter(in_gtex)
bt <- binom.test(sum(v$same_direction), nrow(v), 0.5)
cat(sprintf("\nOf %d genes found in GTEx: %d (%.0f%%) same direction (sign test p = %.2g); %d (%.0f%%) same direction and p < 0.05.\n",
            nrow(v), sum(v$same_direction), 100 * mean(v$same_direction), bt$p.value, sum(v$validated), 100 * mean(v$validated)))
vr <- v |> filter(replicated)
cat(sprintf("Replicated genes only: %d of %d same direction; %d same direction and p < 0.05.\n",
            sum(vr$same_direction), nrow(vr), sum(vr$validated)))
ch <- chrom_map(); au <- tt |> inner_join(ch |> filter(chr_class == "autosome"), by = "gene")
cat(sprintf("Background: %.0f%% of all %d autosomal genes in GTEx have p < 0.05 for sex (large sample), so 'same direction and p < 0.05' is expected for about %.0f%% of unrelated genes.\n",
            100 * mean(au$p_gtex < 0.05), nrow(au), 50 * mean(au$p_gtex < 0.05)))
cmall <- readRDS("stage3_de.rds")$de$cardiomyocyte |> filter(chr_class == "autosome") |> inner_join(tt, by = "gene")
rho <- cor.test(cmall$z_comb, cmall$t_gtex, method = "spearman", exact = FALSE)
cat(sprintf("All %d autosomal genes: Spearman correlation between the single-cell cardiomyocyte sex effect and the GTEx sex effect = %.3f.\n",
            nrow(cmall), rho$estimate))

dir.create("tables", showWarnings = FALSE)
write.csv(res, "tables/S3G_gtex_age_bands.csv", row.names = FALSE)
write.csv(etg, "tables/S3G_gtex_endothelin.csv", row.names = FALSE)
write.csv(et_band, "tables/S3G_gtex_endothelin_age.csv", row.names = FALSE)
write.csv(cmg, "tables/S3G_gtex_cardiomyocyte_genes.csv", row.names = FALSE)
saveRDS(list(age_bands = res, endothelin = etg, endothelin_age = et_band, cm_genes = cmg, sex_all_genes = tt, meta = meta), "stage3_gtex.rds")
