# K3  Plan section Q: the tests that need limma and the Hallmark gene sets (Q1, Q2, Q4)
# Input: per-donor count tables in ~/Documents/heart_capillary_disease/kidney/for_R (written by Claude from the K2 export).
# Run:  source("~/Desktop/heart-capillary-sex-disease/kidney/K3_limma_tests.R")
# Makes: kidney/tables/S3Q1_*.csv, S3Q2_limma_*.csv, S3Q4_*.csv and K3_output.txt

OUT <- path.expand("~/Documents/heart_capillary_disease/kidney"); IN <- file.path(OUT, "for_R")
suppressPackageStartupMessages({ library(limma); library(edgeR); library(metafor); library(dplyr); library(tidyr) })
options(width = 220, dplyr.summarise.inform = FALSE); set.seed(2061)
while (sink.number() > 0) sink()
sink(file.path(OUT, "tables", "K3_output.txt"), split = TRUE)
r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))
don <- read.csv(file.path(IN, "donors.csv"), stringsAsFactors = FALSE)
don$sex <- relevel(factor(don$sex), ref = "male"); don$dis <- factor(don$dis, levels = c("healthy", "CKD", "AKI"))
getct <- function(name, ds, min_cells = 30, independent = TRUE) {       # counts and donor table for one cell type and dataset
  f <- file.path(IN, sprintf("counts_%s_%s.csv.gz", name, ds)); if (!file.exists(f)) return(NULL)
  m <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE)); nc <- read.csv(file.path(IN, sprintf("cells_%s_%s.csv", name, ds)), colClasses = c("character", "numeric"))
  d <- don[don$dataset == ds, ]; d <- d[match(colnames(m), d$donor_id), ]; d$cells <- nc$cells[match(d$donor_id, nc$donor_id)]
  keep <- !is.na(d$donor_id) & d$cells >= min_cells & !is.na(d$age); if (ds == "scRNA" && independent) keep <- keep & !d$in_snRNA
  list(m = m[, keep, drop = FALSE], d = droplevels(d[keep, ]))
}
msig <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
  m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  split(m$gene_symbol, m$gs_name)
}
H <- msig(); IFNG <- H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]; ANDR <- H[["HALLMARK_ANDROGEN_RESPONSE"]]
AMB <- c("LRP2", "CUBN", "SLC34A1", "UMOD")

# ================= Q1: sex x CKD interaction in interferon-gamma signaling =================
q1 <- function(name, ds, min_cells, independent, ambient = FALSE) {
  x <- getct(name, ds, min_cells, independent); if (is.null(x)) return(NULL)
  d <- x$d[x$d$dis %in% c("healthy", "CKD"), ]; d$dis <- droplevels(d$dis); tab <- table(d$sex, d$dis)
  base <- tibble(population = name, dataset = ds, donors = ifelse(ds == "scRNA" && independent, "not in single-nucleus set", "all"), min_cells, ambient_adjusted = ambient,
                 healthy_F = tab["female", "healthy"], healthy_M = tab["male", "healthy"], CKD_F = tab["female", "CKD"], CKD_M = tab["male", "CKD"])
  if (any(dim(tab) < 2) || min(tab) < 3) return(mutate(base, note = "not testable: fewer than 3 donors in a sex x disease group"))
  y <- DGEList(x$m[, d$donor_id]); amb <- log2(colSums(cpm(y)[intersect(AMB, rownames(y)), , drop = FALSE]) + 1)
  y <- calcNormFactors(y[filterByExpr(y, group = interaction(d$sex, d$dis)), , keep.lib.sizes = FALSE]); d$amb <- amb
  X <- if (ambient) model.matrix(~ sex * dis + age + amb, d) else model.matrix(~ sex * dis + age, d)
  v <- voom(y, X); idx <- which(rownames(v) %in% IFNG); co <- "sexfemale:disCKD"
  fr <- fry(v, index = list(IFNG = idx), design = X, contrast = which(colnames(X) == co))
  cam <- camera(v, index = list(IFNG = idx), design = X, contrast = which(colnames(X) == co))
  fit <- eBayes(lmFit(v, X)); tt <- topTable(fit, coef = co, number = Inf, sort.by = "none")
  lc <- cpm(y, log = TRUE, prior.count = 1); z <- t(scale(t(lc[idx, , drop = FALSE]))); d$score <- colMeans(z, na.rm = TRUE)
  sc <- summary(lm(as.formula(paste("scale(score) ~ sex * dis + age", if (ambient) "+ amb" else "")), d))$coefficients
  mutate(base, genes_in_set = length(idx), fry_direction = fr$Direction, fry_p = fr$PValue, fry_p_mixed = fr$PValue.Mixed, camera_direction = cam$Direction, camera_p = cam$PValue,
         mean_t_set = mean(tt$t[idx]), mean_t_other = mean(tt$t[-idx]), score_interaction_SD = sc[co, 1], score_p = sc[co, 4],
         score_sex_in_healthy = sc["sexfemale", 1], score_sex_in_CKD = sc["sexfemale", 1] + sc[co, 1], note = "")
}
cat("== Q1. Sex x CKD interaction, Hallmark interferon-gamma response (interaction < 0 means the CKD-associated change is lower in women) ==\n")
cat("   Predicted direction: Down. The single-nucleus rows re-run the earlier analysis; the single-cell rows use donors not in the single-nucleus set.\n")
grid <- expand.grid(name = c("EC_glomerular_capillary", "EC_peritubular_capillary", "EC"), ds = c("snRNA", "scRNA"), min_cells = c(30, 20), ambient = c(FALSE, TRUE), stringsAsFactors = FALSE)
Q1 <- bind_rows(lapply(seq_len(nrow(grid)), function(i) tryCatch(q1(grid$name[i], grid$ds[i], grid$min_cells[i], TRUE, grid$ambient[i]), error = function(e) tibble(population = grid$name[i], dataset = grid$ds[i], note = conditionMessage(e)))))
Q1 |> select(population, dataset, min_cells, ambient_adjusted, healthy_F, healthy_M, CKD_F, CKD_M, fry_direction, fry_p, camera_p, score_interaction_SD, score_p, score_sex_in_healthy, score_sex_in_CKD, note) |> r3() |> print(n = Inf, width = Inf)
write.csv(Q1, file.path(OUT, "tables", "S3Q1_interferon_interaction.csv"), row.names = FALSE)

# ================= Q2: replicated autosomal sex-biased genes per cell type (limma-voom, as in the heart work) =================
chr <- tryCatch({ gg <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, columns = c("gene_name", "seq_name"), return.type = "data.frame")
  setNames(as.character(gg$seq_name), gg$gene_name)[!duplicated(gg$gene_name)] }, error = function(e) NULL)
if (is.null(chr)) cat("\nQ2 skipped: EnsDb.Hsapiens.v86 not available.\n") else {
  TYPES <- c("PT", "TAL", "DCT_CNT", "PC", "IC", "EC", "stroma", "immune", "EC_peritubular_capillary", "EC_glomerular_capillary")
  de1 <- function(x, groups, perm = FALSE) {                         # sex effect within disease group, inverse-variance pooled across groups
    per <- lapply(groups, function(g) { d <- x$d[x$d$dis == g, ]; if (sum(d$sex == "female") < 3 || sum(d$sex == "male") < 3) return(NULL)
      if (perm) d$sex <- sample(d$sex)
      y <- DGEList(x$m[, d$donor_id]); y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE]); X <- model.matrix(~ sex + age, d)
      tt <- topTable(eBayes(lmFit(voom(y, X), X)), coef = "sexfemale", number = Inf, sort.by = "none"); data.frame(gene = rownames(tt), logFC = tt$logFC, se = tt$logFC / tt$t, n = nrow(d)) })
    per <- per[!vapply(per, is.null, TRUE)]; if (!length(per)) return(NULL)
    a <- bind_rows(per) |> group_by(gene) |> filter(n() == length(per)) |> summarise(est = sum(logFC / se^2) / sum(1 / se^2), se = sqrt(1 / sum(1 / se^2)), n = sum(n))
    a |> mutate(z = est / se, p = 2 * pnorm(-abs(z)), FDR = p.adjust(p, "BH"), cls = ifelse(chr[gene] %in% c("X", "Y"), "XY", ifelse(chr[gene] %in% as.character(1:22), "autosome", "other")))
  }
  one <- function(ct, groups, perm = FALSE) {
    a <- getct(ct, "snRNA"); b <- getct(ct, "scRNA"); if (is.null(a) || is.null(b)) return(NULL)
    da <- de1(a, groups, perm); db <- de1(b, groups, perm); if (is.null(da) || is.null(db)) return(NULL)
    j <- inner_join(da, db, by = "gene", suffix = c("_a", "_b")); au <- j[j$cls_a == "autosome", ]; same <- sign(au$z_a) == sign(au$z_b)
    rep <- (au$FDR_a < 0.05 & au$p_b < 0.05 & same) | (au$FDR_b < 0.05 & au$p_a < 0.05 & same); xy <- j[j$cls_a == "XY", ]
    list(row = tibble(cell_type = ct, n_sn = max(da$n), n_sc = max(db$n), autosomal_genes = nrow(au), rho = cor(au$z_a, au$z_b, method = "spearman"), FDR05_sn = sum(au$FDR_a < 0.05), FDR05_sc = sum(au$FDR_b < 0.05),
                      replicated = sum(rep), XY_FDR05_both = sum(xy$FDR_a < 0.05 & xy$FDR_b < 0.05 & sign(xy$z_a) == sign(xy$z_b))), genes = au[rep, c("gene", "est_a", "FDR_a", "est_b", "FDR_b")])
  }
  for (gs in list(c("healthy"), c("healthy", "CKD"))) {
    lab <- paste(gs, collapse = " + ")
    cat("\n== Q2. Replicated autosomal sex-biased genes (limma-voom), donors:", lab, if (length(gs) > 1) "[not pre-specified: added for power]" else "[as planned]", "==\n")
    real <- lapply(TYPES, one, groups = gs); names(real) <- TYPES; real <- real[!vapply(real, is.null, TRUE)]
    tab <- bind_rows(lapply(real, `[[`, "row")); NP <- 20
    nul <- sapply(names(real), function(ct) replicate(NP, { r <- one(ct, gs, perm = TRUE); if (is.null(r)) NA else r$row$replicated }))
    tab$null_replicated_max <- apply(nul, 2, max, na.rm = TRUE); tab$null_replicated_median <- apply(nul, 2, median, na.rm = TRUE)
    tab |> arrange(desc(replicated)) |> r3() |> print(n = Inf, width = Inf)
    for (ct in names(real)) if (nrow(real[[ct]]$genes)) { cat("  ", ct, ":", paste(sprintf("%s(%s)", real[[ct]]$genes$gene, ifelse(real[[ct]]$genes$est_a > 0, "F", "M")), collapse = ", "), "\n") }
    write.csv(tab, file.path(OUT, "tables", paste0("S3Q2_limma_", gsub("[^A-Za-z]", "", lab), ".csv")), row.names = FALSE)
    write.csv(bind_rows(lapply(names(real), function(ct) if (nrow(real[[ct]]$genes)) cbind(cell_type = ct, real[[ct]]$genes))), file.path(OUT, "tables", paste0("S3Q2_limma_genes_", gsub("[^A-Za-z]", "", lab), ".csv")), row.names = FALSE)
  }
}

# ================= Q4: androgen-response score by sex =================
cat("\n== Q4. Hallmark androgen-response score, women minus men (SD units; age-adjusted; within dataset x disease; pooled) ==\n")
q4 <- function(name) {
  per <- bind_rows(lapply(c("snRNA", "scRNA"), function(ds) { x <- getct(name, ds); if (is.null(x) || ncol(x$m) < 6) return(NULL)
    y <- calcNormFactors(DGEList(x$m)); lc <- cpm(y, log = TRUE, prior.count = 1); g <- intersect(ANDR, rownames(lc)[rowMeans(lc) > 1]); z <- t(scale(t(lc[g, , drop = FALSE]))); x$d$score <- colMeans(z, na.rm = TRUE)
    bind_rows(lapply(levels(x$d$dis), function(gr) { d <- x$d[x$d$dis == gr, ]; if (sum(d$sex == "female") < 3 || sum(d$sex == "male") < 3) return(NULL)
      co <- summary(lm(scale(score) ~ sex + age, d))$coefficients; tibble(stratum = paste(gr, ds), est = co["sexfemale", 1], se = co["sexfemale", 2], F = sum(d$sex == "female"), M = sum(d$sex == "male"), genes = length(g)) })) }))
  if (!nrow(per)) return(NULL)
  m <- if (nrow(per) >= 3) tryCatch(rma(yi = per$est, sei = per$se, method = "REML"), error = function(e) rma(yi = per$est, sei = per$se, method = "DL")) else rma(yi = per$est, sei = per$se, method = "FE")
  tibble(population = name, strata = nrow(per), F = sum(per$F), M = sum(per$M), genes = max(per$genes), est = m$b[1], lo = m$ci.lb, hi = m$ci.ub, p = m$pval)
}
Q4 <- bind_rows(lapply(c("EC_peritubular_capillary", "EC_glomerular_capillary", "EC_arteriole", "EC_ascending_vasa_recta", "EC", "PT"), q4))
Q4$p_holm <- NA; i <- Q4$population != "PT"; Q4$p_holm[i] <- p.adjust(Q4$p[i], "holm")
Q4 |> r3() |> print(n = Inf, width = Inf); cat("   (PT, proximal tubule, is the positive control and is not part of the Holm correction.)\n")
write.csv(Q4, file.path(OUT, "tables", "S3Q4_androgen_response.csv"), row.names = FALSE)
sink()
message("\nFinished. Output saved in ", file.path(OUT, "tables", "K3_output.txt"))
