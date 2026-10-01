# =============================================================================
# SETUP - run first in every new R session
# Paper 2: sex differences in human heart capillaries in DISEASE (CELLxGENE Census)
# Same helpers, gene sets and thresholds as paper 1 (heart-capillary-sex-differences).
# =============================================================================

# Data folder for downloads and saved results. Change here or set the env variable HCAP_DATA.
DATA_DIR   <- Sys.getenv("HCAP2_DATA", "~/Documents/heart_capillary_disease")
PAPER1_DIR <- path.expand(Sys.getenv("HCAP_DATA", "~/Documents/heart_capillary_project"))  # healthy donors (paper 1)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(DATA_DIR)
dir.create("figures", showWarnings = FALSE)
CENSUS_VERSION <- "2025-11-08"                 # pinned; same release as the earlier project

suppressPackageStartupMessages({
  library(cellxgene.census)
  library(Seurat); library(Matrix)
  library(edgeR); library(limma)
  library(metafor)                              # random-effects meta-analysis
  library(lme4); library(lmerTest)
  library(ggplot2); library(tidyr); library(stringr); library(dplyr)
})
filter <- dplyr::filter; select <- dplyr::select; count <- dplyr::count
first  <- dplyr::first;  rename <- dplyr::rename; desc <- dplyr::desc

census <- open_soma(census_version = CENSUS_VERSION)

# ---- Pre-specified settings (see ANALYSIS_PLAN.md) ---------------------------
MIN_EC      <- 100     # donor eligibility: >= 100 endothelial cells in ventricular tissue (Census labels)
MIN_CAP     <- 30      # >= 30 capillary ECs for a donor's capillary pseudobulk
MIN_CT      <- 30      # >= 30 cells for any other cell-type pseudobulk
MIN_PER_SEX <- 3       # a dataset enters a meta-analysis only with >= 3 women and >= 3 men
MARGIN      <- 0.8     # smallest effect of interest (primary equivalence bounds +/- 0.8 donor-level SD = "large")
MARGIN2     <- 0.5     # secondary, stricter bounds (+/- 0.5 SD), reported alongside
# Ventricular myocardium only (atria and valves excluded). Checked against the inventory in 01.
VENT <- paste0("^(heart left ventricle|left ventricle|left cardiac ventricle|anterior wall of left ventricle|",
               "apical region of left ventricle|apex of heart|interventricular septum)$")

# ---- Gene IDs -> symbols ------------------------------------------------------
if (!file.exists("gene_map.rds")) {
  genes <- census$get("census_data")$get("homo_sapiens")$ms$get("RNA")$var$read(
    column_names = c("feature_id", "feature_name"))$concat() |> as.data.frame()
  saveRDS(genes, "gene_map.rds")
}
gene_map <- readRDS("gene_map.rds")

to_symbols <- function(m) {                     # sum Ensembl IDs that share a symbol
  sym <- gene_map$feature_name[match(rownames(m), gene_map$feature_id)]
  f <- factor(sym)
  out <- fac2sparse(f) %*% m
  rownames(out) <- levels(f)
  out
}

# ---- Age from Census development_stage ----------------------------------------
dec_mid <- c("third decade stage" = 25, "fourth decade stage" = 35, "fifth decade stage" = 45,
             "sixth decade stage" = 55, "seventh decade stage" = 65, "eighth decade stage" = 75,
             "ninth decade stage" = 85)
parse_age <- function(x) {
  x <- as.character(x)
  coalesce(as.numeric(str_extract(x, "^\\d+(?=[- ]year)")), unname(dec_mid[x]))
}

# ---- Broad cell classes from Census cell_type (order matters) -------------------
cell_class <- function(ct) {
  ct <- tolower(as.character(ct))
  case_when(
    str_detect(ct, "lymphatic")                          ~ "lymphatic EC",
    str_detect(ct, "endocardial")                        ~ "endocardial",
    str_detect(ct, "endothelial")                        ~ "blood EC",
    str_detect(ct, "pericyte")                           ~ "pericyte",
    str_detect(ct, "smooth muscle")                      ~ "smooth muscle",
    str_detect(ct, "cardiac muscle|cardiomyocyte|myocyte") ~ "cardiomyocyte",
    str_detect(ct, "fibroblast")                         ~ "fibroblast",
    str_detect(ct, "macrophage|monocyte|dendritic|myeloid|mast") ~ "myeloid",
    str_detect(ct, "t cell|b cell|natural killer|lymphocyte|plasma") ~ "lymphoid",
    str_detect(ct, "adipocyte|fat cell")                 ~ "adipocyte",
    str_detect(ct, "neur|glial|schwann")                 ~ "neural",
    str_detect(ct, "mesothelial")                        ~ "mesothelial",
    TRUE                                                 ~ "other")
}

# ---- EC subtype markers (same rules in every dataset; see 03) -------------------
# Capillary markers deliberately exclude genes used in the outcome programs (e.g. CD36, FABP4).
ec_markers <- list(capillary = c("CA4", "RGCC", "BTNL9"),
                   arterial  = c("GJA5", "HEY1", "SEMA3G", "GJA4", "DKK2", "FBLN5"),
                   venous    = c("ACKR1", "NR2F2", "CPE"))

# ---- Outcome programs (pre-specified; ANALYSIS_PLAN.md) ------------------------
programs <- list(
  NO_eNOS        = c("NOS3", "KLF2", "KLF4", "CAV1", "GCH1", "SLC7A1", "DDAH1", "DDAH2"),
  Endothelin_ACE = c("EDN1", "ECE1", "EDNRB", "ACE"),
  Prostacyclin   = c("PTGS1", "PTGS2", "PTGIS", "PLA2G4A"),
  Barrier        = c("CLDN5", "CDH5", "ESAM", "OCLN", "TJP1", "JAM2", "PECAM1"),
  FA_transport   = c("CD36", "FABP4", "FABP5", "LPL", "GPIHBP1"),
  Angiogenic_tip = c("ESM1", "APLN", "DLL4", "ANGPT2", "PGF", "KCNE3"))
PRIMARY <- names(programs)

# Positive control: X-inactivation escape genes (must come out higher in women)
pos_control <- list(X_escape = c("KDM6A", "KDM5C", "DDX3X", "EIF1AX", "ZFX", "USP9X", "JPX"))

# Dissociation / ischemia stress (sensitivity covariate)
CM_AMBIENT <- c("MYH7", "MYH6", "MYL2", "MYL3", "TNNT2", "TNNI3", "ACTC1", "MB", "TTN", "RYR2", "NPPA", "NPPB",
                "CKM", "COX6A2", "MYL7")
peri_programs <- list(
  Pericyte_identity     = c("RGS5", "NOTCH3", "PDGFRB", "CSPG4", "HIGD1B"),
  Contractile           = c("ACTA2", "TAGLN", "MYL9", "CNN1", "MYH11"),
  KATP_channel          = c("KCNJ8", "ABCC9"),
  Constrictor_receptors = c("EDNRA", "AGTR1"),
  NO_cGMP_response      = c("GUCY1A1", "GUCY1B1", "PRKG1", "PDE5A"))
STRESS <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "IER2", "IER3", "DUSP1",
            "ZFP36", "HSPA1A", "HSPA1B", "HSPA8", "HSPH1", "DNAJB1", "HSP90AA1", "HSPE1", "SOCS3", "NR4A1")

# Capillary EC -> pericyte (and back) signaling pairs: sender gene, receiver gene
lr_pairs <- tibble::tribble(
  ~pair,             ~ligand,  ~sender,       ~receptor, ~receiver,
  "PDGFB-PDGFRB",    "PDGFB",  "capillary",   "PDGFRB",  "pericyte",
  "JAG1-NOTCH3",     "JAG1",   "capillary",   "NOTCH3",  "pericyte",
  "EDN1-EDNRA",      "EDN1",   "capillary",   "EDNRA",   "pericyte",
  "ANGPT1-TEK",      "ANGPT1", "pericyte",    "TEK",     "capillary",
  "ANGPT2-TEK",      "ANGPT2", "capillary",   "TEK",     "capillary")

# X-Y paralog pairs
xy_pairs <- tibble::tibble(
  X = c("KDM6A", "KDM5C", "USP9X", "DDX3X", "EIF1AX", "ZFX", "RPS4X", "NLGN4X"),
  Y = c("UTY", "KDM5D", "USP9Y", "DDX3Y", "EIF1AY", "ZFY", "RPS4Y1", "NLGN4Y"))

# ---- Shared helpers -----------------------------------------------------------
logcpm <- function(m) edgeR::cpm(m, log = TRUE, prior.count = 1)

# program score per sample: mean z-score (across samples in one dataset) of the program genes
prog_score <- function(lc, genes) {
  g <- intersect(genes, rownames(lc))
  if (length(g) < 2) return(rep(NA_real_, ncol(lc)))
  z <- t(scale(t(lc[g, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

# women - men difference, in SD units of the outcome, adjusted for age (and assay when it varies)
sex_fit <- function(d, y, extra = NULL) {
  d <- droplevels(d[is.finite(d[[y]]) & !is.na(d$age), ])
  if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX) return(NULL)
  d$yy <- as.numeric(scale(d[[y]]))
  rhs <- c("sex", "age", extra)
  if ("assay" %in% names(d) && n_distinct(d$assay) > 1) rhs <- c(rhs, "assay")
  fit <- lm(as.formula(paste("yy ~", paste(rhs, collapse = " + "))), data = d)
  co <- summary(fit)$coefficients
  if (!"sexfemale" %in% rownames(co)) return(NULL)
  tibble(est = co["sexfemale", 1], se = co["sexfemale", 2], p = co["sexfemale", 4],
         n_w = sum(d$sex == "female"), n_m = sum(d$sex == "male"))
}

# random-effects meta-analysis + equivalence (TOST) against +/- MARGIN
meta_one <- function(r) {
  r <- r[is.finite(r$est) & is.finite(r$se), ]
  if (nrow(r) == 0) return(tibble())
  m <- if (nrow(r) >= 3) tryCatch(rma(yi = r$est, sei = r$se, method = "REML"),
                                  error = function(e) rma(yi = r$est, sei = r$se, method = "DL"))
       else rma(yi = r$est, sei = r$se, method = "FE")
  pr <- predict(m)
  ci90 <- m$b[1] + c(-1, 1) * qnorm(0.95) * m$se
  tibble(k = nrow(r), n_w = sum(r$n_w), n_m = sum(r$n_m),
         est = m$b[1], lo = m$ci.lb, hi = m$ci.ub, p = m$pval,
         pi_lo = if (!is.null(pr$pi.lb)) pr$pi.lb else NA, pi_hi = if (!is.null(pr$pi.ub)) pr$pi.ub else NA,
         I2 = if (!is.null(m$I2)) m$I2 else NA, tau2 = m$tau2,
         lo90 = ci90[1], hi90 = ci90[2],
         p_tost = max(pnorm((m$b[1] + MARGIN) / m$se, lower.tail = FALSE),
                      pnorm((m$b[1] - MARGIN) / m$se)),
         within_0.5 = ci90[1] > -MARGIN2 & ci90[2] < MARGIN2,
         same_dir = sum(sign(r$est) == sign(m$b[1])))
}

verdict <- function(p_adj, lo90, hi90, same_dir = 1, k = 1) case_when(
  p_adj < 0.05 & same_dir > k / 2       ~ "difference",
  lo90 > -MARGIN & hi90 < MARGIN        ~ paste0("equivalent (within +/-", MARGIN, " SD)"),
  TRUE                                  ~ "inconclusive")

# ---- Stage 3 helpers (scripts 08-11) ------------------------------------------
# Both cohorts: stage 1 = healthy (paper 1 data folder), stage 2 = disease (this folder).
# Each sample gets a "stratum": dataset (healthy) or disease (stage 2); normal donors of stage 2 are dropped.
CELLTYPES <- list(capillary = c("capillary", "capillary"), arterial = c("arterial", "arterial"),
                  venous = c("venous", "venous"), mural = c("pericyte", "mural"),
                  fibroblast = c("fibroblast", "fibroblast"), cardiomyocyte = c("cardiomyocyte", "cardiomyocyte"),
                  myeloid = c("myeloid", "myeloid"), lymphoid = c("lymphoid", "lymphoid"))
load_cohorts <- function() {
  h <- readRDS(file.path(PAPER1_DIR, "pb.rds")); d <- readRDS("pb.rds")
  h$samples <- h$samples |> mutate(stratum = substr(as.character(dataset_title), 1, 40), cohort = "healthy")
  d$samples <- d$samples |> filter(!grepl("^normal", stratum)) |> mutate(cohort = "disease")
  list(healthy = h, disease = d)
}
chrom_map <- function() {
  if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)
  gg <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, columns = c("gene_name", "seq_name"), return.type = "data.frame")
  tibble(gene = gg$gene_name, chr = as.character(gg$seq_name)) |>
    filter(chr %in% c(1:22, "X", "Y")) |> distinct(gene, .keep_all = TRUE) |>
    mutate(chr_class = case_when(chr == "X" ~ "X", chr == "Y" ~ "Y", TRUE ~ "autosome"))
}
# limma-voom per stratum (sex + age [+ assay]), inverse-variance meta of logFC across strata
de_meta <- function(mat, S, min_k = 2) {
  S <- S |> filter(sample %in% colnames(mat))
  used <- list()
  per <- bind_rows(lapply(split(S, S$stratum), function(d) {
    d <- droplevels(d |> filter(!is.na(age)))
    if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX) return(NULL)
    y <- DGEList(mat[, d$sample, drop = FALSE]); y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE])
    X <- if (n_distinct(d$assay) > 1) model.matrix(~ sex + age + assay, d) else model.matrix(~ sex + age, d)
    X <- X[, colSums(abs(X)) > 0 & !duplicated(t(X)), drop = FALSE]
    if (ncol(X) >= nrow(X)) return(NULL)
    fit <- eBayes(lmFit(voom(y, X), X)); tt <- topTable(fit, coef = "sexfemale", number = Inf, sort.by = "none")
    used[[d$stratum[1]]] <<- d |> select(sample, stratum, sex)
    tibble(stratum = d$stratum[1], n = nrow(d), gene = rownames(tt), logFC = tt$logFC, se = tt$logFC / tt$t)
  }))
  if (!nrow(per)) return(tibble())
  out <- per |> filter(is.finite(se), se > 0) |> group_by(gene) |> filter(n() >= min(min_k, n_distinct(per$stratum))) |>
    summarise(k = n(), logFC = sum(logFC / se^2) / sum(1 / se^2), se = sqrt(1 / sum(1 / se^2)), .groups = "drop") |>
    mutate(z = logFC / se, p = 2 * pnorm(-abs(z)), FDR = p.adjust(p, "BH"))
  attr(out, "donors") <- bind_rows(used)            # donors that actually entered the model
  out
}
hallmark_sets <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
  if (packageVersion("msigdbr") >= "10.0.0" && !requireNamespace("msigdbdf", quietly = TRUE))
    install.packages("msigdbdf", repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org"))
  m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  split(m$gene_symbol, m$gs_name)
}
r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))
core_cols <- function(S) S |> transmute(sample, stratum, cohort, sex = as.character(sex), age, assay = as.character(assay),
  n_capillary, n_EC, n_ratio_cells = if ("n_mural" %in% names(S)) n_mural else n_pericyte) |>
  mutate(sex = relevel(factor(sex), ref = "male"))

message("Setup complete (paper 2). Census ", CENSUS_VERSION, " open. Data folder: ", getwd(),
        "\nPaper 1 folder (healthy donors): ", PAPER1_DIR, if (file.exists(file.path(PAPER1_DIR, "pb.rds"))) " [found]" else " [NOT FOUND]")

# ---- One-time installation (run manually if a package is missing) -------------
# install.packages(c("BiocManager", "dplyr", "tidyr", "stringr", "ggplot2", "Matrix", "Seurat",
#                    "lme4", "lmerTest", "metafor", "patchwork"))
# install.packages("cellxgene.census",
#                  repos = c("https://chanzuckerberg.r-universe.dev", "https://cloud.r-project.org"))
# BiocManager::install(c("edgeR", "limma", "EnsDb.Hsapiens.v86"), update = FALSE, ask = FALSE); install.packages("msigdbr")
