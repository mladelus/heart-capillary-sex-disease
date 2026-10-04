# Plan section O: MAGMA gene-set test of the endothelial state against full GWAS summary statistics.
# Run in R on a Mac:  source("~/Desktop/heart-capillary-sex-disease/genetics/O1_magma.R")
# The script can be re-run: finished downloads and finished steps are skipped.

repo <- path.expand("~/Desktop/heart-capillary-sex-disease")
proj <- path.expand("~/Documents/heart_capillary_disease")
base <- file.path(proj, "gwas", "magma")
dir.create(base, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(proj, "tables"), showWarnings = FALSE)
options(timeout = 6 * 3600)
for (p in c("data.table", "R.utils")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table)

# ---- 1. downloads ---------------------------------------------------------------------------------------------------
ebi <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
vu <- function(id) paste0("https://vu.data.surfsara.nl/index.php/s/", id, "/download")
files <- list(
  magma_mac.zip = vu("1M1d9vHtVidEwvU"),
  NCBI37.3.zip = vu("Pj2orwuF2JYyKxq"),
  g1000_eur.zip = vu("VZNByNwpD8qqINe"),
  cad.tsv.gz = paste0(ebi, "GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz"),
  scad.tsv.gz = paste0(ebi, "GCST90245001-GCST90246000/GCST90245878/GCST90245878_buildGRCh37.tsv.gz"),
  hf.tsv.gz = paste0(ebi, "GCST009001-GCST010000/GCST009541/ShahS_31919418_HeartFailure.gz"),
  sbp.txt.gz = paste0(ebi, "GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz"),
  height.txt.gz = paste0(ebi, "GCST006001-GCST007000/GCST006901/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz"))
for (f in names(files)) {
  dest <- file.path(base, f)
  if (file.exists(dest) && file.size(dest) > 1000) next
  message("Downloading ", f, " ...")
  tmp <- paste0(dest, ".part")
  if (download.file(files[[f]], tmp, mode = "wb", method = "libcurl") != 0) stop("Download failed: ", f)
  file.rename(tmp, dest)
}
for (z in c("magma_mac.zip", "NCBI37.3.zip", "g1000_eur.zip")) {
  d <- file.path(base, sub(".zip", "", z, fixed = TRUE))
  if (!dir.exists(d)) unzip(file.path(base, z), exdir = d)
}
find1 <- function(pattern) {
  x <- list.files(base, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(x) < 1) stop("File not found after unzip: ", pattern)
  x[1]
}
magma <- find1("^magma$"); Sys.chmod(magma, "755")
bim <- find1("^g1000_eur\\.bim$"); bfile <- sub("\\.bim$", "", bim)
geneloc_raw <- find1("^NCBI37\\.3\\.gene\\.loc$")
run <- function(args, tag) {
  message("MAGMA: ", tag)
  st <- system2(magma, args, stdout = file.path(base, paste0("log_", tag, ".txt")), stderr = file.path(base, paste0("err_", tag, ".txt")))
  if (st != 0) stop("MAGMA step failed: ", tag, ". See log_", tag, ".txt and err_", tag, ".txt in ", base,
                    "\nIf the program cannot start on an Apple Silicon Mac, run this once in Terminal: softwareupdate --install-rosetta")
}

# ---- 2. gene locations (build 37, autosomes, extended MHC removed) and annotation --------------------------------
gl <- fread(geneloc_raw, header = FALSE, col.names = c("entrez", "chr", "start", "stop", "strand", "symbol"))
gl <- gl[chr %in% as.character(1:22)]
gl <- gl[!(duplicated(symbol) | duplicated(symbol, fromLast = TRUE))]
gl <- gl[!(chr == "6" & stop >= 25e6 & start <= 34e6)]
loc <- file.path(base, "genes_b37_symbol.loc")
fwrite(gl[, .(symbol, chr, start, stop, strand)], loc, sep = "\t", col.names = FALSE)
annot <- file.path(base, "annot")
if (!file.exists(paste0(annot, ".genes.annot")))
  run(c("--annotate", "window=35,10", "--snp-loc", bim, "--gene-loc", loc, "--out", annot), "annotate")

# ---- 3. p-value files: SNP, P, N -------------------------------------------------------------------------------------
ref <- fread(bim, header = FALSE, select = c(1, 2, 4), col.names = c("chr", "SNP", "pos"))
clean <- function(d) {
  d <- d[!is.na(SNP) & grepl("^rs", SNP) & !is.na(P) & P >= 0 & P <= 1 & !is.na(N) & N > 0]
  d[P < 1e-300, P := 1e-300]
  d[, N := as.integer(round(N))]
  d[!duplicated(SNP)]
}
readers <- list(
  cad = function(f) { d <- fread(f, select = c("rsid", "p_value", "n")); d[, .(SNP = rsid, P = p_value, N = n)] },
  scad = function(f) { d <- fread(f, select = c("variant_id", "p_value", "n_cases", "n_ctrls")); d[, .(SNP = variant_id, P = p_value, N = n_cases + n_ctrls)] },
  hf = function(f) { d <- fread(f, select = c("SNP", "p", "N")); d[, .(SNP, P = p, N)] },
  sbp = function(f) {
    d <- fread(f, select = c("MarkerName", "P", "TotalSampleSize"))
    d <- d[grepl(":SNP$", MarkerName)]
    d[, c("chr", "pos") := tstrsplit(MarkerName, ":", fixed = TRUE, keep = 1:2)]
    d[, `:=`(chr = as.character(chr), pos = as.integer(pos))]
    r <- copy(ref); r[, chr := as.character(chr)]
    d <- merge(d, r, by = c("chr", "pos"))
    d[, .(SNP, P, N = TotalSampleSize)]
  },
  height = function(f) { d <- fread(f, select = c("SNP", "P", "N")); d[, .(SNP, P = as.numeric(P), N)] })
src <- c(cad = "cad.tsv.gz", scad = "scad.tsv.gz", hf = "hf.tsv.gz", sbp = "sbp.txt.gz", height = "height.txt.gz")
acc <- c(cad = "GCST90132314", scad = "GCST90245878", hf = "GCST009541", sbp = "GCST006624", height = "GCST006901")
step0 <- list()
for (t in names(src)) {
  out <- file.path(base, paste0(t, ".pval.txt"))
  if (!file.exists(out)) {
    message("Preparing ", t, " ...")
    d <- clean(readers[[t]](file.path(base, src[[t]])))
    d[, P := as.numeric(P)]
    fwrite(d, out, sep = "\t", scipen = 0)
    rm(d); gc()
  }
  d <- fread(out, select = c("SNP", "N"))
  step0[[t]] <- data.table(trait = t, accession = acc[[t]], file = src[[t]], variants_with_rsid = nrow(d),
                           variants_in_reference = sum(d$SNP %in% ref$SNP), median_N = median(d$N), max_N = max(d$N))
  rm(d); gc()
}
step0 <- rbindlist(step0); print(step0)
fwrite(step0, file.path(proj, "tables", "S3O_step0_files.csv"))

# ---- 4. gene analysis (all autosomal genes; SNP-wise mean) --------------------------------------------------------
for (t in names(src)) {
  out <- file.path(base, t)
  if (file.exists(paste0(out, ".genes.raw"))) next
  run(c("--bfile", bfile, "--gene-annot", paste0(annot, ".genes.annot"), "--pval", file.path(base, paste0(t, ".pval.txt")),
        "use=SNP,P", "ncol=N", "--gene-model", "snp-wise=mean", "--out", out), paste0("genes_", t))
}

# ---- 5. gene sets, universe and covariate (frozen list from plan section M) ---------------------------------------
sg <- fread(file.path(repo, "genetics", "state_genes.csv"))
sg <- sg[gene %in% gl$symbol]
locus <- gl[chr == "6" & stop >= 12903957 - 1e6 & start <= 12903957 + 1e6, symbol]   # rs9349379, GRCh37
sets <- list(state_all = sg[state == 1, gene], state_up = sg[state_up == 1, gene], state_down = sg[state_down == 1, gene],
             state_all_without_EDN1_locus = setdiff(sg[state == 1, gene], locus))
writeLines(vapply(names(sets), function(s) paste(c(s, sets[[s]]), collapse = " "), ""), file.path(base, "sets.txt"))
writeLines(sg$gene, file.path(base, "universe.txt"))
fwrite(sg[, .(GENE = gene, mean_log2cpm)], file.path(base, "covar.txt"), sep = "\t")
message("Universe genes with a build-37 location: ", nrow(sg), "; state genes: ", length(sets$state_all),
        " (up ", length(sets$state_up), ", down ", length(sets$state_down), "); state genes within 1 Mb of rs9349379 removed in the sensitivity set: ",
        length(sets$state_all) - length(sets$state_all_without_EDN1_locus))

# ---- 6. competitive gene-set analysis --------------------------------------------------------------------------------
S <- file.path(base, "sets.txt"); U <- paste0("gene-include=", file.path(base, "universe.txt")); C <- file.path(base, "covar.txt")
for (t in names(src)) {
  raw <- paste0(file.path(base, t), ".genes.raw")
  run(c("--gene-results", raw, "--set-annot", S, "--gene-covar", C, "--model", "condition-hide=mean_log2cpm", "--settings", U,
        "--out", file.path(base, paste0(t, ".A_primary"))), paste0("sets_A_", t))
  run(c("--gene-results", raw, "--set-annot", S, "--settings", U, "--out", file.path(base, paste0(t, ".B_no_expression"))), paste0("sets_B_", t))
  run(c("--gene-results", raw, "--set-annot", S, "--out", file.path(base, paste0(t, ".C_all_genes"))), paste0("sets_C_", t))
}

# ---- 7. collect -------------------------------------------------------------------------------------------------------
res <- list()
for (f in list.files(base, pattern = "\\.gsa\\.out$", full.names = TRUE)) {
  x <- readLines(f); x <- x[!grepl("^#", x)]
  d <- fread(text = x)
  nm <- strsplit(sub("\\.gsa\\.out$", "", basename(f)), ".", fixed = TRUE)[[1]]
  lg <- readLines(file.path(base, paste0("log_sets_", substr(nm[2], 1, 1), "_", nm[1], ".txt")))
  ng <- grep("genes", lg, value = TRUE, ignore.case = TRUE)
  res[[f]] <- cbind(data.table(trait = nm[1], analysis = nm[2]), d, log_gene_lines = paste(trimws(ng), collapse = " || "))
}
res <- rbindlist(res, fill = TRUE)
labs <- c(cad = "Coronary artery disease", scad = "Spontaneous coronary artery dissection", hf = "Heart failure",
          sbp = "Systolic blood pressure", height = "CONTROL height")
res[, trait_label := labs[trait]]
setorder(res, analysis, trait, VARIABLE)
fwrite(res, file.path(proj, "tables", "S3O_magma_gene_sets.csv"))
print(res[, .(analysis, trait_label, VARIABLE, NGENES, BETA, SE, P)])
message("Finished. Results: ", file.path(proj, "tables", "S3O_magma_gene_sets.csv"))
