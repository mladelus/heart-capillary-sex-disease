# K2  EXPORT from the two KPMP source files already downloaded by K1 (no analysis is done here)
# For each dataset it writes, into ~/Documents/heart_capillary_disease/kidney/export/ :
#   <ds>_obs.rds          the full cell table
#   <ds>_genes.rds        gene identifiers and names
#   <ds>_pb_<label>.rds   raw counts summed per donor x cell label (pseudobulks), one file per label column
#   <ds>_EC_counts.rds    raw counts of every endothelial cell (genes x cells, sparse)
# The files are read in blocks, so memory use stays low. Finished datasets are skipped on a re-run.
# Run:  source("~/Desktop/heart-capillary-sex-disease/kidney/K2_export.R")

OUT <- path.expand("~/Documents/heart_capillary_disease/kidney")
EXP <- file.path(OUT, "export"); dir.create(EXP, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(rhdf5); library(Matrix) })
while (sink.number() > 0) sink()
sink(file.path(OUT, "tables", "K2_export_output.txt"), split = TRUE)
BLOCK <- 15000                                   # cells per block
LABELS <- c("subclass.l1", "subclass.l2", "cell_type")

read_table <- function(f, grp) {                 # reads /obs or /var (categoricals, nullable and plain columns)
  ls <- h5ls(f, recursive = TRUE); ls <- ls[ls$group == grp | startsWith(ls$group, paste0(grp, "/")), ]
  top <- ls[ls$group == grp, ]; out <- list()
  for (i in seq_len(nrow(top))) {
    nm <- top$name[i]; path <- paste0(grp, "/", nm)
    v <- tryCatch({
      if (top$otype[i] == "H5I_GROUP") {
        kids <- ls$name[ls$group == path]
        if (all(c("categories", "codes") %in% kids)) {
          cats <- as.character(h5read(f, paste0(path, "/categories"), bit64conversion = "double")); codes <- as.integer(h5read(f, paste0(path, "/codes")))
          ifelse(codes < 0, NA, cats[codes + 1])
        } else if ("values" %in% kids) {
          x <- as.vector(h5read(f, paste0(path, "/values"), bit64conversion = "double"))
          if ("mask" %in% kids) x[as.logical(h5read(f, paste0(path, "/mask")))] <- NA
          x
        } else NULL
      } else as.vector(h5read(f, path, bit64conversion = "double"))
    }, error = function(e) NULL)
    if (!is.null(v)) out[[nm]] <- v
  }
  n <- max(lengths(out)); out <- out[lengths(out) == n]
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

export_one <- function(ds) {
  f <- file.path(OUT, "source", paste0("KPMP_", ds, ".h5ad"))
  if (!file.exists(f)) stop("Source file not found (run K1 first): ", f)
  done <- file.path(EXP, paste0(ds, "_DONE.txt"))
  if (file.exists(done)) { cat("\n", ds, ": already exported, skipped\n"); return(invisible()) }
  cat("\n==================== ", ds, " ====================\n", sep = "")
  obs <- read_table(f, "/obs"); cat("cells:", nrow(obs), "| columns:", ncol(obs), "\n")
  saveRDS(obs, file.path(EXP, paste0(ds, "_obs.rds")))
  tree <- h5ls(f, recursive = TRUE)
  xg <- if (any(tree$group == "/raw/X")) "/raw/X" else "/X"; vg <- if (xg == "/raw/X") "/raw/var" else "/var"
  at <- h5readAttributes(f, xg, bit64conversion = "double"); enc <- as.character(at[["encoding-type"]]); shp <- as.numeric(at[["shape"]])
  cat("count matrix:", xg, "| encoding:", enc, "| shape:", shp[1], "x", shp[2], "\n")
  if (!identical(enc, "csr_matrix")) stop("The matrix is stored as '", enc, "', which this script does not read. Please send this message to Claude.")
  var <- read_table(f, vg); var0 <- if (vg != "/var") read_table(f, "/var") else var
  id <- if ("_index" %in% names(var)) var[["_index"]] else if ("feature_id" %in% names(var)) var$feature_id else var[[1]]
  id0 <- if ("_index" %in% names(var0)) var0[["_index"]] else id
  nm <- if ("feature_name" %in% names(var0)) as.character(var0$feature_name)[match(id, id0)] else id
  genes <- data.frame(feature_id = id, feature_name = nm, stringsAsFactors = FALSE)
  saveRDS(genes, file.path(EXP, paste0(ds, "_genes.rds")))
  ng <- nrow(genes); nc <- nrow(obs); stopifnot(shp[1] == nc, shp[2] == ng)
  cat("genes:", ng, "| with a name:", sum(!is.na(genes$feature_name)), "\n")

  ctl <- tolower(as.character(obs$cell_type))
  is_ec <- grepl("endothelial|vasa recta", ctl)
  if ("class" %in% names(obs)) is_ec <- is_ec | grepl("endothelial", tolower(as.character(obs$class)))
  if ("subclass.l1" %in% names(obs)) is_ec <- is_ec | as.character(obs$subclass.l1) == "EC"
  cat("endothelial cells to export:", sum(is_ec), "\n")
  labs <- intersect(LABELS, names(obs)); cat("label columns used for pseudobulks:", paste(labs, collapse = ", "), "\n")
  grp <- lapply(labs, function(l) { k <- paste(obs$donor_id, obs[[l]], sep = "||"); fct <- factor(k); list(f = fct, n = nlevels(fct)) }); names(grp) <- labs
  pb <- lapply(grp, function(g) Matrix(0, nrow = ng, ncol = g$n, sparse = TRUE))
  ecl <- list(); indptr <- as.numeric(h5read(f, paste0(xg, "/indptr"), bit64conversion = "double")); stopifnot(length(indptr) == nc + 1)
  starts <- seq(1, nc, by = BLOCK); checked <- FALSE; tot <- numeric(nc)
  for (b in seq_along(starts)) {
    r1 <- starts[b]; r2 <- min(nc, r1 + BLOCK - 1); p1 <- indptr[r1]; p2 <- indptr[r2 + 1]; n <- p2 - p1
    if (n == 0) next
    x <- as.numeric(h5read(f, paste0(xg, "/data"), start = p1 + 1, count = n, bit64conversion = "double"))
    j <- as.integer(h5read(f, paste0(xg, "/indices"), start = p1 + 1, count = n, bit64conversion = "double"))
    if (!checked) { cat("first block: values are whole numbers =", all(abs(x - round(x)) < 1e-6), "| max =", max(x), "\n"); checked <- TRUE
      if (!all(abs(x - round(x)) < 1e-6)) stop("The matrix does not hold raw counts. Please send this message to Claude.") }
    len <- diff(indptr[r1:(r2 + 1)]); i <- rep.int(seq_len(r2 - r1 + 1), len)
    M <- sparseMatrix(i = j + 1L, j = i, x = x, dims = c(ng, r2 - r1 + 1))            # genes x cells of this block
    tot[r1:r2] <- colSums(M)
    for (l in labs) { ind <- sparseMatrix(i = seq_len(r2 - r1 + 1), j = as.integer(grp[[l]]$f[r1:r2]), x = 1, dims = c(r2 - r1 + 1, grp[[l]]$n)); pb[[l]] <- pb[[l]] + M %*% ind }
    e <- which(is_ec[r1:r2]); if (length(e)) ecl[[length(ecl) + 1]] <- M[, e, drop = FALSE]
    cat(sprintf("  block %d of %d (cells %d-%d)\n", b, length(starts), r1, r2)); rm(M, x, j, i); gc(verbose = FALSE)
  }
  for (l in labs) { m <- pb[[l]]; rownames(m) <- genes$feature_id; colnames(m) <- levels(grp[[l]]$f)
    saveRDS(as(m, "CsparseMatrix"), file.path(EXP, paste0(ds, "_pb_", gsub("[^A-Za-z0-9]", "", l), ".rds"))); cat("pseudobulk", l, ":", ncol(m), "donor x label groups\n") }
  ec <- do.call(cbind, ecl); rownames(ec) <- genes$feature_id; colnames(ec) <- obs[["_index"]][is_ec]
  half <- ceiling(ncol(ec) / 2)                                                       # two files, to keep each one small
  saveRDS(as(ec[, seq_len(half), drop = FALSE], "CsparseMatrix"), file.path(EXP, paste0(ds, "_EC_counts_1.rds")))
  saveRDS(as(ec[, (half + 1):ncol(ec), drop = FALSE], "CsparseMatrix"), file.path(EXP, paste0(ds, "_EC_counts_2.rds")))
  obs$total_counts_raw <- tot; saveRDS(obs, file.path(EXP, paste0(ds, "_obs.rds")))
  cat("endothelial matrix:", nrow(ec), "genes x", ncol(ec), "cells | median counts per endothelial cell:", median(colSums(ec)), "\n")
  writeLines(format(Sys.time()), done)
}
for (ds in c("snRNA", "scRNA")) export_one(ds)
cat("\nFiles written:\n"); print(data.frame(file = list.files(EXP), MB = round(file.size(list.files(EXP, full.names = TRUE)) / 1e6, 1)))
sink()
message("\nFinished. Files are in ", EXP)
