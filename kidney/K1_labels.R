# K1  AUTHOR LABELS (metadata only, no expression values are read)
# The Census labels almost all kidney endothelial cells simply "endothelial cell". The atlas authors' own labels
# (peritubular, glomerular, arteriolar, vasa recta; injury states) are in the source files of the two KPMP datasets.
# This script downloads those two files, extracts only the cell table (obs), and summarises the labels by donor and sex.
# It also finishes the Tabula Sapiens inventory that stopped with a printing error in K0.
# Run:  source("~/Desktop/heart-capillary-sex-disease/kidney/K1_labels.R")
# The two downloads are large (several GB in total); finished downloads are skipped on a re-run.

OUT <- path.expand("~/Documents/heart_capillary_disease/kidney")
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "source"), showWarnings = FALSE)
CENSUS_VERSION <- "2025-11-08"
options(timeout = 6 * 3600, width = 200, dplyr.summarise.inform = FALSE)
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("rhdf5", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({ library(cellxgene.census); library(rhdf5); library(dplyr); library(tidyr); library(stringr) })
while (sink.number() > 0) sink()          # close any output file left open by an earlier run
sink(file.path(OUT, "tables", "K1_labels_output_kidney.txt"), split = TRUE)
show <- function(x, n = Inf) print(as_tibble(x), n = n, width = Inf)

DS <- c(snRNA = "a12ccb9b-4fbe-457d-8590-ac78053259ef", scRNA = "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d")

# ---- read the obs table of an h5ad file without loading expression ---------------------------------------------------
read_obs_h5ad <- function(f) {
  ls <- h5ls(f, recursive = 3); ls <- ls[grepl("^/obs", ls$group), ]; top <- ls[ls$group == "/obs", ]
  if (!nrow(top)) stop("This file stores obs in an old format that this script does not read: ", f)
  out <- list()
  for (i in seq_len(nrow(top))) {
    nm <- top$name[i]; path <- paste0("/obs/", nm)
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
  as_tibble(as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE))
}

obs_all <- list()
for (k in names(DS)) {
  f <- file.path(OUT, "source", paste0("KPMP_", k, ".h5ad"))
  if (!file.exists(f) || file.size(f) < 1e8) {
    message("Downloading the ", k, " source file (large) ...")
    tmp <- paste0(f, ".part")
    download_source_h5ad(DS[[k]], tmp, overwrite = TRUE, census_version = CENSUS_VERSION)
    file.rename(tmp, f)
  }
  message("Reading the cell table of ", k, " ...")
  o <- read_obs_h5ad(f); o$dataset <- k
  saveRDS(o, file.path(OUT, "source", paste0("KPMP_", k, "_obs.rds")))
  obs_all[[k]] <- o
  cat("\n==================== ", k, ": ", nrow(o), " cells, ", ncol(o), " columns ====================\n", sep = "")
  cat("\n-- Columns, with number of distinct values and examples --\n")
  show(tibble(column = names(o), type = vapply(o, function(x) class(x)[1], ""), distinct = vapply(o, n_distinct, 1L),
              examples = vapply(o, function(x) paste(head(unique(as.character(x)), 6), collapse = " | "), "")) |>
         mutate(examples = substr(examples, 1, 110)))
  # tabulate every label-like column with few levels (cell types, states, conditions, clinical categories)
  lab <- names(o)[vapply(o, function(x) (is.character(x) || is.logical(x)) && n_distinct(x) <= 120, TRUE)]
  lab <- setdiff(lab, c("dataset"))
  for (cl in lab) {
    cat("\n-- ", k, " | ", cl, " (cells; donors with >= 30 cells, by sex) --\n", sep = "")
    if (all(c("donor_id", "sex") %in% names(o))) {
      d <- o |> count(.data[[cl]], donor_id, sex, name = "n") |> filter(n >= 30) |> count(.data[[cl]], sex, name = "donors") |>
        pivot_wider(names_from = sex, values_from = donors, values_fill = 0, names_prefix = "donors30_")
      show(o |> count(.data[[cl]], name = "cells") |> left_join(d, by = cl) |> arrange(desc(cells)), n = 130)
    } else show(o |> count(.data[[cl]], name = "cells") |> arrange(desc(cells)), n = 130)
  }
  num <- names(o)[vapply(o, is.numeric, TRUE)]
  cat("\n-- ", k, " | numeric columns (range) --\n", sep = "")
  show(tibble(column = num, min = vapply(o[num], function(x) as.numeric(suppressWarnings(min(x, na.rm = TRUE))), 1), median = vapply(o[num], function(x) as.numeric(median(x, na.rm = TRUE)), 1),
              max = vapply(o[num], function(x) as.numeric(suppressWarnings(max(x, na.rm = TRUE))), 1), missing = vapply(o[num], function(x) sum(is.na(x)), 1L)))
}

# ---- donors: overlap between the two datasets, and donor-level table ---------------------------------------------------
if (all(vapply(obs_all, function(o) all(c("donor_id", "sex", "disease") %in% names(o)), TRUE))) {
  dd <- bind_rows(lapply(obs_all, function(o) o |> group_by(dataset, donor_id, sex, disease) |> summarise(cells = n()) |> ungroup()))
  cat("\n==================== Donors ====================\n")
  cat("donors in snRNA:", n_distinct(dd$donor_id[dd$dataset == "snRNA"]), "| in scRNA:", n_distinct(dd$donor_id[dd$dataset == "scRNA"]),
      "| in both:", length(intersect(dd$donor_id[dd$dataset == "snRNA"], dd$donor_id[dd$dataset == "scRNA"])), "\n")
  both <- intersect(dd$donor_id[dd$dataset == "snRNA"], dd$donor_id[dd$dataset == "scRNA"])
  cat("\n-- scRNA donors NOT in snRNA, by disease and sex --\n")
  show(dd |> filter(dataset == "scRNA", !donor_id %in% both) |> count(disease, sex) |> pivot_wider(names_from = sex, values_from = n, values_fill = 0))
  write.csv(dd, file.path(OUT, "tables", "K1_KPMP_donors.csv"), row.names = FALSE)
}

# (The Tabula Sapiens inventory was completed in the first run and is not repeated.)
sink()
message("\nFinished. Output saved in ", file.path(OUT, "tables", "K1_labels_output_kidney.txt"))
