# =============================================================================
# 02b  DOWNLOAD MURAL CELLS (added before any pericyte result was seen; see Deviation 1 in ANALYSIS_PLAN.md)
# Two of the three datasets label pericytes and smooth muscle cells together as "mural cell", so script 02
# put them in "other". Here every dataset gets the same rule: mural = pericyte + mural cell + smooth muscle.
# Downloads only those cells, in chunks of about 20,000 cells, and makes one pseudobulk per donor.
# Skips chunks already saved, so if it stops just run it again.  Makes: dl_mural/chunk_XX.rds
# =============================================================================
MURAL <- "^(pericyte|mural cell|smooth muscle.*)$"
donors <- readRDS("inv_donors.rds") |> filter(included)
mc <- readRDS("inv_cells.rds") |> filter(str_detect(cell_type, MURAL)) |>
  semi_join(donors, by = c("dataset_id", "donor_id")) |> arrange(dataset_id, donor_id)
cat("Mural cells to download:", nrow(mc), "from", n_distinct(paste(mc$dataset_id, mc$donor_id)), "donors\n")
print(mc |> count(dataset = substr(dataset_title, 1, 30), cell_type))

keys <- mc |> count(dataset_id, donor_id) |> mutate(chunk = 1 + cumsum(n) %/% 20000)
dir.create("dl_mural", showWarnings = FALSE)
for (ch in sort(unique(keys$chunk))) {
  out <- file.path("dl_mural", sprintf("chunk_%02d.rds", ch))
  if (file.exists(out)) next
  m <- mc |> semi_join(keys |> filter(chunk == ch), by = c("dataset_id", "donor_id"))
  message(format(Sys.time(), "%H:%M"), "  chunk ", ch, " of ", max(keys$chunk), "  (", nrow(m), " cells)")
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = m$soma_joinid, obs_column_names = "soma_joinid")
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  m <- m[match(seu$soma_joinid, m$soma_joinid), ]; rm(seu); gc(verbose = FALSE)
  grp <- split(seq_len(ncol(counts)), paste(substr(m$dataset_id, 1, 8), m$donor_id, sep = "|"))
  saveRDS(list(pb = lapply(grp, function(j) Matrix::rowSums(counts[, j, drop = FALSE])),
               n  = sapply(grp, length),
               n_pericyte_label = tapply(m$cell_type == "pericyte", paste(substr(m$dataset_id, 1, 8), m$donor_id, sep = "|"), sum)),
          out)
  rm(counts); gc(verbose = FALSE)
}
message("Saved mural chunks: ", length(list.files("dl_mural")), " of ", max(keys$chunk))
