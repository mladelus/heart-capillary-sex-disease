# K0  INVENTORY (metadata only, no expression values are read)
# What is in the CELLxGENE Census for (a) adult human kidney, including the KPMP atlas, and (b) Tabula Sapiens, all organs:
# donors by sex, endothelial cells per donor, endothelial labels, assays, disease labels.
# The analysis plan for the kidney / multi-organ endothelium work is written from this inventory, before any expression data
# are downloaded.
# Run in a fresh R window:  source("~/Desktop/heart-capillary-sex-disease/kidney/K0_inventory.R")

OUT <- path.expand("~/Documents/heart_capillary_disease/kidney")
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)
CENSUS_VERSION <- "2025-11-08"      # same pinned release as the heart work
MIN_EC <- 100                       # same donor rule as the heart work
suppressPackageStartupMessages({ library(cellxgene.census); library(dplyr); library(tidyr); library(stringr) })
options(width = 200, dplyr.summarise.inform = FALSE)
sink(file.path(OUT, "tables", "K0_inventory_output.txt"), split = TRUE)

census <- open_soma(census_version = CENSUS_VERSION)
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat()) |>
  select(dataset_id, dataset_title, collection_name, collection_doi, dataset_total_cell_count)
cols <- c("dataset_id", "donor_id", "sex", "development_stage", "tissue", "tissue_general", "disease", "cell_type",
          "assay", "suspension_type", "is_primary_data")
read_obs <- function(filter) {
  census$get("census_data")$get("homo_sapiens")$obs$read(value_filter = filter, column_names = cols)$concat() |>
    as.data.frame() |> as_tibble() |> mutate(across(where(is.factor), as.character))
}
dec_mid <- c("third decade stage" = 25, "fourth decade stage" = 35, "fifth decade stage" = 45, "sixth decade stage" = 55,
             "seventh decade stage" = 65, "eighth decade stage" = 75, "ninth decade stage" = 85)
parse_age <- function(x) coalesce(as.numeric(str_extract(x, "^\\d+(?=[- ]year)")), unname(dec_mid[x]))
is_blood_ec <- function(ct) {
  ct <- tolower(ct)
  (str_detect(ct, "endothelial|vasa recta|aerocyte")) & !str_detect(ct, "lymphatic|endocardial")
}
prep <- function(o) o |> left_join(ds, by = "dataset_id") |>
  mutate(age = parse_age(development_stage), blood_ec = is_blood_ec(cell_type))

# ---- A. kidney -------------------------------------------------------------------------------------------------------
message("Reading kidney metadata ...")
kid <- prep(read_obs("tissue_general == 'kidney' & is_primary_data == TRUE"))
cat("\n== A1. Kidney datasets in the Census (primary cells) ==\n")
kid |> group_by(collection_name, dataset_title, dataset_id) |>
  summarise(cells = n(), donors = n_distinct(donor_id), blood_ECs = sum(blood_ec),
            sexes = paste(sort(unique(sex)), collapse = "/"), assays = paste(sort(unique(assay)), collapse = "; ")) |>
  arrange(desc(cells)) |> mutate(dataset_title = substr(dataset_title, 1, 60), collection_name = substr(collection_name, 1, 50)) |>
  print(n = Inf, width = Inf)

kd <- kid |> filter(sex %in% c("female", "male")) |>
  group_by(collection_name, dataset_title, dataset_id, donor_id, sex) |>
  summarise(age = suppressWarnings(max(age, na.rm = TRUE)), stage = paste(sort(unique(development_stage)), collapse = "; "),
            disease = paste(sort(unique(disease)), collapse = "; "), tissues = paste(sort(unique(tissue)), collapse = "; "),
            assay = paste(sort(unique(assay)), collapse = "; "), suspension = paste(sort(unique(suspension_type)), collapse = "; "),
            cells = n(), blood_ECs = sum(blood_ec)) |> ungroup() |>
  mutate(age = ifelse(is.finite(age), age, NA), eligible = blood_ECs >= MIN_EC & (is.na(age) | age >= 18))
write.csv(kd, file.path(OUT, "tables", "K0_kidney_donors.csv"), row.names = FALSE)

cat("\n== A2. Kidney: donors with >= ", MIN_EC, " blood endothelial cells, by dataset, disease, assay and sex ==\n", sep = "")
kd |> filter(eligible) |> count(collection_name, dataset_title, disease, assay, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |>
  mutate(dataset_title = substr(dataset_title, 1, 45), collection_name = substr(collection_name, 1, 40)) |> print(n = Inf, width = Inf)

cat("\n== A3. Kidney: the same donor in more than one dataset? (donor_id seen in > 1 dataset) ==\n")
kd |> count(donor_id, name = "datasets") |> count(datasets, name = "donors") |> print()

cat("\n== A4. Kidney: age recorded for eligible donors ==\n")
kd |> filter(eligible) |> group_by(dataset_title = substr(dataset_title, 1, 45), sex) |>
  summarise(donors = n(), with_age = sum(!is.na(age)), median_age = median(age, na.rm = TRUE), min = min(age, na.rm = TRUE), max = max(age, na.rm = TRUE),
            example_stage = first(stage)) |> print(n = Inf, width = Inf)

cat("\n== A5. Kidney: endothelial and related cell-type labels (cells; donors with >= 30 such cells, by sex) ==\n")
lab <- kid |> filter(str_detect(tolower(cell_type), "endothelial|vasa recta|capillary|glomerul|podocyte|mesangial|pericyte|smooth muscle"))
lab_d <- lab |> filter(sex %in% c("female", "male")) |> count(dataset_id, donor_id, sex, cell_type) |> filter(n >= 30) |>
  count(cell_type, sex, name = "donors") |> pivot_wider(names_from = sex, values_from = donors, values_fill = 0, names_prefix = "donors30_")
lab |> count(cell_type, name = "cells") |> left_join(lab_d, by = "cell_type") |> arrange(desc(cells)) |> print(n = Inf, width = Inf)
write.csv(kid |> count(collection_name, dataset_title, dataset_id, cell_type, sex, disease, assay, name = "cells"),
          file.path(OUT, "tables", "K0_kidney_celltypes.csv"), row.names = FALSE)

# ---- B. Tabula Sapiens, all organs ------------------------------------------------------------------------------------
ts_ds <- ds |> filter(str_detect(tolower(collection_name), "tabula sapiens"))
cat("\n== B1. Tabula Sapiens datasets in the Census ==\n")
ts_ds |> mutate(dataset_title = substr(dataset_title, 1, 70)) |> select(collection_name, dataset_title, dataset_total_cell_count) |> print(width = Inf)
message("Reading Tabula Sapiens metadata ...")
ts <- bind_rows(lapply(ts_ds$dataset_id, function(id) read_obs(paste0("dataset_id == '", id, "' & is_primary_data == TRUE")))) |> prep()

td <- ts |> filter(sex %in% c("female", "male")) |>
  group_by(collection_name, donor_id, sex, tissue_general) |>
  summarise(age = suppressWarnings(max(age, na.rm = TRUE)), stage = first(development_stage),
            assay = paste(sort(unique(assay)), collapse = "; "), cells = n(), blood_ECs = sum(blood_ec)) |> ungroup() |>
  mutate(age = ifelse(is.finite(age), age, NA))
write.csv(td, file.path(OUT, "tables", "K0_tabula_sapiens_donor_by_organ.csv"), row.names = FALSE)

cat("\n== B2. Tabula Sapiens donors (sex, age, assays, organs, blood endothelial cells in total) ==\n")
td |> group_by(donor_id, sex, age, stage) |>
  summarise(organs = n_distinct(tissue_general), cells = sum(cells), blood_ECs = sum(blood_ECs), assay = paste(sort(unique(unlist(strsplit(assay, "; ")))), collapse = "; ")) |>
  arrange(sex, donor_id) |> print(n = Inf, width = Inf)

cat("\n== B3. Tabula Sapiens: donors per organ by sex, at three endothelial-cell thresholds ==\n")
td |> group_by(tissue_general) |>
  summarise(blood_ECs = sum(blood_ECs),
            F_any = sum(sex == "female" & blood_ECs > 0), M_any = sum(sex == "male" & blood_ECs > 0),
            F_30 = sum(sex == "female" & blood_ECs >= 30), M_30 = sum(sex == "male" & blood_ECs >= 30),
            F_100 = sum(sex == "female" & blood_ECs >= 100), M_100 = sum(sex == "male" & blood_ECs >= 100)) |>
  arrange(desc(blood_ECs)) |> print(n = Inf, width = Inf)

cat("\n== B4. Tabula Sapiens: endothelial cell-type labels ==\n")
ts |> filter(blood_ec | str_detect(tolower(cell_type), "lymphatic|endocardial")) |> count(cell_type, name = "cells") |> arrange(desc(cells)) |> print(n = Inf, width = Inf)
write.csv(ts |> filter(blood_ec) |> count(donor_id, sex, tissue_general, tissue, cell_type, assay, name = "cells"),
          file.path(OUT, "tables", "K0_tabula_sapiens_ec_labels.csv"), row.names = FALSE)

cat("\n== B5. Tabula Sapiens kidney donors also present in the kidney datasets above? ==\n")
print(intersect(unique(ts$donor_id[ts$tissue_general == "kidney"]), unique(kid$donor_id[!kid$dataset_id %in% ts_ds$dataset_id])))

sink()
message("\nFinished. Tables are in ", file.path(OUT, "tables"))
