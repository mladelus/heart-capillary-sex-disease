# Runs the whole paper 2 analysis in order. Set the working directory to this folder first:
#   setwd("~/Desktop/heart-capillary-sex-disease"); source("run_all.R")
# Needs the paper 1 outputs (pb.rds, aim1.rds, aim4.rds) in ~/Documents/heart_capillary_project
# (or set HCAP_DATA). Paper 2 data go to ~/Documents/heart_capillary_disease (or set HCAP2_DATA).
repo <- getwd()
for (f in c("00_setup.R", "01_inventory.R", "02_download.R", "03_ec_subtypes_pseudobulk.R",
            "04_primary_sex_within_disease.R", "05_sex_by_disease.R", "06_pericytes_xy.R", "07_figures.R")) {
  message("\n########## ", f, " ##########")
  source(file.path(repo, "R", f))
}
