###source
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
load("other_files/HMDB/MS1/hmdb_ms1.rda")
load("other_files/KEGG/kegg_ms1.rda")
source("R/read_msp_data.R")
source("R/3_utils.R")
source("R/7_GNPS.R")
setwd("other_files/GNPS")

library(tidyverse)
library(tidymass)
library(plyr)

# ###BILELIB19
# data <-
#   read_msp_data(file = "BILELIB19.msp",
#                 source = "gnps",
#                 threads = 5)
#
# gnps_bilelib19_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "BILELIB19",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_bilelib19_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_bilelib19_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_bilelib19_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_bilelib19_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_bilelib19_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_bilelib19_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_bilelib19_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
#
# gnps_bilelib19_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_bilelib19_ms2, file = "gnps_bilelib19_ms2.rda")
#
# ###IOBA-NHC
# data <-
#   read_msp_data(file = "GNPS-IOBA-NHC.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data)
#
# gnps_iobanhc_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "IOBA-NHC",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_iobanhc_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_iobanhc_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_iobanhc_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_iobanhc_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_iobanhc_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_iobanhc_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_iobanhc_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
#
# gnps_iobanhc_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_iobanhc_ms2, file = "gnps_iobanhc_ms2.rda")
#
#
#
# ###MSMLS
# data <-
#   read_msp_data(file = "GNPS-MSMLS.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data)
#
# gnps_msmls_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "MSMLS",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_msmls_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_msmls_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_msmls_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_msmls_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_msmls_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_msmls_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_msmls_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_msmls_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_msmls_ms2, file = "gnps_msmls_ms2.rda")
#
#
# ###NIH-CLINICALCOLLECTION1
# data <-
#   read_msp_data(file = "GNPS-NIH-CLINICALCOLLECTION1.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data)
#
# gnps_nihclinicalcollection1_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "NIH-CLINICALCOLLECTION1",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_nihclinicalcollection1_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_nihclinicalcollection1_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_nihclinicalcollection1_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_nihclinicalcollection1_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_nihclinicalcollection1_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_nihclinicalcollection1_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_nihclinicalcollection1_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_nihclinicalcollection1_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_nihclinicalcollection1_ms2, file = "gnps_nihclinicalcollection1_ms2.rda")
#
#
# ###NIH-CLINICALCOLLECTION2
# data <-
#   read_msp_data(file = "GNPS-NIH-CLINICALCOLLECTION2.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data)
#
# gnps_nihclinicalcollection2_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "NIH-CLINICALCOLLECTION2",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_nihclinicalcollection2_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_nihclinicalcollection2_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_nihclinicalcollection2_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_nihclinicalcollection2_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_nihclinicalcollection2_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_nihclinicalcollection2_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_nihclinicalcollection2_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_nihclinicalcollection2_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_nihclinicalcollection2_ms2, file = "gnps_nihclinicalcollection2_ms2.rda")
#
#
# ###NIH-NUTRI-METAB-FEM
# data1 <-
#   read_msp_data(file = "GNPS-NUTRI-METAB-FEM-POS.msp",
#                 source = "gnps",
#                 threads = 5)
#
# data2 <-
#   read_msp_data(file = "GNPS-NUTRI-METAB-FEM-NEG.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data1)
# length(data2)
#
# data <- c(data1, data2)
#
# gnps_nutrimetabfem_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "NUTRI-METAB-FEM",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_nutrimetabfem_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_nutrimetabfem_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_nutrimetabfem_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_nutrimetabfem_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_nutrimetabfem_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_nutrimetabfem_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_nutrimetabfem_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_nutrimetabfem_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_nutrimetabfem_ms2, file = "gnps_nutrimetabfem_ms2.rda")
#
#
#
# ###SCIEX-LIBRARY
# data <-
#   read_msp_data(file = "GNPS-SCIEX-LIBRARY.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data)
#
# gnps_sciexlibrary_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "SCIEX-LIBRARY",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_sciexlibrary_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_sciexlibrary_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_sciexlibrary_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_sciexlibrary_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_sciexlibrary_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_sciexlibrary_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_sciexlibrary_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_sciexlibrary_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_sciexlibrary_ms2, file = "gnps_sciexlibrary_ms2.rda")
#
#
# ###SELLECKCHEM-FDA-PART1
# data <-
#   read_msp_data(file = "GNPS-SELLECKCHEM-FDA-PART1.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data)
#
# gnps_selleckchemfdapart1_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "SELLECKCHEM-FDA-PART1",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_selleckchemfdapart1_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_selleckchemfdapart1_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_selleckchemfdapart1_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_selleckchemfdapart1_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_selleckchemfdapart1_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_selleckchemfdapart1_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_selleckchemfdapart1_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_selleckchemfdapart1_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_selleckchemfdapart1_ms2, file = "gnps_selleckchemfdapart1_ms2.rda")
#
#
# ###HCE-CELL-LYSATE-LIPIDS
# data <-
#   read_msp_data(file = "HCE-CELL-LYSATE-LIPIDS.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data)
#
# gnps_hcecelllysatelipids_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "HCE-CELL-LYSATE-LIPIDS",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_hcecelllysatelipids_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_hcecelllysatelipids_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_hcecelllysatelipids_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_hcecelllysatelipids_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_hcecelllysatelipids_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_hcecelllysatelipids_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_hcecelllysatelipids_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_hcecelllysatelipids_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_hcecelllysatelipids_ms2, file = "gnps_hcecelllysatelipids_ms2.rda")
#
#
# ###HCE-CELL-LYSATE-LIPIDS
# data1 <-
#   read_msp_data(file = "PNNL-LIPIDS-POSITIVE.msp",
#                 source = "gnps",
#                 threads = 5)
#
# data2 <-
#   read_msp_data(file = "PNNL-LIPIDS-NEGATIVE.msp",
#                 source = "gnps",
#                 threads = 5)
#
# length(data1)
# length(data2)
#
# data <- c(data1, data2)
#
# gnps_pnnllipids_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "PNNL-LIPIDS-POSITIVE",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_pnnllipids_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_pnnllipids_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_pnnllipids_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_pnnllipids_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_pnnllipids_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_pnnllipids_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_pnnllipids_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_pnnllipids_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_pnnllipids_ms2, file = "gnps_pnnllipids_ms2.rda")
#
#
# ###PSU-MSMLS
# data <-
#   read_msp_data(file = "PSU-MSMLS.msp",
#                 source = "gnps",
#                 threads = 5)
# length(data)
#
# gnps_pusmsmls_ms2 <-
#   convert_gnps2metid(data = data,
#                      path = "PSU-MSMLS",
#                      threads = 5)
#
# idx1 <-
#   match(
#     gnps_pusmsmls_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     gnps_pusmsmls_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     gnps_pusmsmls_ms2@spectra.info$CAS.ID,
#     hmdb_ms1@spectra.info$CAS.ID,
#     incomparables = NA
#   )
#
# idx4 <-
#   match(
#     gnps_pusmsmls_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3, idx4) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# cbind(
#   x = gnps_pusmsmls_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(gnps_pusmsmls_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   gnps_pusmsmls_ms2@spectra.info
#
# spectra.info <-
#   spectra.info %>%
#   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID, Formula))
#
# spectra.info <-
#   data.frame(spectra.info,
#              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
#
# gnps_pusmsmls_ms2@spectra.info <-
#   spectra.info
#
# save(gnps_pusmsmls_ms2, file = "gnps_pusmsmls_ms2.rda")
#

####combine all the gnps
load("gnps_bilelib19_ms2.rda")
load("gnps_hcecelllysatelipids_ms2.rda")
load("gnps_iobanhc_ms2.rda")
load("gnps_msmls_ms2.rda")
load("gnps_nihclinicalcollection1_ms2.rda")
load("gnps_nihclinicalcollection2_ms2.rda")
load("gnps_nutrimetabfem_ms2.rda")
load("gnps_pnnllipids_ms2.rda")
load("gnps_pusmsmls_ms2.rda")
load("gnps_sciexlibrary_ms2.rda")
load("gnps_selleckchemfdapart1_ms2.rda")

colnames(gnps_bilelib19_ms2@spectra.info)
colnames(gnps_hcecelllysatelipids_ms2@spectra.info)
colnames(gnps_iobanhc_ms2@spectra.info)
colnames(gnps_msmls_ms2@spectra.info)
colnames(gnps_nihclinicalcollection1_ms2@spectra.info)
colnames(gnps_nihclinicalcollection2_ms2@spectra.info)
colnames(gnps_nutrimetabfem_ms2@spectra.info)
colnames(gnps_pnnllipids_ms2@spectra.info)
colnames(gnps_pusmsmls_ms2@spectra.info)
colnames(gnps_sciexlibrary_ms2@spectra.info)
colnames(gnps_selleckchemfdapart1_ms2@spectra.info)

gnps_bilelib19_ms2@spectra.info$Submitter_team <- "BILELIB19"
gnps_hcecelllysatelipids_ms2@spectra.info$Submitter_team <- "HCE-CELL-LYSATE-LIPIDS"
gnps_iobanhc_ms2@spectra.info$Submitter_team <- "IOBA-NHC"
gnps_msmls_ms2@spectra.info$Submitter_team <- "MSMLS"
gnps_nihclinicalcollection1_ms2@spectra.info$Submitter_team <- "NIH-CLINICALCOLLECTION1"
gnps_nihclinicalcollection2_ms2@spectra.info$Submitter_team <- "NIH-CLINICALCOLLECTION2"
gnps_nutrimetabfem_ms2@spectra.info$Submitter_team <- "NUTRI-METAB-FEM"
gnps_pnnllipids_ms2@spectra.info$Submitter_team <- "PNNL-LIPIDS-POSITIVE"
gnps_pusmsmls_ms2@spectra.info$Submitter_team <- "PSU-MSMLS"
gnps_sciexlibrary_ms2@spectra.info$Submitter_team <- "SCIEX-LIBRARY"
gnps_selleckchemfdapart1_ms2@spectra.info$Submitter_team <- "SELLECKCHEM-FDA-PART1"

gnps_ms2 <- gnps_bilelib19_ms2

intersect(gnps_bilelib19_ms2@spectra.info$Lab.ID,gnps_hcecelllysatelipids_ms2@spectra.info$Lab.ID)
intersect(gnps_bilelib19_ms2@spectra.info$Lab.ID,gnps_iobanhc_ms2@spectra.info$Lab.ID)
intersect(gnps_bilelib19_ms2@spectra.info$Lab.ID,gnps_iobanhc_ms2@spectra.info$Lab.ID)
intersect(gnps_bilelib19_ms2@spectra.info$Lab.ID,gnps_msmls_ms2@spectra.info$Lab.ID)
intersect(gnps_bilelib19_ms2@spectra.info$Lab.ID,gnps_nihclinicalcollection1_ms2@spectra.info$Lab.ID)
intersect(gnps_bilelib19_ms2@spectra.info$Lab.ID,gnps_nihclinicalcollection1_ms2@spectra.info$Lab.ID)

gnps_ms2@spectra.info <-
  rbind(
    gnps_bilelib19_ms2@spectra.info,
    gnps_hcecelllysatelipids_ms2@spectra.info,
    gnps_iobanhc_ms2@spectra.info,
    gnps_msmls_ms2@spectra.info,
    gnps_nihclinicalcollection1_ms2@spectra.info,
    gnps_nihclinicalcollection2_ms2@spectra.info,
    gnps_nutrimetabfem_ms2@spectra.info,
    gnps_pnnllipids_ms2@spectra.info,
    gnps_pusmsmls_ms2@spectra.info,
    gnps_sciexlibrary_ms2@spectra.info,
    gnps_selleckchemfdapart1_ms2@spectra.info
  )

gnps_ms2@spectra.data$Spectra.positive <-
  c(
    gnps_bilelib19_ms2@spectra.data$Spectra.positive,
    gnps_hcecelllysatelipids_ms2@spectra.data$Spectra.positive,
    gnps_iobanhc_ms2@spectra.data$Spectra.positive,
    gnps_msmls_ms2@spectra.data$Spectra.positive,
    gnps_nihclinicalcollection1_ms2@spectra.data$Spectra.positive,
    gnps_nihclinicalcollection2_ms2@spectra.data$Spectra.positive,
    gnps_nutrimetabfem_ms2@spectra.data$Spectra.positive,
    gnps_pnnllipids_ms2@spectra.data$Spectra.positive,
    gnps_pusmsmls_ms2@spectra.data$Spectra.positive,
    gnps_sciexlibrary_ms2@spectra.data$Spectra.positive,
    gnps_selleckchemfdapart1_ms2@spectra.data$Spectra.positive
  )

gnps_ms2@spectra.data$Spectra.negative <-
  c(
    gnps_bilelib19_ms2@spectra.data$Spectra.negative,
    gnps_hcecelllysatelipids_ms2@spectra.data$Spectra.negative,
    gnps_iobanhc_ms2@spectra.data$Spectra.negative,
    gnps_msmls_ms2@spectra.data$Spectra.negative,
    gnps_nihclinicalcollection1_ms2@spectra.data$Spectra.negative,
    gnps_nihclinicalcollection2_ms2@spectra.data$Spectra.negative,
    gnps_nutrimetabfem_ms2@spectra.data$Spectra.negative,
    gnps_pnnllipids_ms2@spectra.data$Spectra.negative,
    gnps_pusmsmls_ms2@spectra.data$Spectra.negative,
    gnps_sciexlibrary_ms2@spectra.data$Spectra.negative,
    gnps_selleckchemfdapart1_ms2@spectra.data$Spectra.negative
  )

gnps_ms2@database.info

gnps_ms2@spectra.info$Compound.name <-
gnps_ms2@spectra.info$Compound.name %>%
  stringr::str_replace_all('\"', "")

gnps_ms2@spectra.info$Compound.name[grep(";", gnps_ms2@spectra.info$Compound.name)] <-
  gnps_ms2@spectra.info$Compound.name[grep(";", gnps_ms2@spectra.info$Compound.name)] %>%
  stringr::str_split(";") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

grep("[0-9]{1,6}-[0-9]{1,6}-[0-9]{1,6}",
     gnps_ms2@spectra.info$Compound.name, value = TRUE)

idx =grep("[0-9]{1,6}-[0-9]{1,6}-[0-9]{1,6}",
          gnps_ms2@spectra.info$Compound.name, value = FALSE)
idx
gnps_ms2@spectra.info$Compound.name <-
  gnps_ms2@spectra.info$Compound.name %>%
  stringr::str_replace_all("[0-9]{1,6}-[0-9]{1,6}-[0-9]{1,6}", "")

sum(gnps_ms2@spectra.info$Compound.name == "")
sum(is.na(gnps_ms2@spectra.info$Compound.name))

gnps_ms2 <-
gnps_ms2 %>%
  dplyr::filter(!is.na(Compound.name)) %>%
  dplyr::filter(Compound.name != "")

unique(gnps_ms2@spectra.info$Compound.name)

grep("name|Name",
     gnps_ms2@spectra.info$Compound.name, value = TRUE)

gnps_ms2@spectra.info$Compound.name[grep("Name|name", gnps_ms2@spectra.info$Compound.name)] <-
  gnps_ms2@spectra.info$Compound.name[grep("Name|name", gnps_ms2@spectra.info$Compound.name)] %>%
  stringr::str_split("Name\\:") %>%
  purrr::map(function(x){
    x[2]
  }) %>%
  unlist() %>%
  stringr::str_trim()

grep("low|poor",
     gnps_ms2@spectra.info$Compound.name, value = TRUE) %>% unique()

gnps_ms2@spectra.info$Compound.name[grep("low|poor", gnps_ms2@spectra.info$Compound.name)] <-
  gnps_ms2@spectra.info$Compound.name[grep("low|poor", gnps_ms2@spectra.info$Compound.name)] %>%
  stringr::str_split(" ", n = 2) %>%
  purrr::map(function(x){
    x[!stringr::str_detect(x, "low|poor")]
  }) %>%
  unlist() %>%
  stringr::str_trim()


grep("homodimer",
     gnps_ms2@spectra.info$Compound.name, value = TRUE) %>% unique()

gnps_ms2@spectra.info$Compound.name[grep("homodimer", gnps_ms2@spectra.info$Compound.name)] <-
  gnps_ms2@spectra.info$Compound.name[grep("homodimer", gnps_ms2@spectra.info$Compound.name)] %>%
  stringr::str_replace_all("\\[homodimer\\]", "") %>%
  stringr::str_trim()

idx <- which(!is.na(gnps_ms2@spectra.info$HMDB.ID))

gnps_ms2@spectra.info$Compound.name[idx] <-
hmdb_ms1@spectra.info$Compound.name[match(gnps_ms2@spectra.info$HMDB.ID[idx],
                                          hmdb_ms1@spectra.info$HMDB.ID)]

grep("\\[[0-9]{0,1}M",
     gnps_ms2@spectra.info$Compound.name, value = TRUE) %>% unique()

gnps_ms2@spectra.info$Compound.name[grep("\\[[0-9]{0,1}M", gnps_ms2@spectra.info$Compound.name)] <-
  gnps_ms2@spectra.info$Compound.name[grep("\\[[0-9]{0,1}M", gnps_ms2@spectra.info$Compound.name)] %>%
  stringr::str_split("\\[", n = 2) %>%
  purrr::map(function(x){
    x[1]
  }) %>%
  unlist() %>%
  stringr::str_trim()

grep("\\[",
     gnps_ms2@spectra.info$Compound.name, value = TRUE) %>% unique()

gnps_ms2@spectra.info$Compound.name[grep("\\[mostly precursor ion\\]", gnps_ms2@spectra.info$Compound.name)] <-
  gnps_ms2@spectra.info$Compound.name[grep("\\[mostly precursor ion\\]", gnps_ms2@spectra.info$Compound.name)] %>%
  stringr::str_replace_all("\\[mostly precursor ion\\]", "") %>%
  stringr::str_trim()

grep("\\[37Cl\\]",
     gnps_ms2@spectra.info$Compound.name, value = TRUE) %>% unique()

gnps_ms2@spectra.info$Compound.name[grep("\\[37Cl\\]", gnps_ms2@spectra.info$Compound.name)] <-
  gnps_ms2@spectra.info$Compound.name[grep("\\[37Cl\\]", gnps_ms2@spectra.info$Compound.name)] %>%
  stringr::str_replace_all("\\[37Cl\\]", "") %>%
  stringr::str_trim()

unique(gnps_ms2@spectra.info$Compound.name)

gnps_ms2@spectra.info$Compound.name <-
  gnps_ms2@spectra.info$Compound.name %>%
  stringr::str_trim()

dim(gnps_ms2)

gnps_ms2

table(gnps_ms2@spectra.info$Compound.name)[1:5]

gnps_ms2@spectra.info$Lab.ID[which(gnps_ms2@spectra.info$Compound.name == "(-)-Isopulegol")]
gnps_ms2@spectra.info$Polarity[which(gnps_ms2@spectra.info$Compound.name == "(-)-Isopulegol")]

ms2_spectrum1 <-
metid::get_ms2_spectrum(lab.id = "CCMSLIB00006674743",
                        database = gnps_ms2,
                        polarity = "positive", ce = "Unknown_1")

ms2_spectrum2 <-
  metid::get_ms2_spectrum(lab.id = "CCMSLIB00006674745",
                          database = gnps_ms2,
                          polarity = "positive", ce = "Unknown_1")

ms2_spectrum2 <-
  metid::get_ms2_spectrum(lab.id = "CCMSLIB00006674745",
                          database = gnps_ms2,
                          polarity = "positive", ce = "Unknown_1")

ms2_spectrum3 <-
  metid::get_ms2_spectrum(lab.id = "CCMSLIB00006674746",
                          database = gnps_ms2,
                          polarity = "positive", ce = "Unknown_1")

masstools::ms2_plot(spectrum1 = ms2_spectrum1,
                    spectrum2 = ms2_spectrum3)

get_spectra_match_score(
  exp.spectrum = ms2_spectrum1,
  lib.spectrum = ms2_spectrum2,
  remove_fragment_intensity_coutoff = 0.01
)

rownames(gnps_ms2@spectra.info) <- NULL

gnps_ms2

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

source(here::here("R/3_utils.R"))

intersect(colnames(gnps_ms2@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(gnps_ms2@spectra.info))

source(here::here("R/3_utils.R"))

gnps_ms2 <-
  update_metid_database_info(
    database = gnps_ms2,
    ref_database = kegg_ms1,
    by = c(
      "INCHIKEY.ID",
      "SMILES.ID",
      "INCHI.ID",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "IUPAC_name",
      "Traditional_IUPAC_name",
      "CHEMSPIDER.ID",
      "DRUGBANK.ID",
      "FOODB.ID",
      "PUBCHEM.ID",
      "CHEBI.ID",
      "BIOCYC.ID",
      "BIGG.ID",
      "WIKIPEDIA.ID",
      "METLIN.ID"
    ),
    combine_columns = c(
      "INCHIKEY.ID",
      "SMILES.ID",
      "INCHI.ID",
      "Synonyms",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "status",
      "IUPAC_name",
      "Traditional_IUPAC_name",
      "Kingdom",
      "Super_class",
      "Class",
      "Sub_class",
      "State",
      "Biospecimen_locations",
      "Cellular_locations",
      "Tissue_locations",
      "CHEMSPIDER.ID",
      "DRUGBANK.ID",
      "FOODB.ID",
      "PUBCHEM.ID",
      "CHEBI.ID",
      "BIOCYC.ID",
      "BIGG.ID",
      "WIKIPEDIA.ID",
      "METLIN.ID",
      "From_human"
    ),
    new_columns = c(
      "CHEMBL.ID",
      "LIPIDMAPS.ID",
      "LIPIDBANK.ID",
      "From_drug",
      "KEGG_DRUG.ID"
    )
  )

gnps_ms2

load(here::here("other_files/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

gnps_ms2 <-
  update_metid_database_source_system(
    database = gnps_ms2,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(gnps_ms2, file = "gnps_ms2.rda")

