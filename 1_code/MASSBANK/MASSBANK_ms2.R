###source
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
load("2_data/HMDB/MS1/hmdb_ms1.rda")
load("2_data/KEGG/kegg_ms1.rda")
source("R/read_msp_data.R")
source("1_code/3_utils.R")
source("R/10_MASSBANK.R")
setwd("2_data/MASSBANK")

library(tidyverse)
library(tidymass)
library(plyr)

# data1 <-
#   read_msp_data_massbank(file = "MassBank_NIST.msp", threads = 5)
#
# massbank_nist_ms2 <-
#   convert_massbank2metid(data = data1,
#                          path = "NIST",
#                          source = "nist",
#                          threads = 5)
#
# load("NIST/massbank_ms2")
# massbank_nist_ms2 <- massbank_ms2
#
# colnames(massbank_nist_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# idx1 <-
#   match(
#     massbank_nist_ms2@spectra.info$INCHI.ID,
#     hmdb_ms1@spectra.info$INCHI.ID,
#     incomparables = NA
#   )
#
# idx2 <-
#   match(
#     massbank_nist_ms2@spectra.info$INCHIKEY.ID,
#     hmdb_ms1@spectra.info$INCHIKEY.ID,
#     incomparables = NA
#   )
#
# idx3 <-
#   match(
#     massbank_nist_ms2@spectra.info$SMILES.ID,
#     hmdb_ms1@spectra.info$SMILES.ID,
#     incomparables = NA
#   )
#
# idx <-
#   data.frame(idx1, idx2, idx3) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
#
# cbind(
#   x = massbank_nist_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(massbank_nist_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   massbank_nist_ms2@spectra.info
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
# # cbind(spectra.info$SMILES.ID,
# #       hmdb_ms1@spectra.info$SMILES.ID[idx])
#
#
# massbank_nist_ms2@spectra.info <-
#   spectra.info
#
# save(massbank_nist_ms2, file = "massbank_nist_ms2.rda")
#
# data2 <-
#   read_msp_data_massbank(file = "MassBank_RIKEN.msp", threads = 5)
#
# massbank_riken_ms2 <-
#   convert_massbank2metid(data = data2,
#                          threads = 5,
#                          source = "riken")
#
# new_spectra_info <-
# seq_along(massbank_riken_ms2@spectra.info$Comment) %>%
#   purrr::map(function(i){
#     # cat(i, " ")
#     x <- massbank_riken_ms2@spectra.info$Comment[i]
#     x <-
#     stringr::str_split(x, ";")[[1]] %>%
#       stringr::str_trim()
#     DB <- x[stringr::str_detect(x, "^DB")] %>%
#       stringr::str_replace("^DB#=", "") %>%
#       stringr::str_trim()
#     DB <- ifelse(length(DB) == 0, NA, DB)
#     origin <- x[stringr::str_detect(x, "^origin=")] %>%
#       stringr::str_replace("^origin=", "") %>%
#       stringr::str_trim()
#     origin <- ifelse(length(origin) == 0, NA, origin)
#     Comment_confidence <- x[stringr::str_detect(x, "^Annotation")] %>%
#       stringr::str_replace("^Annotation", "") %>%
#       stringr::str_trim()
#     Comment_confidence <- ifelse(length(Comment_confidence) == 0, NA, Comment_confidence)
#     data.frame(DB, origin, Comment_confidence)
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# massbank_riken_ms2@spectra.info$Lab.ID <- new_spectra_info$DB
# massbank_riken_ms2@spectra.info$MASSBANK.ID <- new_spectra_info$DB
# massbank_riken_ms2@spectra.info$Submitter_team <- new_spectra_info$origin
# massbank_riken_ms2@spectra.info$Comment_confidence <- new_spectra_info$Comment_confidence
#
# colnames(massbank_riken_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# idx1 <-
#   match(massbank_riken_ms2@spectra.info$INCHI.ID,
#         hmdb_ms1@spectra.info$INCHI.ID,
#         incomparables = NA)
#
# idx2 <-
#   match(massbank_riken_ms2@spectra.info$INCHIKEY.ID,
#         hmdb_ms1@spectra.info$INCHIKEY.ID,
#         incomparables = NA)
#
# idx3 <-
#   match(massbank_riken_ms2@spectra.info$SMILES.ID,
#         hmdb_ms1@spectra.info$SMILES.ID,
#         incomparables = NA)
#
# idx <-
#   data.frame(idx1, idx2, idx3) %>%
#   apply(1, function(x) {
#     x <- as.numeric(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
#
# cbind(
#   x = massbank_riken_ms2@spectra.info$Compound.name,
#   y = hmdb_ms1@spectra.info$Compound.name[idx]
# ) %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(y))
#
# colnames(massbank_riken_ms2@spectra.info)
# colnames(hmdb_ms1@spectra.info)
#
# spectra.info <-
#   massbank_riken_ms2@spectra.info
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
# massbank_riken_ms2@spectra.info <-
#   spectra.info
#
# save(massbank_riken_ms2, file = "massbank_riken_ms2.rda")



# ####combine riken and nist
# load("massbank_riken_ms2.rda")
# load("massbank_nist_ms2.rda")
#
# intersect_compound_name <-
#   intersect(
#     massbank_nist_ms2@spectra.info$Compound.name,
#     massbank_riken_ms2@spectra.info$Compound.name
#   )
#
# intersect_compound_name[1]
# Lab.ID1 <-
#   massbank_nist_ms2@spectra.info %>%
#   dplyr::filter(Compound.name == intersect_compound_name[1]) %>%
#   pull(Lab.ID)
#
# Lab.ID2 <-
#   massbank_riken_ms2@spectra.info %>%
#   dplyr::filter(Compound.name == intersect_compound_name[1]) %>%
#   pull(Lab.ID)
#
# Lab.ID1
# Lab.ID2
#
# which(names(massbank_nist_ms2@spectra.data$Spectra.positive) == Lab.ID1[1])
# which(names(massbank_riken_ms2@spectra.data$Spectra.positive) == Lab.ID2[1])
#
# masstools::ms2_plot(
#   spectrum1 = massbank_nist_ms2@spectra.data$Spectra.positive[[1]][[1]],
#   spectrum2 = massbank_riken_ms2@spectra.data$Spectra.positive[[1]][[1]]
# )
#
# which(names(massbank_nist_ms2@spectra.data$Spectra.negative) == Lab.ID1[4])
# which(names(massbank_riken_ms2@spectra.data$Spectra.negative) == Lab.ID2[4])
#
#
# ##remove redundant spectra
# intersect_compound_name
# massbank_riken_ms2 <-
#   massbank_riken_ms2 %>%
#   dplyr::filter(!Compound.name %in% intersect_compound_name)
#
# massbank_ms2 <-
#   massbank_nist_ms2
#
# grep("[0-9]{1,5}\\.[0-9]{1,5}_[0-9]{1,5}\\.[0-9]{1,5}",
#      massbank_ms2@spectra.info$Compound.name,
#      value = FALSE)
#
# grep("[0-9]{1,5}\\.[0-9]{1,5}_[0-9]{1,5}\\.[0-9]{1,5}",
#      massbank_ms2@spectra.info$Compound.name,
#      value = TRUE)
#
# massbank_ms2@spectra.info$Compound.name[!is.na(massbank_ms2@spectra.info$HMDB.ID)] <-
#   hmdb_ms1@spectra.info$Compound.name[match(massbank_ms2@spectra.info$HMDB[!is.na(massbank_ms2@spectra.info$HMDB.ID)],
#                                             hmdb_ms1@spectra.info$HMDB.ID)]
#
#
# rownames(massbank_ms2@spectra.info) <- NULL
#
# save(massbank_ms2, file = "massbank_ms2.rda")
load("massbank_ms2.rda")
massbank_ms2

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

intersect(colnames(massbank_ms2@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(massbank_ms2@spectra.info))

massbank_ms2 <-
  update_metid_database_info(
    database = massbank_ms2,
    ref_database = hmdb_ms1,
    by = c(
      "Compound.name",
      "INCHI.ID",
      "INCHIKEY.ID",
      "SMILES.ID",
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
      "Compound.name",
      "INCHI.ID",
      "INCHIKEY.ID",
      "SMILES.ID",
      "Synonyms",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "status",
      "Description",
      "monisotopic_molecular_weight",
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
    new_columns = c())

intersect(colnames(massbank_ms2@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(massbank_ms2@spectra.info))

source(here::here("1_code/3_utils.R"))

massbank_ms2 <-
  update_metid_database_info(
    database = massbank_ms2,
    ref_database = kegg_ms1,
    by = c(
      "Compound.name",
      "IUPAC_name",
      "Traditional_IUPAC_name",
      "CAS.ID",
      "SMILES.ID",
      "INCHI.ID",
      "INCHIKEY.ID",
      "FOODB.ID",
      "HMDB.ID",
      "PUBCHEM.ID",
      "CHEMSPIDER.ID",
      "KEGG.ID",
      "CHEBI.ID",
      "BIOCYC.ID",
      "WIKIPEDIA.ID",
      "DRUGBANK.ID",
      "BIGG.ID",
      "METLIN.ID"
    ),
    combine_columns = c(
      "Compound.name",
      "Synonyms",
      "IUPAC_name",
      "Traditional_IUPAC_name",
      "CAS.ID",
      "SMILES.ID",
      "INCHI.ID",
      "INCHIKEY.ID",
      "Kingdom",
      "Super_class",
      "Class",
      "Sub_class",
      "State",
      "FOODB.ID",
      "HMDB.ID",
      "PUBCHEM.ID",
      "CHEMSPIDER.ID",
      "KEGG.ID",
      "CHEBI.ID",
      "BIOCYC.ID",
      "WIKIPEDIA.ID",
      "status",
      "Biospecimen_locations",
      "Cellular_locations",
      "Tissue_locations",
      "DRUGBANK.ID",
      "BIGG.ID",
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

massbank_ms2

load(here::here("2_data/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

massbank_ms2 <-
  update_metid_database_source_system(
    database = massbank_ms2,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(massbank_ms2, file = "massbank_ms2.rda")
