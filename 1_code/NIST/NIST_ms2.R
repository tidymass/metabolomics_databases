###no source
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())

load("2_data/KEGG/kegg_ms1.rda")
load("2_data/HMDB/MS1/hmdb_ms1.rda")
load("2_data/NIST/nistDatabase0.0.3")

source("1_code/3_utils.R")

setwd("2_data/NIST/")

library(tidyverse)
library(tidymass)

# data <-
#   nistDatabase0.0.3@spectra.info
#
# data[which(data == "", arr.ind = TRUE)] <- NA
# data[which(data == "NA", arr.ind = TRUE)] <- NA
#
# sum(is.na(data$KEGG.ID))
# sum(is.na(data$HMDB.ID))
#
# # data_2 <-
# #   metid::read_msp(file = "MSMS-Pos-VS14.msp", threads = 5)
# #
# # data_2 %>%
# #   purrr::map(function(x){
# #     names(x$info)
# #   }) %>%
# #   unlist() %>%
# #   unique()
# #
# # spectra.info <-
# #   1:length(data_2) %>%
# #   purrr::map(function(i){
# #     cat(i, " ")
# #     x <- data_2[[i]]
# #     result <-
# #     matrix(x$info, nrow = 1) %>%
# #   as.data.frame()
# #     colnames(result) <- names(x$info)
# #     result
# #   }) %>%
# #   dplyr::bind_rows() %>%
# #   as.data.frame()
#
# sum(is.na(data$HMDB.ID))
# sum(is.na(data$KEGG.ID))
#
# write.csv(data, "data.csv", row.names = FALSE)
#
data <-
  readr::read_csv("data_manual.csv")

data$CAS.ID <-
data$CAS.ID %>%
  purrr::map(function(x) {
    if (is.na(x)) {
      return(x)
    }
    if (stringr::str_detect(x, "\\)")) {
      x <-
        stringr::str_split(x, "\\(")[[1]] %>%
        stringr::str_replace_all("\\)", "") %>%
        stringr::str_trim() %>%
        paste(collapse = "{}")
      return(x)
    } else{
      return(x)
    }
  }) %>%
  unlist()

data$HMDB.ID <-
  data$HMDB.ID %>%
  purrr::map(function(x) {
    if (is.na(x)) {
      return(x)
    }
    if (stringr::str_detect(x, "\\/")) {
      x <-
        stringr::str_split(x, "\\/")[[1]] %>%
        stringr::str_trim()

      x <-
      lapply(x, function(y){
        if(nchar(y) == 9){
          stringr::str_replace(y, "HMDB", "HMDB00")
          }
      }) %>%
        unlist()

      x <-
        paste(x, collapse = "{}")
      return(x)
    } else{
      if(nchar(x) == 9){
        x <-
          stringr::str_replace(x, "HMDB", "HMDB00")
      }
      return(x)
    }
  }) %>%
  unlist()

idx <-
which(is.na(data$KEGG.ID))

kegg_id <-
hmdb_ms1@spectra.info$KEGG.ID[match(unlist(lapply(stringr::str_split(data$HMDB.ID[idx], "\\{\\}"), function(x)x[1])),
      hmdb_ms1@spectra.info$HMDB.ID)]

data$KEGG.ID[idx] <- kegg_id

head(data[idx,])

idx <-
match(unlist(lapply(stringr::str_split(data$HMDB.ID, "\\{\\}"), function(x)x[1])),
      hmdb_ms1@spectra.info$HMDB.ID)

new_name <-
setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(data))

new_data <-
hmdb_ms1@spectra.info[idx,] %>%
  dplyr::select(new_name) %>%
  dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))

new_data$Submitter <- "NIST"

data <-
  cbind(data, new_data)

nist_ms2 <- nistDatabase0.0.3

nist_ms2@spectra.info <- data

nist_ms2@database.info$Version <- "20220425"
nist_ms2@database.info$Email <- "shenxt@stanford.edu"


load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

intersect(colnames(nist_ms2@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(nist_ms2@spectra.info))

nist_ms2 <-
  update_metid_database_info(
    database = nist_ms2,
    ref_database = hmdb_ms1,
    by = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "METLIN.ID"
    ),
    combine_columns = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "METLIN.ID"
    ),
    new_columns = c())

intersect(colnames(nist_ms2@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(nist_ms2@spectra.info))

source(here::here("1_code/3_utils.R"))

nist_ms2 <-
  update_metid_database_info(
    database = nist_ms2,
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

load(here::here("2_data/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

nist_ms2 <-
  update_metid_database_source_system(
    database = nist_ms2,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(nist_ms2, file = "nist_ms2.rda")



