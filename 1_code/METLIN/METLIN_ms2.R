###no source
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())
source("R/3_utils.R")

load("other_files/KEGG/kegg_ms1.rda")
load("other_files/HMDB/MS1/hmdb_ms1.rda")
load("other_files/METLIN/metlinDatabase0.0.3")

setwd("other_files/METLIN/")

library(tidyverse)
library(tidymass)
#
# data <-
#   metlinDatabase0.0.3@spectra.info
#
# data[which(data == "NA", arr.ind = TRUE)] <- NA
#
# sum(is.na(data$KEGG.ID))
#
# hmdb_id1 <-
#   data$HMDB.ID
#
# hmdb_id2 <-
#   kegg_ms1@spectra.info$HMDB.ID[match(data$KEGG.ID, kegg_ms1@spectra.info$KEGG.ID)]
#
# hmdb_id3 <-
#   hmdb_ms1@spectra.info$HMDB.ID[match(data$KEGG.ID, hmdb_ms1@spectra.info$KEGG.ID)]
#
# hmdb_id4 <-
#   1:nrow(data) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x = data$KEGG.ID[i]
#     y = c(hmdb_id1[i],
#           hmdb_id2[i],
#           hmdb_id3[i])
#     y = y[!is.na(y)] %>%
#       unique()
#     if (is.na(x)) {
#       return(NA)
#     } else{
#       if(length(y) != 0){
#         return(NA)
#       }
#       masstools::trans_ID(
#         query = x,
#         from = "KEGG",
#         to = "Human Metabolome Database",
#         top = 1
#       )$`Human Metabolome Database`
#     }
#   })
#
# hmdb_id4 <- unlist(hmdb_id4)
#
# length(hmdb_id1)
# length(hmdb_id2)
# length(hmdb_id3)
# length(hmdb_id4)
#
# hmdb_id <-
#   data.frame(hmdb_id1,
#              hmdb_id2,
#              hmdb_id3,
#              hmdb_id4)
#
# head(hmdb_id)
#
# hmdb_id <-
# hmdb_id %>%
#   apply(1, function(x){
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     x <- x[x!="NA"]
#     if(length(x) == 0){
#       return(NA)
#     }else{
#       return(x[1])
#     }
#   })
#
# data$HMDB.ID <- hmdb_id
#
#
# write.csv(data, file = "data.csv", row.names = FALSE)

data <-
  readr::read_csv("data_manual.csv")

which(is.na(data$HMDB.ID) & is.na(data$KEGG.ID))

idx <-
  which(!is.na(data$HMDB.ID) & is.na(data$KEGG.ID))

kegg_id <-
hmdb_ms1@spectra.info$KEGG.ID[match(data$HMDB.ID[idx], hmdb_ms1@spectra.info$HMDB.ID)]

data$KEGG.ID[idx] <- kegg_id

metlin_ms2 <- metlinDatabase0.0.3

metlin_ms2@spectra.info <- data

metlin_ms2@database.info$Version <- "20220425"
metlin_ms2@database.info$Email <- "shenxt@stanford.edu"

metlin_ms2

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

intersect(colnames(metlin_ms2@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(metlin_ms2@spectra.info))

metlin_ms2 <-
  update_metid_database_info(
    database = metlin_ms2,
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
    new_columns = c(
      "status",
      "secondary_accessions",
      "Description",
      "Synonyms",
      "monisotopic_molecular_weight",
      "IUPAC_name",
      "Traditional_IUPAC_name",
      "SMILES.ID",
      "INCHI.ID",
      "INCHIKEY.ID",
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
      "From_human"
    ))

# metlin_ms2@spectra.info$Synonyms

intersect(colnames(metlin_ms2@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(metlin_ms2@spectra.info))

source(here::here("R/3_utils.R"))

metlin_ms2 <-
  update_metid_database_info(
    database = metlin_ms2,
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

load(here::here("other_files/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

metlin_ms2 <-
  update_metid_database_source_system(
    database = metlin_ms2,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(metlin_ms2, file = "metlin_ms2.rda")
