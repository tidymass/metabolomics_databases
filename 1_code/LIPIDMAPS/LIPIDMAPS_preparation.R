no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

source("R/3_utils.R")

setwd('other_files/LIPIDMAPS/')

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
library(cinf)

library(ChemmineR)

# lipidmaps <-
#   read.SDFset(sdfstr = "structures.sdf", skipErrors = TRUE)
#
# lipidmaps_result <-
#   1:length(lipidmaps) %>%
#   purrr::map(function(i){
#     cat(i, " ")
#     x <- lipidmaps[[i]]
#   result <-
#     tryCatch(matrix(x[[4]], nrow = 1) %>%
#                as.data.frame(), error = NULL)
#   if(is.null(result)){
#     return(NULL)
#   }
#   colnames(result) <- names(x[[4]])
#   result
#   })
#
# all_column_name <-
# lipidmaps_result %>%
#   lapply(colnames) %>%
#   unlist() %>%
#   unique()
#
#
# lipidmaps_result <-
#   lipidmaps_result <-
#   1:length(lipidmaps_result) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     x <- lipidmaps_result[[i]]
#
#     diff_names <- setdiff(all_column_name, colnames(x))
#     if (length(diff_names) == 0) {
#       return(x[, all_column_name])
#     }
#     add_info <-
#       matrix(NA, nrow = 1, ncol = length(diff_names)) %>%
#       as.data.frame()
#     colnames(add_info) <- diff_names
#     cbind(x, add_info) %>%
#       dplyr::select(all_column_name)
#   })
#
#
# lipidmaps_result <-
# lipidmaps_result %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# lipidmaps_result <-
#   lipidmaps_result %>%
#   dplyr::rename(
#     Lab.ID = LM_ID,
#     PUBCHEM.ID = PUBCHEM_CID,
#     Compound.name = NAME,
#     Systematic.name = SYSTEMATIC_NAME,
#     Main_class_lipidmaps = MAIN_CLASS,
#     Sub_class_lipidmaps = SUB_CLASS,
#     Abbreviation = ABBREVIATION,
#     Class_level4 = CLASS_LEVEL4,
#     HMDB.ID = HMDB_ID,
#     Category = CATEGORY,
#     Formula = FORMULA,
#     mz = EXACT_MASS,
#     KEGG.ID = KEGG_ID,
#     CHEBI.ID = CHEBI_ID,
#     SWISSLIPIDS.ID = SWISSLIPIDS_ID,
#     LIPIDBANK.ID = LIPIDBANK_ID,
#     PLANTFA.ID = PLANTFA_ID,
#     SMILES.ID = SMILES,
#     INCHI.ID = INCHI,
#     INCHIKEY.ID = INCHI_KEY,
#     Synonyms = SYNONYMS
#   ) %>%
#   dplyr::mutate(
#     LIPIDMAPS.ID = Lab.ID,
#     CAS.ID = NA,
#     RT = NA,
#     mz.pos = NA,
#     mz.neg = NA,
#     Submitter = "LIPIDMAPS"
#   ) %>%
#   dplyr::select(Lab.ID,
#                 Compound.name,
#                 mz,
#                 RT,
#                 CAS.ID,
#                 HMDB.ID,
#                 KEGG.ID,
#                 Formula,
#                 mz.pos,
#                 mz.neg,
#                 Submitter,
#                 everything())
#
# lipidmaps_result$mz <-
#   as.numeric(lipidmaps_result$mz)
#
#
# ###remove mz is na
# sum(is.na(lipidmaps_result$mz))
#
# lipidmaps_result$Synonyms <-
#   lipidmaps_result$Synonyms %>%
#   stringr::str_replace_all('; ', "{}")
#
# lipidmaps_result[which(lipidmaps_result == "", arr.ind = TRUE)] <- NA
#
# load("../KEGG/kegg_ms1.rda")
# load("../HMDB/hmdb_ms1.rda")
#
# lipidmaps_result$KEGG.ID
#
# lipidmaps_result <-
#   lipidmaps_result %>%
#   dplyr::left_join(kegg_ms1@spectra.info[, c("KEGG.ID", "From_human", "From_drug")],
#                    by = c("KEGG.ID"),
#                    na_matches = "never")
#
# lipidmaps_result <-
#   lipidmaps_result %>%
#   dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "From_human")],
#                    by = c("HMDB.ID"),
#                    na_matches = "never")
#
# From_human <-
# lipidmaps_result[,c("From_human.x", "From_human.y")] %>%
#   apply(1, function(x){
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if(length(x) == 0){
#       return(NA)
#     }
#     return(x[1])
#   })
#
# lipidmaps_result <-
# lipidmaps_result %>%
#   dplyr::mutate(From_human = From_human) %>%
#   dplyr::select(-c(From_human.x, From_human.y))
#
# openxlsx::write.xlsx(lipidmaps_result, file = "lipidmaps_result.xlsx", asTable = TRUE)


library(metid)

lipidmaps_ms1 =
  construct_database(
    path = ".",
    version = "2022-04-19",
    metabolite.info.name = "lipidmaps_result.xlsx",
    source = "LIPIDMAPS",
    link = "https://www.lipidmaps.org/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

source(here::here("R/3_utils.R"))

intersect(colnames(lipidmaps_ms1@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(lipidmaps_ms1@spectra.info))

lipidmaps_ms1 <-
  update_metid_database_info(
    database = lipidmaps_ms1,
    ref_database = hmdb_ms1,
    by = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "INCHIKEY.ID",
      "INCHI.ID",
      "SMILES.ID",
      "PUBCHEM.ID",
      "CHEBI.ID"
    ),
    combine_columns = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "INCHIKEY.ID",
      "INCHI.ID",
      "SMILES.ID",
      "Synonyms",
      "PUBCHEM.ID",
      "CHEBI.ID",
      "From_human"
    ),
    new_columns = c(
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
      "BIOCYC.ID",
      "BIGG.ID",
      "WIKIPEDIA.ID",
      "METLIN.ID"
    )
  )

intersect(colnames(lipidmaps_ms1@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(lipidmaps_ms1@spectra.info))

source(here::here("R/3_utils.R"))

lipidmaps_ms1 <-
  update_metid_database_info(
    database = lipidmaps_ms1,
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
      "KEGG_DRUG.ID"
    )
  )

lipidmaps_ms1


load(here::here("other_files/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

lipidmaps_ms1 <-
  update_metid_database_source_system(
    database = lipidmaps_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(lipidmaps_ms1, file = "lipidmaps_ms1.rda")

