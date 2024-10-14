###
no_source
rm(list = ls())
setwd(masstools::get_project_wd())
source("R/3_utils.R")
setwd('other_files/T3DB')
library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)

library(rjson)

# t3db_2 <-
#   readr::read_csv("toxins.csv")
#
# t3db <-
#   jsonlite::fromJSON(txt = "toxins.json")
#
# t3db <-
#   t3db %>%
#   dplyr::left_join(t3db_2[, c("T3DB ID", "Class")], by = c("title" = "T3DB ID")) %>%
#   dplyr::filter(Class == "SmallMolecule")
#
# colnames(t3db)
#
# t3db <-
#   t3db %>%
#   dplyr::filter(Class == "SmallMolecule") %>%
#   dplyr::select(-id) %>%
#   dplyr::rename(
#     Lab.ID = title,
#     PUBCHEM.ID = pubchem_id,
#     Compound.name = common_name,
#     HMDB.ID = hmdb_id,
#     CAS.ID = cas,
#     Formula = chemical_formula,
#     mz = moldb_mono_mass,
#     WIKIPEDIA.ID = wikipedia,
#     KEGG.ID = kegg_compound_id,
#     UNIPROT.ID = uniprot_id,
#     OMIN.ID = omim_id,
#     CHEBI.ID = chebi_id,
#     BIOCYC.ID = biocyc_id,
#     CTD.ID = ctd_id,
#     STITCH.ID = stitch_id,
#     DRUGBANK.ID = drugbank_id,
#     PDB.ID = pdb_id,
#     CHEMBL.ID = chembl_id,
#     CHEMSPIDER.ID = chemspider_id,
#     BIODB.ID = biodb_id,
#     SMILES.ID = moldb_smiles,
#     INCHI.ID = moldb_inchi,
#     INCHIKEY.ID = moldb_inchikey,
#     Create_date = created_at,
#     Updated_date = updated_at,
#     Average.mass = moldb_average_mass,
#     Synonyms = synonyms_list
#   ) %>%
#   dplyr::mutate(
#     T3DB.ID = Lab.ID,
#     RT = NA,
#     mz.pos = NA,
#     mz.neg = NA,
#     Submitter = "T3DB"
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
# # t3db$Formula
# t3db$From_environment <- "Yes"
#
# t3db$toxicity <-
#   t3db$toxicity %>%
#   stringr::str_replace_all('\r\n', "{}")
#
# t3db$Synonyms <-
# t3db$Synonyms %>%
#   stringr::str_replace_all('\r\n', "{}")
#
# t3db$types <-
#   t3db$types %>%
#   lapply(function(x){
#     if(nrow(x) == 0){
#       return(NA)
#     }
#     paste(x$type_name, collapse = "{}")
#   }) %>%
#   unlist()
#
# t3db$cellular_locations <-
#   t3db$cellular_locations %>%
#   lapply(function(x){
#     if(nrow(x) == 0){
#       return(NA)
#     }
#     paste(x$name, collapse = "{}")
#   }) %>%
#   unlist()
#
# t3db <-
# t3db %>%
#   dplyr::rename(Cellular_locations = cellular_locations)
#
# t3db$tissues <-
#   t3db$tissues %>%
#   lapply(function(x){
#     if(nrow(x) == 0){
#       return(NA)
#     }
#     paste(x$name, collapse = "{}")
#   }) %>%
#   unlist()
#
# t3db <-
#   t3db %>%
#   dplyr::rename(Tissues = tissues)
#
# t3db$Pathway_name <-
#   t3db$pathways %>%
#   lapply(function(x){
#     if(nrow(x) == 0){
#       return(NA)
#     }
#     paste(x$name, collapse = "{}")
#   }) %>%
#   unlist()
#
# t3db$Pathway_KEGG.ID <-
#   t3db$pathways %>%
#   lapply(function(x){
#     if(nrow(x) == 0){
#       return(NA)
#     }
#     paste(x$kegg_id[!is.na(x$kegg_id)], collapse = "{}")
#   }) %>%
#   unlist()
#
# t3db$Pathway_SMPDB.ID <-
#   t3db$pathways %>%
#   lapply(function(x){
#     if(nrow(x) == 0){
#       return(NA)
#     }
#     x$kegg_id[x$kegg_id == ""] <- NA
#     x$smpdb_id[x$smpdb_id == ""] <- NA
#     if(all(is.na(x$smpdb_id))){
#       return(NA)
#     }
#     paste(x$smpdb_id[!is.na(x$smpdb_id)], collapse = "{}")
#   }) %>%
#   unlist()
#
# t3db <-
#   t3db %>%
#   dplyr::select(-pathways)
#
# t3db[which(t3db == "", arr.ind = TRUE)] <- NA
#
# load("../KEGG/kegg_ms1.rda")
#
# t3db$KEGG.ID
#
# t3db <-
#   t3db %>% dplyr::left_join(kegg_ms1@spectra.info[, c("KEGG.ID", "From_human", "From_drug")],
#                             by = c("KEGG.ID"),
#                             na_matches = "never")
#
# t3db$From_drug[is.na(t3db$From_drug)] <- "No"
# t3db$From_human[is.na(t3db$From_human)] <- "No"
#
# openxlsx::write.xlsx(t3db, file = "t3db.xlsx", asTable = TRUE)

library(metid)

t3db_ms1 =
  construct_database(
    path = ".",
    version = "2022-04-08",
    metabolite.info.name = "t3db.xlsx",
    source = "T3DB",
    link = "http://www.t3db.ca/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

intersect(colnames(t3db_ms1@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(t3db_ms1@spectra.info))

t3db_ms1 <-
  update_metid_database_info(
    database = t3db_ms1,
    ref_database = hmdb_ms1,
    by = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "PUBCHEM.ID",
      "WIKIPEDIA.ID",
      "CHEBI.ID",
      "BIOCYC.ID",
      "DRUGBANK.ID",
      "SMILES.ID",
      "INCHI.ID",
      "INCHIKEY.ID",
      "CHEMSPIDER.ID"
    ),
    combine_columns = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "PUBCHEM.ID",
      "WIKIPEDIA.ID",
      "CHEBI.ID",
      "BIOCYC.ID",
      "DRUGBANK.ID",
      "SMILES.ID",
      "INCHI.ID",
      "INCHIKEY.ID",
      "CHEMSPIDER.ID",
      "Synonyms",
      "Cellular_locations"
    ),
    new_columns = c(
      "status",
      "Description",
      "monisotopic_molecular_weight",
      "IUPAC_name",
      "Traditional_IUPAC_name",
      "Kingdom",
      "Super_class",
      "Sub_class",
      "State",
      "Biospecimen_locations",
      "Tissue_locations",
      "FOODB.ID",
      "BIGG.ID",
      "METLIN.ID"
    ))

intersect(colnames(t3db_ms1@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(t3db_ms1@spectra.info))

source(here::here("R/3_utils.R"))

t3db_ms1 <-
  update_metid_database_info(
    database = t3db_ms1,
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
      "LIPIDMAPS.ID",
      "LIPIDBANK.ID",
      "KEGG_DRUG.ID"
    )
  )

load(here::here("other_files/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

t3db_ms1 <-
  update_metid_database_source_system(
    database = t3db_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(t3db_ms1, file = "t3db_ms1.rda")
