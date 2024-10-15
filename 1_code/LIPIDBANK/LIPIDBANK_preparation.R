no_source()

setwd(masstools::get_project_wd())

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
library(cinf)
rm(list = ls())

source("R/9_LIPIDBANK.R")
source("1_code/3_utils.R")
setwd('2_data/LIPIDBANK/')
library(ChemmineR)

# lipidbank <-
#   request_lipidbank(lipid_class = "All data")
#
# openxlsx::write.xlsx(lipidbank, file = "lipidbank.xlsx", asTable = TRUE)

library(metid)

lipidbank_ms1 =
  construct_database(
    path = ".",
    version = "2022-04-19",
    metabolite.info.name = "lipidbank.xlsx",
    source = "lipidbank",
    link = "https://www.lipidbank.org/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

source(here::here("1_code/3_utils.R"))

intersect(colnames(lipidbank_ms1@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(lipidbank_ms1@spectra.info))

lipidbank_ms1 <-
  update_metid_database_info(
    database = lipidbank_ms1,
    ref_database = hmdb_ms1,
    by = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID"
    ),
    combine_columns = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "Formula",
      "Synonyms",
      "From_human"
    ),
    new_columns = c(
      "status",
      "Description",
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
      "METLIN.ID"
    )
  )

intersect(colnames(lipidbank_ms1@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(lipidbank_ms1@spectra.info))

source(here::here("1_code/3_utils.R"))

lipidbank_ms1 <-
  update_metid_database_info(
    database = lipidbank_ms1,
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
      "KEGG_DRUG.ID"
    )
  )

lipidbank_ms1

load(here::here("2_data/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

lipidbank_ms1 <-
  update_metid_database_source_system(
    database = lipidbank_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(lipidbank_ms1, file = "lipidbank_ms1.rda")
