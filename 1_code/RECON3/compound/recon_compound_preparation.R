no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("2_data/RECON3/compound/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
library(R.matlab)

rm(list = ls())

result <-
  readr::read_tsv("recon-store-metabolites-1.tsv")

recon3_compound_data <-
  result %>%
  dplyr::rename(
    Compound.name = fullName,
    KEGG_ID = keggId,
    PUBCHEM_ID = pubChemId,
    CHEBl_ID = cheBlId,
    HMDB_ID = hmdb,
    FOODB_ID = food_db,
    CHEMSPIDER_ID = chemspider,
    BIGG_ID = biggId,
    wikipedia_ID = wikipedia,
    DRUGBANK_ID = drugbank,
    METLIN_ID = metlin,
    CAS_ID = casRegistry,
    INCHIKEY = inchiKey,
    SMILE = smile
  )

save(recon3_compound_data, file = "recon3_compound_data")
