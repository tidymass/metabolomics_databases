####HMDB database: https://hmdb.ca/
setwd(masstools::get_project_wd())
rm(list = ls())

source("1_code/3_utils.R")

###load data

load("2_data/HMDB/MS1/hmdb_metabolite.rda")

dir.create("3_data_analysis/HMDB/MS1", showWarnings = FALSE)
setwd("3_data_analysis/HMDB/MS1")

library(tidyverse)
library(xml2)
library(stringr)

dim(hmdb_metabolite)

hmdb_metabolite <-
  hmdb_metabolite %>%
  dplyr::filter(!is.na(monisotopic_molecular_weight))

hmdb_metabolite[which(hmdb_metabolite == "", arr.ind = TRUE)] <- NA

hmdb_metabolite <-
  hmdb_metabolite %>%
  dplyr::rename(
    mz = average_molecular_weight,
    Create_date = creation_date,
    Updated_date = update_date,
    Lab.ID = accession,
    Compound.name = name,
    Description = description,
    Synonyms = synonyms,
    Formula = chemical_formula,
    IUPAC_name = iupac_name,
    Traditional_IUPAC_name = traditional_iupac,
    CAS.ID = cas_registry_number,
    SMILES.ID = smiles,
    INCHI.ID = inchi,
    INCHIKEY.ID = inchikey,
    Kingdom = kingdom,
    Super_class = super_class,
    Class = class,
    Sub_class = sub_class,
    State = state,
    Biospecimen_locations = biospecimen_locations,
    Cellular_locations = cellular_locations,
    Tissue_locations = tissue_locations,
    CHEMSPIDER.ID = chemspider_id,
    DRUGBANK.ID = drugbank_id,
    FOODB.ID = foodb_id,
    PUBCHEM.ID = pubchem_compound_id,
    CHEBI.ID = chebi_id,
    KEGG.ID = kegg_id,
    BIOCYC.ID = biocyc_id,
    BIGG.ID = bigg_id,
    WIKIPEDIA.ID = wikipedia_id,
    METLIN.ID = metlin_id
  ) %>%
  dplyr::mutate(
    HMDB.ID = Lab.ID,
    RT = NA,
    mz.pos = NA,
    mz.neg = NA,
    Submitter = "HMDB",
    From_human = "Yes"
  ) %>%
  dplyr::select(
    Lab.ID,
    Compound.name,
    mz,
    RT,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Formula,
    mz.pos,
    mz.neg,
    Submitter,
    everything()
  )

load("../../KEGG/kegg_ms1.rda")

idx1 <-
  match(hmdb_metabolite$CAS.ID,
        kegg_ms1@spectra.info$CAS.ID,
        incomparables = NA)

idx2 <-
  match(hmdb_metabolite$KEGG.ID,
        kegg_ms1@spectra.info$KEGG.ID,
        incomparables = NA)

idx <-
  data.frame(idx1, idx2) %>%
  keep_one_from_multiple()

KEGG.ID <-
  data.frame(hmdb_metabolite$KEGG.ID, kegg_ms1@spectra.info$KEGG.ID[idx]) %>%
  keep_one_from_multiple()
hmdb_metabolite$KEGG.ID <- KEGG.ID

CAS.ID <-
  data.frame(hmdb_metabolite$CAS.ID, kegg_ms1@spectra.info$CAS.ID[idx]) %>%
  keep_one_from_multiple()
hmdb_metabolite$CAS.ID <- CAS.ID

openxlsx::write.xlsx(hmdb_metabolite,
                     file = "hmdb_metabolite.xlsx",
                     asTable = TRUE, overwrite = TRUE)

library(metid)

hmdb_ms1 <-
  construct_database(
    path = ".",
    version = "2022-04-11",
    metabolite.info.name = "hmdb_metabolite.xlsx",
    source = "HMDB",
    link = "https://hmdb.ca/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

hmdb_ms1@spectra.info$mz <-
hmdb_ms1@spectra.info$monisotopic_molecular_weight

load(here::here("2_data/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

hmdb_ms1 <-
  update_metid_database_source_system(
    database = hmdb_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(hmdb_ms1, file = "hmdb_ms1.rda")
