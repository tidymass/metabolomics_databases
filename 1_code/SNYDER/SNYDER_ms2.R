###no source
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())
load("other_files/SNYDER/MS1/mpsnyder_hilic_ms1.rda")
load("other_files/SNYDER/MS1/mpsnyder_rplc_ms1.rda")

load("other_files/SNYDER/MS1/snyder_database_hilic0.0.3.rda")
load("other_files/SNYDER/MS1/snyder_database_rplc0.0.3.rda")

load("other_files/HMDB/MS1/hmdb_ms1.rda")

source("R/3_utils.R")

setwd("other_files/SNYDER/MS2/")

library(tidyverse)
library(tidymass)

dim(snyder_database_rplc0.0.3@spectra.info)

dim(mpsnyder_rplc_ms1@spectra.info)

mpsnyder_rplc_ms2 <- snyder_database_rplc0.0.3

mpsnyder_rplc_ms2@database.info$Version <- "20220424"
mpsnyder_rplc_ms2@database.info$Source <-
  mpsnyder_rplc_ms1@database.info$Source
mpsnyder_rplc_ms2@database.info$Link <-
  mpsnyder_rplc_ms1@database.info$Link
mpsnyder_rplc_ms2@database.info$Email <-
  mpsnyder_rplc_ms1@database.info$Email
mpsnyder_rplc_ms2@database.info$RT <-
  mpsnyder_rplc_ms1@database.info$RT

mpsnyder_rplc_ms2@spectra.info <- mpsnyder_rplc_ms1@spectra.info

mpsnyder_hilic_ms2 <- snyder_database_hilic0.0.3

mpsnyder_hilic_ms2@database.info$Version <- "20220424"
mpsnyder_hilic_ms2@database.info$Source <-
  mpsnyder_hilic_ms1@database.info$Source
mpsnyder_hilic_ms2@database.info$Link <-
  mpsnyder_hilic_ms1@database.info$Link
mpsnyder_hilic_ms2@database.info$Email <-
  mpsnyder_hilic_ms1@database.info$Email
mpsnyder_hilic_ms2@database.info$RT <-
  mpsnyder_hilic_ms1@database.info$RT

mpsnyder_hilic_ms2@spectra.info <- mpsnyder_hilic_ms1@spectra.info

load("../../HMDB/MS1/hmdb_ms1.rda")
load("../../KEGG/kegg_ms1.rda")

intersect(colnames(mpsnyder_hilic_ms2@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(mpsnyder_hilic_ms2@spectra.info))

mpsnyder_hilic_ms2 <-
  update_metid_database_info(
    database = mpsnyder_hilic_ms2,
    ref_database = hmdb_ms1,
    by = c("Compound.name",
           "CAS.ID",
           "HMDB.ID",
           "KEGG.ID",
           "METLIN.ID"),
    combine_columns = c("Compound.name",
                        "CAS.ID",
                        "HMDB.ID",
                        "KEGG.ID",
                        "METLIN.ID"),
    new_columns = c()
  )

intersect(colnames(mpsnyder_hilic_ms2@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(mpsnyder_hilic_ms2@spectra.info))

source(here::here("R/3_utils.R"))

mpsnyder_hilic_ms2 <-
  update_metid_database_info(
    database = mpsnyder_hilic_ms2,
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



####change formula and mz
spectra.info <-
  mpsnyder_hilic_ms2@spectra.info

idx <-
  match(spectra.info$HMDB.ID, hmdb_ms1@spectra.info$HMDB.ID)

temp <-
data.frame(spectra.info[,c("Lab.ID","HMDB.ID","mz", "Formula", "Compound.name")],
           hmdb_ms1@spectra.info[idx,c("mz", "Formula", "Compound.name")]) %>%
  dplyr::mutate(mz_error = abs(mz - mz.1)*10^6/mz) %>%
  dplyr::filter(is.na(mz_error) | mz_error > 10)

write.csv(temp, file = "temp.csv", row.names = FALSE)

temp <- readr::read_csv("temp_manual.csv")

temp$mz
temp$Formula

spectra.info$mz[match(temp$Lab.ID, spectra.info$Lab.ID)] <-
temp$mz
spectra.info$Formula[match(temp$Lab.ID, spectra.info$Lab.ID)] <-
  temp$Formula

mpsnyder_hilic_ms2@spectra.info <-
  spectra.info

write.csv(spectra.info, file = "spectra.info2.csv", row.names = FALSE)

save(mpsnyder_hilic_ms2, file = "mpsnyder_hilic_ms2.rda")



load("mpsnyder_rplc_ms2.rda")

load("../../HMDB/MS1/hmdb_ms1.rda")
load("../../KEGG/kegg_ms1.rda")

intersect(colnames(mpsnyder_rplc_ms2@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(mpsnyder_rplc_ms2@spectra.info))

mpsnyder_rplc_ms2 <-
  update_metid_database_info(
    database = mpsnyder_rplc_ms2,
    ref_database = hmdb_ms1,
    by = c("Compound.name",
           "CAS.ID",
           "HMDB.ID",
           "KEGG.ID",
           "METLIN.ID"),
    combine_columns = c("Compound.name",
                        "CAS.ID",
                        "HMDB.ID",
                        "KEGG.ID",
                        "METLIN.ID"),
    new_columns = c()
  )

intersect(colnames(mpsnyder_rplc_ms2@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(mpsnyder_rplc_ms2@spectra.info))

source(here::here("R/3_utils.R"))

mpsnyder_rplc_ms2 <-
  update_metid_database_info(
    database = mpsnyder_rplc_ms2,
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



####change formula and mz
spectra.info <-
  mpsnyder_rplc_ms2@spectra.info


idx <-
  match(spectra.info$HMDB.ID, hmdb_ms1@spectra.info$HMDB.ID)

temp2 <-
  data.frame(spectra.info[,c("Lab.ID","HMDB.ID","mz", "Formula", "Compound.name")],
             hmdb_ms1@spectra.info[idx,c("mz", "Formula", "Compound.name")]) %>%
  dplyr::mutate(mz_error = abs(mz - mz.1)*10^6/mz) %>%
  dplyr::filter(is.na(mz_error) | mz_error > 10)

write.csv(temp2, file = "temp2.csv", row.names = FALSE)

temp2 <- readr::read_csv("temp2_manual.csv")

temp2$mz
temp2$Formula
temp2$HMDB.ID

spectra.info$mz[match(temp2$Lab.ID, spectra.info$Lab.ID)] <-
  temp2$mz
spectra.info$Formula[match(temp2$Lab.ID, spectra.info$Lab.ID)] <-
  temp2$Formula

spectra.info$HMDB.ID[match(temp2$Lab.ID, spectra.info$Lab.ID)] <-
  temp2$HMDB.ID

mpsnyder_rplc_ms2@spectra.info <-
  spectra.info

write.csv(spectra.info, file = "spectra.info.csv", row.names = FALSE)

spectra.info <- readr::read_csv("spectra.info_manual.csv")
mpsnyder_rplc_ms2@spectra.info <-
  spectra.info

save(mpsnyder_rplc_ms2, file = "mpsnyder_rplc_ms2.rda")
