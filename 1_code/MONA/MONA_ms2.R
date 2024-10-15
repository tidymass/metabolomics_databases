###source
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
load("2_data/HMDB/MS1/hmdb_ms1.rda")
load("2_data/KEGG/kegg_ms1.rda")
source("R/read_msp_data.R")
source("1_code/3_utils.R")
source("R/11_MONA.R")
setwd("2_data/MONA")

library(tidyverse)
library(tidymass)

# # #####CASMI2016
# # data1 <-
# #   read_msp_data(file = "MoNA-export-CASMI_2016.msp",
# #                 source = "mona",
# #                 threads = 5)
# #
# # mona_casmi2016_ms2 <-
# #   convert_mona2metid(data = data1,
# #                      path = "CASMI2016",
# #                      threads = 5)
# #
# # new_spectra_info <-
# #   mona_casmi2016_ms2@spectra.info$Comments %>%
# #   purrr::map(function(x) {
# #     x <-
# #       x %>%
# #       stringr::str_split(pattern = '\"') %>%
# #       `[[`(1) %>%
# #       stringr::str_trim()
# #     x <- x[x != ""]
# #     x <- x[x != " "]
# #     SMILES.ID =
# #       x[stringr::str_detect(x, "^SMILES")] %>% stringr::str_replace("SMILES=", "")
# #     SMILES.ID <- ifelse(length(SMILES.ID) == 0, NA, SMILES.ID)
# #     CAS.ID = x[stringr::str_detect(x, "^cas")] %>% stringr::str_replace("cas=", "")
# #     CAS.ID <- ifelse(length(CAS.ID) == 0, NA, CAS.ID)
# #     PUBCHEM.ID = x[stringr::str_detect(x, "^pubchem cid=")] %>% stringr::str_replace("^pubchem cid=", "")
# #     PUBCHEM.ID <- ifelse(length(PUBCHEM.ID) == 0, NA, PUBCHEM.ID)
# #     CHEMSPIDER.ID = x[stringr::str_detect(x, "^chemspider=")] %>% stringr::str_replace("^chemspider=", "")
# #     CHEMSPIDER.ID <-
# #       ifelse(length(CHEMSPIDER.ID) == 0, NA, CHEMSPIDER.ID)
# #     Author = x[stringr::str_detect(x, "^author=")] %>% stringr::str_replace("^author=", "")
# #     Author <- ifelse(length(Author) == 0, NA, Author)
# #     Comment_confidence <-
# #       x[stringr::str_detect(x, "^comment=CONFIDENCE")] %>% stringr::str_replace("^comment=", "")
# #     Comment_confidence <-
# #       ifelse(length(Comment_confidence) == 0, NA, Comment_confidence)
# #     Submitter_team <-
# #       x[stringr::str_detect(x, "^submitter=")] %>% stringr::str_replace("^submitter=", "")
# #     Submitter_team <-
# #       ifelse(length(Submitter_team) == 0, NA, Submitter_team)
# #     data.frame(
# #       SMILES.ID = SMILES.ID,
# #       CAS.ID = CAS.ID,
# #       PUBCHEM.ID = PUBCHEM.ID,
# #       CHEMSPIDER.ID = CHEMSPIDER.ID,
# #       Author = Author,
# #       Comment_confidence = Comment_confidence,
# #       Submitter_team = Submitter_team
# #     )
# #   }) %>%
# #   dplyr::bind_rows() %>%
# #   as.data.frame()
# #
# # colnames(mona_casmi2016_ms2@spectra.info)
# # colnames(new_spectra_info)
# # mona_casmi2016_ms2@spectra.info$SMILES.ID <-
# #   new_spectra_info$SMILES.ID
# # mona_casmi2016_ms2@spectra.info$CAS.ID <- new_spectra_info$CAS.ID
# # mona_casmi2016_ms2@spectra.info$PUBCHEM.ID <-
# #   new_spectra_info$PUBCHEM.ID
# # mona_casmi2016_ms2@spectra.info$CHEMSPIDER.ID <-
# #   new_spectra_info$CHEMSPIDER.ID
# # mona_casmi2016_ms2@spectra.info$Author <- new_spectra_info$Author
# # mona_casmi2016_ms2@spectra.info$Comment_confidence <-
# #   new_spectra_info$Comment_confidence
# # mona_casmi2016_ms2@spectra.info$Submitter_team <-
# #   new_spectra_info$Submitter_team
# #
# # idx1 <-
# #   match(
# #     mona_casmi2016_ms2@spectra.info$INCHI.ID,
# #     hmdb_ms1@spectra.info$INCHI.ID,
# #     incomparables = NA
# #   )
# #
# # idx2 <-
# #   match(
# #     mona_casmi2016_ms2@spectra.info$INCHIKEY.ID,
# #     hmdb_ms1@spectra.info$INCHIKEY.ID,
# #     incomparables = NA
# #   )
# #
# # idx3 <-
# #   match(
# #     mona_casmi2016_ms2@spectra.info$CAS.ID,
# #     hmdb_ms1@spectra.info$CAS.ID,
# #     incomparables = NA
# #   )
# #
# # idx4 <-
# #   match(
# #     mona_casmi2016_ms2@spectra.info$SMILES.ID,
# #     hmdb_ms1@spectra.info$SMILES.ID,
# #     incomparables = NA
# #   )
# #
# # idx <-
# #   data.frame(idx1, idx2, idx3, idx4) %>%
# #   apply(1, function(x) {
# #     x <- as.numeric(x)
# #     x <- x[!is.na(x)]
# #     if (length(x) == 0) {
# #       return(NA)
# #     }
# #     return(x[1])
# #   })
# #
# # colnames(mona_casmi2016_ms2@spectra.info)
# # colnames(hmdb_ms1@spectra.info)
# #
# # spectra.info <-
# #   mona_casmi2016_ms2@spectra.info
# #
# # spectra.info <-
# #   spectra.info %>%
# #   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID))
# #
# # spectra.info <-
# #   data.frame(spectra.info,
# #              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
# #   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
# #
# #
# # mona_casmi2016_ms2@spectra.info <-
# #   spectra.info
# #
# # save(mona_casmi2016_ms2, file = "mona_casmi2016_ms2.rda")
# #
# #
# #
# #
# #
# #
# # #####CASMI2012
# # data1 <-
# #   read_msp_data(file = "MoNA-export-CASMI_2012.msp",
# #                 source = "mona",
# #                 threads = 5)
# #
# # mona_casmi2012_ms2 <-
# #   convert_mona2metid(data = data1,
# #                      path = "CASMI2012",
# #                      threads = 5)
# #
# # new_spectra_info <-
# #   mona_casmi2012_ms2@spectra.info$Comments %>%
# #   purrr::map(function(x) {
# #     x <-
# #       x %>%
# #       stringr::str_split(pattern = '\"') %>%
# #       `[[`(1) %>%
# #       stringr::str_trim()
# #     x <- x[x != ""]
# #     x <- x[x != " "]
# #     SMILES.ID =
# #       x[stringr::str_detect(x, "^SMILES")] %>% stringr::str_replace("SMILES=", "")
# #     SMILES.ID <- ifelse(length(SMILES.ID) == 0, NA, SMILES.ID)
# #     CAS.ID = x[stringr::str_detect(x, "^cas")] %>% stringr::str_replace("cas=", "")
# #     CAS.ID <- ifelse(length(CAS.ID) == 0, NA, CAS.ID)
# #     PUBCHEM.ID = x[stringr::str_detect(x, "^pubchem cid=")] %>% stringr::str_replace("^pubchem cid=", "")
# #     PUBCHEM.ID <- ifelse(length(PUBCHEM.ID) == 0, NA, PUBCHEM.ID)
# #     CHEMSPIDER.ID = x[stringr::str_detect(x, "^chemspider=")] %>% stringr::str_replace("^chemspider=", "")
# #     CHEMSPIDER.ID <-
# #       ifelse(length(CHEMSPIDER.ID) == 0, NA, CHEMSPIDER.ID)
# #     Author = x[stringr::str_detect(x, "^author=")] %>% stringr::str_replace("^author=", "")
# #     Author <- ifelse(length(Author) == 0, NA, Author)
# #     Comment_confidence <-
# #       x[stringr::str_detect(x, "^comment=CONFIDENCE")] %>% stringr::str_replace("^comment=", "")
# #     Comment_confidence <-
# #       ifelse(length(Comment_confidence) == 0, NA, Comment_confidence)
# #     Submitter_team <-
# #       x[stringr::str_detect(x, "^submitter=")] %>% stringr::str_replace("^submitter=", "")
# #     Submitter_team <-
# #       ifelse(length(Submitter_team) == 0, NA, Submitter_team)
# #     data.frame(
# #       SMILES.ID = SMILES.ID,
# #       CAS.ID = CAS.ID,
# #       PUBCHEM.ID = PUBCHEM.ID,
# #       CHEMSPIDER.ID = CHEMSPIDER.ID,
# #       Author = Author,
# #       Comment_confidence = Comment_confidence,
# #       Submitter_team = Submitter_team
# #     )
# #   }) %>%
# #   dplyr::bind_rows() %>%
# #   as.data.frame()
# #
# # colnames(mona_casmi2012_ms2@spectra.info)
# # colnames(new_spectra_info)
# # mona_casmi2012_ms2@spectra.info$SMILES.ID <-
# #   new_spectra_info$SMILES.ID
# # mona_casmi2012_ms2@spectra.info$CAS.ID <- new_spectra_info$CAS.ID
# # mona_casmi2012_ms2@spectra.info$PUBCHEM.ID <-
# #   new_spectra_info$PUBCHEM.ID
# # mona_casmi2012_ms2@spectra.info$CHEMSPIDER.ID <-
# #   new_spectra_info$CHEMSPIDER.ID
# # mona_casmi2012_ms2@spectra.info$Author <- new_spectra_info$Author
# # mona_casmi2012_ms2@spectra.info$Comment_confidence <-
# #   new_spectra_info$Comment_confidence
# # mona_casmi2012_ms2@spectra.info$Submitter_team <-
# #   new_spectra_info$Submitter_team
# #
# # idx1 <-
# #   match(
# #     mona_casmi2012_ms2@spectra.info$INCHI.ID,
# #     hmdb_ms1@spectra.info$INCHI.ID,
# #     incomparables = NA
# #   )
# #
# # idx2 <-
# #   match(
# #     mona_casmi2012_ms2@spectra.info$INCHIKEY.ID,
# #     hmdb_ms1@spectra.info$INCHIKEY.ID,
# #     incomparables = NA
# #   )
# #
# # idx3 <-
# #   match(
# #     mona_casmi2012_ms2@spectra.info$CAS.ID,
# #     hmdb_ms1@spectra.info$CAS.ID,
# #     incomparables = NA
# #   )
# #
# # idx4 <-
# #   match(
# #     mona_casmi2012_ms2@spectra.info$SMILES.ID,
# #     hmdb_ms1@spectra.info$SMILES.ID,
# #     incomparables = NA
# #   )
# #
# # idx <-
# #   data.frame(idx1, idx2, idx3, idx4) %>%
# #   apply(1, function(x) {
# #     x <- as.numeric(x)
# #     x <- x[!is.na(x)]
# #     if (length(x) == 0) {
# #       return(NA)
# #     }
# #     return(x[1])
# #   })
# #
# # colnames(mona_casmi2012_ms2@spectra.info)
# # colnames(hmdb_ms1@spectra.info)
# #
# # spectra.info <-
# #   mona_casmi2012_ms2@spectra.info
# #
# # spectra.info <-
# #   spectra.info %>%
# #   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID))
# #
# # spectra.info <-
# #   data.frame(spectra.info,
# #              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
# #   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
# #
# #
# # mona_casmi2012_ms2@spectra.info <-
# #   spectra.info
# #
# # save(mona_casmi2012_ms2, file = "mona_casmi2012_ms2.rda")
# #
# #
# #
# #
# #
# #
# #
# #
# # #####Fiehn HILIC
# # data1 <-
# #   read_msp_data(file = "MoNA-export-Fiehn_HILIC.msp",
# #                 source = "mona",
# #                 threads = 5)
# #
# # mona_fiehnhilic_ms2 <-
# #   convert_mona2metid(data = data1,
# #                      path = "Fiehn_HILIC",
# #                      threads = 5)
# #
# # new_spectra_info <-
# #   mona_fiehnhilic_ms2@spectra.info$Comments %>%
# #   purrr::map(function(x) {
# #     x <-
# #       x %>%
# #       stringr::str_split(pattern = '\"') %>%
# #       `[[`(1) %>%
# #       stringr::str_trim()
# #     x <- x[x != ""]
# #     x <- x[x != " "]
# #     SMILES.ID =
# #       x[stringr::str_detect(x, "^SMILES")] %>% stringr::str_replace("SMILES=", "")
# #     SMILES.ID <- ifelse(length(SMILES.ID) == 0, NA, SMILES.ID)
# #     CAS.ID = x[stringr::str_detect(x, "^cas")] %>% stringr::str_replace("cas=", "")
# #     CAS.ID <- ifelse(length(CAS.ID) == 0, NA, CAS.ID)
# #     PUBCHEM.ID = x[stringr::str_detect(x, "^pubchem cid=")] %>% stringr::str_replace("^pubchem cid=", "")
# #     PUBCHEM.ID <- ifelse(length(PUBCHEM.ID) == 0, NA, PUBCHEM.ID)
# #     CHEMSPIDER.ID = x[stringr::str_detect(x, "^chemspider=")] %>% stringr::str_replace("^chemspider=", "")
# #     CHEMSPIDER.ID <-
# #       ifelse(length(CHEMSPIDER.ID) == 0, NA, CHEMSPIDER.ID)
# #     Author = x[stringr::str_detect(x, "^author=")] %>% stringr::str_replace("^author=", "")
# #     Author <- ifelse(length(Author) == 0, NA, Author)
# #     Comment_confidence <-
# #       x[stringr::str_detect(x, "^comment=CONFIDENCE")] %>% stringr::str_replace("^comment=", "")
# #     Comment_confidence <-
# #       ifelse(length(Comment_confidence) == 0, NA, Comment_confidence)
# #     Submitter_team <-
# #       x[stringr::str_detect(x, "^submitter=")] %>% stringr::str_replace("^submitter=", "")
# #     Submitter_team <-
# #       ifelse(length(Submitter_team) == 0, NA, Submitter_team)
# #     data.frame(
# #       SMILES.ID = SMILES.ID,
# #       CAS.ID = CAS.ID,
# #       PUBCHEM.ID = PUBCHEM.ID,
# #       CHEMSPIDER.ID = CHEMSPIDER.ID,
# #       Author = Author,
# #       Comment_confidence = Comment_confidence,
# #       Submitter_team = Submitter_team
# #     )
# #   }) %>%
# #   dplyr::bind_rows() %>%
# #   as.data.frame()
# #
# # colnames(mona_fiehnhilic_ms2@spectra.info)
# # colnames(new_spectra_info)
# # mona_fiehnhilic_ms2@spectra.info$SMILES.ID <-
# #   new_spectra_info$SMILES.ID
# # mona_fiehnhilic_ms2@spectra.info$CAS.ID <- new_spectra_info$CAS.ID
# # mona_fiehnhilic_ms2@spectra.info$PUBCHEM.ID <-
# #   new_spectra_info$PUBCHEM.ID
# # mona_fiehnhilic_ms2@spectra.info$CHEMSPIDER.ID <-
# #   new_spectra_info$CHEMSPIDER.ID
# # mona_fiehnhilic_ms2@spectra.info$Author <- new_spectra_info$Author
# # mona_fiehnhilic_ms2@spectra.info$Comment_confidence <-
# #   new_spectra_info$Comment_confidence
# # mona_fiehnhilic_ms2@spectra.info$Submitter_team <-
# #   new_spectra_info$Submitter_team
# #
# # idx1 <-
# #   match(
# #     mona_fiehnhilic_ms2@spectra.info$INCHI.ID,
# #     hmdb_ms1@spectra.info$INCHI.ID,
# #     incomparables = NA
# #   )
# #
# # idx2 <-
# #   match(
# #     mona_fiehnhilic_ms2@spectra.info$INCHIKEY.ID,
# #     hmdb_ms1@spectra.info$INCHIKEY.ID,
# #     incomparables = NA
# #   )
# #
# # idx3 <-
# #   match(
# #     mona_fiehnhilic_ms2@spectra.info$CAS.ID,
# #     hmdb_ms1@spectra.info$CAS.ID,
# #     incomparables = NA
# #   )
# #
# # idx4 <-
# #   match(
# #     mona_fiehnhilic_ms2@spectra.info$SMILES.ID,
# #     hmdb_ms1@spectra.info$SMILES.ID,
# #     incomparables = NA
# #   )
# #
# # idx <-
# #   data.frame(idx1, idx2, idx3, idx4) %>%
# #   apply(1, function(x) {
# #     x <- as.numeric(x)
# #     x <- x[!is.na(x)]
# #     if (length(x) == 0) {
# #       return(NA)
# #     }
# #     return(x[1])
# #   })
# #
# # colnames(mona_fiehnhilic_ms2@spectra.info)
# # colnames(hmdb_ms1@spectra.info)
# #
# # spectra.info <-
# #   mona_fiehnhilic_ms2@spectra.info
# #
# # spectra.info <-
# #   spectra.info %>%
# #   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID))
# #
# # spectra.info <-
# #   data.frame(spectra.info,
# #              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
# #   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
# #
# # mona_fiehnhilic_ms2@spectra.info <-
# #   spectra.info
# #
# # save(mona_fiehnhilic_ms2, file = "mona_fiehnhilic_ms2.rda")
# #
# #
# #
# #
# #
# #
# #
# #
# # #####EMBL-MCF
# # data1 <-
# #   read_msp_data(file = "MoNA-export-EMBL-MCF.msp",
# #                 source = "mona",
# #                 threads = 5)
# #
# # mona_emblmcf_ms2 <-
# #   convert_mona2metid(data = data1,
# #                      path = "EMBL-MCF",
# #                      threads = 5)
# #
# # new_spectra_info <-
# #   mona_emblmcf_ms2@spectra.info$Comments %>%
# #   purrr::map(function(x) {
# #     x <-
# #       x %>%
# #       stringr::str_split(pattern = '\"') %>%
# #       `[[`(1) %>%
# #       stringr::str_trim()
# #     x <- x[x != ""]
# #     x <- x[x != " "]
# #     SMILES.ID =
# #       x[stringr::str_detect(x, "^SMILES")] %>% stringr::str_replace("SMILES=", "")
# #     SMILES.ID <- ifelse(length(SMILES.ID) == 0, NA, SMILES.ID)
# #     CAS.ID = x[stringr::str_detect(x, "^cas")] %>% stringr::str_replace("cas=", "")
# #     CAS.ID <- ifelse(length(CAS.ID) == 0, NA, CAS.ID)
# #     PUBCHEM.ID = x[stringr::str_detect(x, "^pubchem cid=")] %>% stringr::str_replace("^pubchem cid=", "")
# #     PUBCHEM.ID <- ifelse(length(PUBCHEM.ID) == 0, NA, PUBCHEM.ID)
# #     CHEMSPIDER.ID = x[stringr::str_detect(x, "^chemspider=")] %>% stringr::str_replace("^chemspider=", "")
# #     CHEMSPIDER.ID <-
# #       ifelse(length(CHEMSPIDER.ID) == 0, NA, CHEMSPIDER.ID)
# #     Author = x[stringr::str_detect(x, "^author=")] %>% stringr::str_replace("^author=", "")
# #     Author <- ifelse(length(Author) == 0, NA, Author)
# #     Comment_confidence <-
# #       x[stringr::str_detect(x, "^comment=CONFIDENCE")] %>% stringr::str_replace("^comment=", "")
# #     Comment_confidence <-
# #       ifelse(length(Comment_confidence) == 0, NA, Comment_confidence)
# #     Submitter_team <-
# #       x[stringr::str_detect(x, "^submitter=")] %>% stringr::str_replace("^submitter=", "")
# #     Submitter_team <-
# #       ifelse(length(Submitter_team) == 0, NA, Submitter_team)
# #     data.frame(
# #       SMILES.ID = SMILES.ID,
# #       CAS.ID = CAS.ID,
# #       PUBCHEM.ID = PUBCHEM.ID,
# #       CHEMSPIDER.ID = CHEMSPIDER.ID,
# #       Author = Author,
# #       Comment_confidence = Comment_confidence,
# #       Submitter_team = Submitter_team
# #     )
# #   }) %>%
# #   dplyr::bind_rows() %>%
# #   as.data.frame()
# #
# # colnames(mona_emblmcf_ms2@spectra.info)
# # colnames(new_spectra_info)
# # mona_emblmcf_ms2@spectra.info$SMILES.ID <-
# #   new_spectra_info$SMILES.ID
# # mona_emblmcf_ms2@spectra.info$CAS.ID <-
# #   new_spectra_info$CAS.ID
# # mona_emblmcf_ms2@spectra.info$PUBCHEM.ID <-
# #   new_spectra_info$PUBCHEM.ID
# # mona_emblmcf_ms2@spectra.info$CHEMSPIDER.ID <-
# #   new_spectra_info$CHEMSPIDER.ID
# # mona_emblmcf_ms2@spectra.info$Author <-
# #   new_spectra_info$Author
# # mona_emblmcf_ms2@spectra.info$Comment_confidence <-
# #   new_spectra_info$Comment_confidence
# # mona_emblmcf_ms2@spectra.info$Submitter_team <-
# #   new_spectra_info$Submitter_team
# #
# # idx1 <-
# #   match(
# #     mona_emblmcf_ms2@spectra.info$INCHI.ID,
# #     hmdb_ms1@spectra.info$INCHI.ID,
# #     incomparables = NA
# #   )
# #
# # idx2 <-
# #   match(
# #     mona_emblmcf_ms2@spectra.info$INCHIKEY.ID,
# #     hmdb_ms1@spectra.info$INCHIKEY.ID,
# #     incomparables = NA
# #   )
# #
# # idx3 <-
# #   match(
# #     mona_emblmcf_ms2@spectra.info$CAS.ID,
# #     hmdb_ms1@spectra.info$CAS.ID,
# #     incomparables = NA
# #   )
# #
# # idx4 <-
# #   match(
# #     mona_emblmcf_ms2@spectra.info$SMILES.ID,
# #     hmdb_ms1@spectra.info$SMILES.ID,
# #     incomparables = NA
# #   )
# #
# # idx <-
# #   data.frame(idx1, idx2, idx3, idx4) %>%
# #   apply(1, function(x) {
# #     x <- as.numeric(x)
# #     x <- x[!is.na(x)]
# #     if (length(x) == 0) {
# #       return(NA)
# #     }
# #     return(x[1])
# #   })
# #
# # colnames(mona_emblmcf_ms2@spectra.info)
# # colnames(hmdb_ms1@spectra.info)
# #
# # spectra.info <-
# #   mona_emblmcf_ms2@spectra.info
# #
# # spectra.info <-
# #   spectra.info %>%
# #   dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID))
# #
# # spectra.info <-
# #   data.frame(spectra.info,
# #              hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
# #   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))
# #
# # mona_emblmcf_ms2@spectra.info <-
# #   spectra.info
# #
# # save(mona_emblmcf_ms2, file = "mona_emblmcf_ms2.rda")
#
#
# load("mona_casmi2012_ms2.rda")
# load("mona_casmi2016_ms2.rda")
# load("mona_emblmcf_ms2.rda")
# load("mona_fiehnhilic_ms2.rda")
#
#
# mona_casmi2012_ms2
# mona_casmi2016_ms2
#
# intersect(
#   mona_casmi2012_ms2@spectra.info$Compound.name,
#   mona_casmi2016_ms2@spectra.info$Compound.name
# )
#
# colnames(mona_casmi2012_ms2@spectra.info) ==
# colnames(mona_casmi2016_ms2@spectra.info)
#
# mona_casmi2012_ms2@spectra.info$Lab.ID
# mona_casmi2016_ms2@spectra.info$Lab.ID
#
# names(mona_casmi2012_ms2@spectra.data$Spectra.positive)
# names(mona_casmi2012_ms2@spectra.data$Spectra.negative)
#
# mona_ms2 <-
#   mona_casmi2012_ms2
#
# mona_casmi2012_ms2@spectra.info$Lab.ID
# mona_casmi2016_ms2@spectra.info$Lab.ID
# mona_emblmcf_ms2@spectra.info$Lab.ID
# mona_fiehnhilic_ms2@spectra.info$Lab.ID
#
# mona_ms2@spectra.info <-
#   rbind(mona_casmi2012_ms2@spectra.info,
#         mona_casmi2016_ms2@spectra.info,
#         mona_emblmcf_ms2@spectra.info,
#         mona_fiehnhilic_ms2@spectra.info)
#
# mona_ms2@spectra.data$Spectra.positive <-
#   c(mona_casmi2012_ms2@spectra.data$Spectra.positive,
#     mona_casmi2016_ms2@spectra.data$Spectra.positive,
#     mona_emblmcf_ms2@spectra.data$Spectra.positive,
#     mona_fiehnhilic_ms2@spectra.data$Spectra.positive)
#
# mona_ms2@spectra.data$Spectra.negative <-
#   c(mona_casmi2012_ms2@spectra.data$Spectra.negative,
#     mona_casmi2016_ms2@spectra.data$Spectra.negative,
#     mona_emblmcf_ms2@spectra.data$Spectra.negative,
#     mona_fiehnhilic_ms2@spectra.data$Spectra.negative)
#
# ####remove labeled compounds
# mona_ms2@spectra.info$Submitter_team
#
# mona_ms2@spectra.info$Compound.name
#
# grep("\\&", mona_ms2@spectra.info$Compound.name, value = TRUE)
# mona_ms2@spectra.info$INCHIKEY.ID[grep("\\&", mona_ms2@spectra.info$Compound.name, value = FALSE)]
#
# mona_ms2@spectra.info$Compound.name <-
# mona_ms2@spectra.info$Compound.name %>%
#   stringr::str_replace_all("&\\#39\\;", "'")
#
# save(mona_ms2, file = "mona_ms2.rda")
load("mona_ms2.rda")

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

intersect(colnames(mona_ms2@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(mona_ms2@spectra.info))

mona_ms2 <-
  update_metid_database_info(
    database = mona_ms2,
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

intersect(colnames(mona_ms2@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(mona_ms2@spectra.info))

source(here::here("1_code/3_utils.R"))

mona_ms2 <-
  update_metid_database_info(
    database = mona_ms2,
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

mona_ms2 <-
  update_metid_database_source_system(
    database = mona_ms2,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(mona_ms2, file = "mona_ms2.rda")







