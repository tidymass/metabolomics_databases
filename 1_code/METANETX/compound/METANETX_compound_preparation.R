setwd(masstools::get_project_wd())
rm(list = ls())

setwd("2_data/METANETX/compound/")

library(tidyverse)
library(xml2)
library(stringr)
library(massdatabase)
library(tibble)

compound_xref <-
  readr::read_tsv("chem_xref.tsv")

compound_xref <-
  compound_xref %>%
  dplyr::filter(ID != "BIOMASS") %>%
  dplyr::select(ID, `#source`) %>%
  dplyr::rename(compound_METANETX_ID = ID,
                other_id = `#source`) %>%
  dplyr::arrange(compound_METANETX_ID) %>%
  dplyr::filter(compound_METANETX_ID != other_id) %>%
  dplyr::filter(!stringr::str_detect(other_id, "mnx:")) %>%
  dplyr::filter(stringr::str_detect(other_id, "(chebi)|(CHEBI)|(hmdb)|(kegg\\.compound)"))

all_database <-
  compound_xref$other_id %>%
  stringr::str_split(pattern = "\\:") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist() %>%
  unique() %>%
  sort()

sort(all_database)

other_ids <-
  1:length(unique(compound_xref$compound_METANETX_ID)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <- unique(compound_xref$compound_METANETX_ID)[i]

    temp <-
      compound_xref %>%
      dplyr::filter(compound_METANETX_ID == x)

    compound_CHEBI_ID <-
      grep("(chebi)|(CHEBI)", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace_all("(chebi\\:)|(CHEBI\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")

    compound_HMDB_ID <-
      grep("(hmdb)|(hmdb)", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace("(hmdb\\:)|(hmdb\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")

    compound_KEGG_ID <-
      grep("kegg\\.compound", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace("(kegg\\.compound\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")


    data.frame(compound_CHEBI_ID,
               compound_HMDB_ID,
               compound_KEGG_ID)
  }) %>%
  dplyr::bind_rows()

other_ids$compound_METANETX_ID <-
  unique(compound_xref$compound_METANETX_ID)

other_ids <-
  other_ids %>%
  tibble::as_tibble()

other_ids[other_ids == ""] <- NA

metanetx_compound_data <-
  other_ids

save(metanetx_compound_data, file = "metanetx_compound_data")

