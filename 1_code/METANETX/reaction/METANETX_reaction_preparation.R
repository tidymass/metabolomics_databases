setwd(masstools::get_project_wd())
rm(list = ls())

setwd("2_data/METANETX/reaction/")

library(tidyverse)
library(xml2)
library(stringr)
library(massdatabase)
library(tibble)

reaction_xref <-
  readr::read_tsv("reac_xref.tsv") %>%
  dplyr::filter(ID != "EMPTY") %>%
  dplyr::select(ID, `#source`) %>%
  dplyr::rename(reaction_METANETX_ID = ID,
                other_id = `#source`) %>%
  dplyr::arrange(reaction_METANETX_ID) %>%
  dplyr::filter(reaction_METANETX_ID != other_id) %>%
  dplyr::filter(!stringr::str_detect(other_id, "mnx:")) %>%
  dplyr::filter(!stringr::str_detect(other_id, "bigg\\.reaction\\:R_")) %>%
  dplyr::filter(!stringr::str_detect(other_id, "biggR\\:R_"))

all_database <-
  reaction_xref$other_id %>%
  stringr::str_split(pattern = "\\:") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist() %>%
  unique()

other_ids <-
  1:length(unique(reaction_xref$reaction_METANETX_ID)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <- unique(reaction_xref$reaction_METANETX_ID)[i]
    temp <-
      reaction_xref %>%
      dplyr::filter(reaction_METANETX_ID == x)

    reaction_BIGG_ID <-
      grep("(bigg\\.reaction)|(biggR)", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace_all("(bigg\\.reaction\\:)|(biggR\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")

    reaction_METACYC_ID <-
      grep("(metacyc\\.reaction)|(metacycR)", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace_all("(metacyc\\.reaction\\:)|(metacycR\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")

    reaction_RHEA_ID <-
      grep("(rhea)|(rheaR)", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace_all("(rhea\\:)|(rheaR\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")

    reaction_MODELSEED_ID <-
      grep("(seed\\.reaction)|(seedR)", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace_all("(seed\\.reaction\\:)|(seedR\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")

    reaction_KEGG_ID <-
      grep("(kegg\\.reaction)|(keggR)", temp$other_id,
           value = TRUE) %>%
      stringr::str_replace_all("(kegg\\.reaction\\:)|(keggR\\:)", "") %>%
      unique() %>%
      paste(collapse = "{}")

    data.frame(
      reaction_BIGG_ID,
      reaction_METACYC_ID,
      reaction_RHEA_ID,
      reaction_MODELSEED_ID,
      reaction_KEGG_ID
    )
  }) %>%
  dplyr::bind_rows()

other_ids$reaction_METANETX_ID <-
  unique(reaction_xref$reaction_METANETX_ID)

other_ids <-
  other_ids %>%
  tibble::as_tibble()

other_ids[other_ids == ""] <- NA

reaction_prop <-
  readr::read_tsv("reac_prop.tsv")

reaction_prop <-
  reaction_prop[-1,] %>%
  tibble::as_tibble() %>%
  dplyr::select(`#ID`, mnx_equation, reference, classifs) %>%
  dplyr::rename(
    reaction_METANETX_ID = `#ID`,
    reaction_METANETX_equation = mnx_equation,
    reaction_EC_number = classifs
  )


reaction_prop %>% dim()
dim(other_ids)

reaction_prop <-
  reaction_prop %>%
  dplyr::left_join(other_ids,
                   by = "reaction_METANETX_ID")

metanetx_reaction_data <-
  reaction_prop

save(metanetx_reaction_data, file = "metanetx_reaction_data")
