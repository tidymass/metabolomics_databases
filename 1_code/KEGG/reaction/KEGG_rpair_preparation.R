no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("2_data/KEGG/reaction/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)

rm(list = ls())

library(tidyverse)

load("../kegg_rpair_database")

load("../kegg_reaction_database")

load("../kegg_reaction_universal_database")

kegg_rpair_database_human <-
  tibble::as_tibble(kegg_rpair_database)

load("../kegg_ms1.rda")

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::rename(
    from_compound_KEGG_ID = from_KEGG.ID,
    to_compound_KEGG_ID = to_KEGG.ID,
    reaction_KEGG_ID = reaction_KEGG
  )

reaction_RHEA_ID <-
  seq_len(nrow(kegg_rpair_database_human)) %>%
  purrr::map(function(i) {
    # cat(i, " ")
    kegg_id <-
      kegg_rpair_database_human$reaction_KEGG_ID[i] %>%
      stringr::str_split("\\{\\}") %>%
      `[[`(1)

    reaction_RHEA_ID <-
      kegg_reaction_universal_database$DBLINKS[match(kegg_id, kegg_reaction_universal_database$KEGG.ID)]

    reaction_RHEA_ID <-
      reaction_RHEA_ID[!is.na(reaction_RHEA_ID)] %>%
      stringr::str_replace_all("RHEA: ", "")

    if (length(reaction_RHEA_ID) > 0) {
      reaction_RHEA_ID <-
        paste(reaction_RHEA_ID, collapse = "{}")
    } else{
      return(NA)
    }

    reaction_RHEA_ID

  }) %>%
  unlist()

kegg_rpair_database_human$reaction_RHEA_ID <-
  reaction_RHEA_ID

reaction_Enzyme_EC_number <-
  seq_len(nrow(kegg_rpair_database_human)) %>%
  purrr::map(function(i) {
    # cat(i, " ")
    reaction_id <-
      stringr::str_split(kegg_rpair_database_human$reaction_KEGG_ID[i], "\\{\\}")[[1]] %>%
      unique()
    kegg_reaction_universal_database$ENZYME[match(reaction_id, kegg_reaction_universal_database$KEGG.ID)] %>%
      unique() %>%
      paste(collapse = "{}")
  }) %>%
  unlist()

kegg_rpair_database_human$reaction_Enzyme_EC_number <-
  reaction_Enzyme_EC_number

spectra.info <-
  kegg_ms1@spectra.info %>%
  tibble::as_tibble() %>%
  dplyr::select(KEGG.ID, Compound.name, HMDB.ID, mz) %>%
  dplyr::distinct(KEGG.ID, .keep_all = TRUE)

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::left_join(spectra.info, by = c("from_compound_KEGG_ID" = 'KEGG.ID')) %>%
  dplyr::rename(
    from_compound_name = Compound.name,
    from_compound_HMDB_ID = HMDB.ID,
    from_compound_mz = mz
  ) %>%
  dplyr::left_join(spectra.info, by = c("to_compound_KEGG_ID" = 'KEGG.ID')) %>%
  dplyr::rename(
    to_compound_name = Compound.name,
    to_compound_HMDB_ID = HMDB.ID,
    to_compound_mz = mz
  )

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::select(
    from_compound_HMDB_ID,
    from_compound_KEGG_ID,
    to_compound_HMDB_ID,
    to_compound_KEGG_ID,
    from_compound_name,
    to_compound_name,
    from_compound_mz,
    to_compound_mz,
    reaction_KEGG_ID,
    reaction_RHEA_ID,
    reaction_Enzyme_EC_number
  )

save(kegg_rpair_database_human, file = "kegg_rpair_database_human")
