no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
source("R/read_sbml_data.R")
source("R/19_REACTOME.R")

load("other_files/CHEBI/chebi_ms1.rda")

chebi_compound_data <-
  chebi_ms1@spectra.info %>%
  dplyr::select(
    compound_name = Compound.name,
    compound_mz = mz,
    compound_HMDB_ID = HMDB.ID,
    compound_KEGG_ID = KEGG.ID,
    compound_CHEBI_ID = CHEBI.ID
  )

chebi_compound_data$compound_KEGG_ID <-
  chebi_compound_data$compound_KEGG_ID %>%
  stringr::str_split("\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

setwd("other_files/REACTOME/reaction/")

load("reactome_reaction_human_database")

reactome_reaction_human_info <-
  readr::read_delim("reactome_reaction_exporter_v81.txt") %>%
  dplyr::select(reaction_id, reaction_name) %>%
  dplyr::rename(reaction_REACTOME_ID = reaction_id,
                reaction_REACTOME_name = reaction_name) %>%
  dplyr::distinct(.keep_all = TRUE)

library(metid)
library(tidyverse)
library(XML)
library(SBMLR)

reactome_rpair_database_human <-
  vector(mode = "list",
         length = length(reactome_reaction_human_database))

for (i in 1:length(reactome_reaction_human_database)) {
  # cat(i, " ")
  if (is.null(reactome_reaction_human_database[[i]])) {
    reactome_rpair_database_human[[i]] <- NULL
    next()
  }

  x <-
    reactome_reaction_human_database[[i]] %>%
    purrr::map(function(y) {
      y %>%
        dplyr::filter(!is.na(compound_CHEBI_ID))
    })

  if (nrow(x[[1]]) == 0 | nrow(x[[2]]) == 0) {
    reactome_rpair_database_human[[i]] <- NULL
    next()
  }

  x <-
    x %>%
    purrr::map(function(y) {
      # y <-
      y %>%
        dplyr::left_join(chebi_compound_data,
                         by = c("compound_CHEBI_ID" = "compound_CHEBI_ID")) %>%
        dplyr::filter(!is.na(compound_HMDB_ID) |
                        !is.na(compound_KEGG_ID)) %>%
        dplyr::select(-c(species, name)) %>%
        dplyr::filter(compound_mz > 18.01056) %>%
        dplyr::filter(
          !compound_name %in% c("calcium(2+)", "dioxygen", "carbon dioxide", "nitric oxide")
        )
    })

  if (nrow(x[[1]]) == 0 | nrow(x[[2]]) == 0) {
    reactome_rpair_database_human[[i]] <- NULL
    next()
  }

  if ((i / 100) %in% seq(1, 100000, 1)) {
    cat(i, " ")
  }

  reactome_rpair_database_human[[i]] <-
    purrr::map(seq_len(nrow(x[[1]])), function(index1) {
      purrr::map(seq_len(nrow(x[[2]])), function(index2) {
        cbind(
          x[[1]][index1,] %>%
            dplyr::rename(
              from_compound_CHEBI_ID = compound_CHEBI_ID,
              from_compound_REACTOME_ID = compound_REACTOME_ID,
              from_compound_name = compound_name,
              from_compound_mz = compound_mz,
              from_compound_HMDB_ID = compound_HMDB_ID,
              from_compound_KEGG_ID = compound_KEGG_ID
            ),
          x[[2]][index2,] %>%
            dplyr::rename(
              to_compound_CHEBI_ID = compound_CHEBI_ID,
              to_compound_REACTOME_ID = compound_REACTOME_ID,
              to_compound_name = compound_name,
              to_compound_mz = compound_mz,
              to_compound_HMDB_ID = compound_HMDB_ID,
              to_compound_KEGG_ID = compound_KEGG_ID
            )
        )
      }) %>%
        dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows()

  reactome_rpair_database_human[[i]] <-
    data.frame(reactome_rpair_database_human[[i]],
               reactome_reaction_human_info[i, ])

}


reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  dplyr::bind_rows()

dim(reactome_rpair_database_human)


#####
library(plyr)

temp <-
  reactome_rpair_database_human %>%
  plyr::dlply(.variables = .(from_compound_REACTOME_ID, to_compound_REACTOME_ID))

temp <-
  temp %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }

    x$reaction_REACTOME_ID <-
      x$reaction_REACTOME_ID[!is.na(x$reaction_REACTOME_ID)] %>%
      paste(collapse = "{}")

    x$reaction_REACTOME_name <-
      x$reaction_REACTOME_name[!is.na(x$reaction_REACTOME_name)] %>%
      paste(collapse = "{}")

    x %>%
      dplyr::distinct(from_compound_REACTOME_ID,
                      to_compound_REACTOME_ID,
                      .keep_all = TRUE)

  }) %>%
  dplyr::bind_rows()

temp <-
  temp %>%
  dplyr::filter(from_compound_REACTOME_ID != to_compound_REACTOME_ID)

rpair_id <-
  seq_len(nrow(temp)) %>%
  purrr::map(function(i) {
    c(temp$from_compound_REACTOME_ID[i],
      temp$to_compound_REACTOME_ID[i]) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

temp <-
  temp %>%
  dplyr::mutate(rpair_id) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE) %>%
  dplyr::select(-rpair_id)

reactome_rpair_database_human <-
  temp %>%
  tibble::as_tibble()

save(reactome_rpair_database_human, file = "reactome_rpair_database_human")
