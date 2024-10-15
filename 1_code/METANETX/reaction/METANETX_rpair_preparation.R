no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

setwd("2_data/METANETX/reaction/")

library(tidyverse)
library(xml2)
library(stringr)
library(massdatabase)
library(tibble)

load("metanetx_reaction_data")
load("../compound/metanetx_compound_data")

metanetx_compound_data$compound_HMDB_ID <-
  metanetx_compound_data$compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

metanetx_compound_data$compound_KEGG_ID <-
  metanetx_compound_data$compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

metanetx_rpair_database <-
  vector(mode = "list", length = nrow(metanetx_reaction_data))

for (i in 1:nrow(metanetx_reaction_data)) {
  # if (!is.na(metanetx_reaction_data$reaction_KEGG_ID[i])) {
  #   metanetx_rpair_database[[i]] <- NULL
  #   next()
  # }

  x <-
    stringr::str_split(metanetx_reaction_data[i, ]$reaction_METANETX_equation, " = ")[[1]] %>%
    stringr::str_replace_all("\\@[a-zA-Z1-9]{4,6}", "")

  if (x[1] == x[2]) {
    metanetx_rpair_database[[i]] <- NULL
    next()
  }

  x <-
    x %>%
    stringr::str_split(" \\+ ") %>%
    purrr::map(function(y) {
      y <-
        stringr::str_replace(y, "^[0-9]{1,2} ", "") %>%
        stringr::str_trim()
      metanetx_compound_data[match(y, metanetx_compound_data$compound_METANETX_ID), ] %>%
        dplyr::filter(!is.na(compound_HMDB_ID) |
                        !is.na(compound_KEGG_ID))
    })

  if (nrow(x[[1]]) == 0 | nrow(x[[2]]) == 0) {
    metanetx_rpair_database[[i]] <- NULL
    next()
  }

  metanetx_rpair_database[[i]] <-
    purrr::map(seq_len(nrow(x[[1]])), function(index1) {
      purrr::map(seq_len(nrow(x[[2]])), function(index2) {
        cbind(
          x[[1]][index1, ] %>%
            dplyr::rename(
              from_compound_CHEBI_ID = compound_CHEBI_ID,
              from_compound_HMDB_ID = compound_HMDB_ID,
              from_compound_KEGG_ID = compound_KEGG_ID,
              from_compound_METANETX_ID = compound_METANETX_ID
            ),
          x[[2]][index2, ] %>%
            dplyr::rename(
              to_compound_CHEBI_ID = compound_CHEBI_ID,
              to_compound_HMDB_ID = compound_HMDB_ID,
              to_compound_KEGG_ID = compound_KEGG_ID,
              to_compound_METANETX_ID = compound_METANETX_ID
            )
        )
      }) %>%
        dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows()

  metanetx_rpair_database[[i]] <-
    data.frame(metanetx_rpair_database[[i]],
               metanetx_reaction_data[i,])
  if ((i / 10) %in% seq(1, 1000000, 1)) {
    cat(i, " ")
  }

}

metanetx_rpair_database[[1]]

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::bind_rows() %>%
  tibble::as_tibble()

metanetx_rpair_database$from_compound_HMDB_ID
metanetx_rpair_database$to_compound_HMDB_ID

####remove some compounds
# temp <-
#   rbind(
#     metanetx_rpair_database[,c("from_compound_name")] %>%
#       dplyr::rename(name = from_compound_name)
#   )

library(plyr)

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  plyr::dlply(.variables = .(from_compound_METANETX_ID, to_compound_METANETX_ID)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }

    x$reaction_METANETX_ID <-
      x$reaction_METANETX_ID[!is.na(x$reaction_METANETX_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_METANETX_equation <-
      x$reaction_METANETX_equation[!is.na(x$reaction_METANETX_equation)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_EC_number <-
      x$reaction_EC_number[!is.na(x$reaction_EC_number)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_BIGG_ID <-
      x$reaction_BIGG_ID[!is.na(x$reaction_BIGG_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_METACYC_ID <-
      x$reaction_METACYC_ID[!is.na(x$reaction_METACYC_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_RHEA_ID <-
      x$reaction_RHEA_ID[!is.na(x$reaction_RHEA_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_MODELSEED_ID <-
      x$reaction_MODELSEED_ID[!is.na(x$reaction_MODELSEED_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_KEGG_ID <-
      x$reaction_KEGG_ID[!is.na(x$reaction_KEGG_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x <-
      x %>%
      dplyr::distinct(from_compound_METANETX_ID,
                      to_compound_METANETX_ID,
                      .keep_all = TRUE)

    x
  })

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::bind_rows() %>%
  tibble::as_tibble()

###remove some duplicated reaction pair
rpair_id <-
  seq_len(nrow(metanetx_rpair_database)) %>%
  purrr::map(function(i) {
    c(
      metanetx_rpair_database$from_compound_METANETX_ID[i],
      metanetx_rpair_database$to_compound_METANETX_ID[i]
    ) %>%
      sort() %>%
      paste(collapse = "{}")
  }) %>%
  unlist()


metanetx_rpair_database <-
metanetx_rpair_database %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE) %>%
  dplyr::select(-rpair_id)

save(metanetx_rpair_database,
     file = "metanetx_rpair_database")
