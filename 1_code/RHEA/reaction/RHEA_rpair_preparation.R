no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

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

setwd("other_files/RHEA/reaction/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)

load("rhea_reaction")

# equation <-
#   rhea_reaction$Equation %>%
#   stringr::str_split("(\\=)|( \\+ )")
#
# CHEBI_name <-
#   rhea_reaction$CHEBI_name %>%
#   stringr::str_split("\\{\\}")
#
# CHEBI_ID <-
#   rhea_reaction$CHEBI_ID %>%
#   stringr::str_split("\\{\\}")
#
# len1 <-
#   equation %>%
#   lapply(length) %>%
#   unlist()
#
# len2 <-
#   CHEBI_name %>%
#   lapply(length) %>%
#   unlist()
#
# len3 <-
#   CHEBI_ID %>%
#   lapply(length) %>%
#   unlist()

new_equation <-
  seq_len(nrow(rhea_reaction)) %>%
  purrr::map(function(i) {
    if((i/10) %in% seq(1, 100000, 1)){
      cat(i, " ")
    }

    Equation <-
      rhea_reaction$Equation[i]

    CHEBI_name <-
      rhea_reaction$CHEBI_name[i]

    CHEBI_ID <-
      rhea_reaction$CHEBI_ID[i]

    chebi_name_id <-
      data.frame(
        name = stringr::str_split(CHEBI_name, "\\{\\}")[[1]],
        id = stringr::str_split(CHEBI_ID, "\\{\\}")[[1]]
      )

    new_equation <-
      stringr::str_split(Equation, " = ")[[1]]

    new_equation <-
      new_equation %>%
      purrr::map(function(x) {
        x <-
          stringr::str_replace_all(x, "\\(in\\)", "") %>%
          stringr::str_replace_all("\\(out\\)", "")
        x <-
          stringr::str_split(x, " \\+ ")[[1]] %>%
          stringr::str_replace_all("^[0-9]{1,2} ", "") %>%
          stringr::str_trim()
        x <-
          chebi_name_id$id[match(x, chebi_name_id$name)]
        x <-
          x[!is.na(x)]
        x %>%
          paste(collapse = " + ")
      }) %>%
      unlist() %>%
      paste(collapse = " = ")
    new_equation
  })

new_equation <-
  unlist(new_equation)

rhea_reaction$new_equation <-
  new_equation

rhea_rpair_database <-
  vector(mode = "list", length = nrow(rhea_reaction))

for (i in 1:nrow(rhea_reaction)) {
  # if (!is.na(rhea_reaction$reaction_KEGG[i])) {
  #   rhea_rpair_database[[i]] <- NULL
  #   next()
  # }

  x <-
    rhea_reaction$new_equation[i] %>%
    stringr::str_split("\\=") %>%
    `[[`(1) %>%
    stringr::str_trim() %>%
    purrr::map(function(y) {
      y <-
        y %>%
        stringr::str_split(" \\+ ") %>%
        `[[`(1) %>%
        stringr::str_trim()
      chebi_compound_data[match(y, chebi_compound_data$compound_CHEBI_ID), ] %>%
        dplyr::filter(!is.na(compound_KEGG_ID) |
                        !is.na(compound_HMDB_ID))
        # dplyr::filter(compound_mz > 18.01056)
    })


  if (nrow(x[[1]]) == 0 | nrow(x[[2]]) == 0) {
    rhea_rpair_database[[i]] <- NULL
    next()
  }

  if((i%in%10) %in% seq(1, 100000, 1)){
    cat(i, " ")
  }

  rhea_rpair_database[[i]] <-
    purrr::map(seq_len(nrow(x[[1]])), function(index1) {
      purrr::map(seq_len(nrow(x[[2]])), function(index2) {
        cbind(
          x[[1]][index1,] %>%
            dplyr::rename(
              from_compound_CHEBI_ID = compound_CHEBI_ID,
              from_compound_name = compound_name,
              from_compound_mz = compound_mz,
              from_compound_HMDB_ID = compound_HMDB_ID,
              from_compound_KEGG_ID = compound_KEGG_ID
            ),
          x[[2]][index2,] %>%
            dplyr::rename(
              to_compound_CHEBI_ID = compound_CHEBI_ID,
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


  rhea_rpair_database[[i]] <-
    data.frame(
      rhea_rpair_database[[i]],
      rhea_reaction[i, ] %>%
        dplyr::select(
          reaction_Enzyme_EC_number = EC_number,
          reaction_ECOCYC_ID,
          reaction_METACYC_ID,
          reaction_KEGG_ID = reaction_KEGG,
          reaction_REACTOME_ID = reaction_REACTOME
        )
    )
}


rhea_rpair_database <-
  rhea_rpair_database %>%
  dplyr::bind_rows()

dim(rhea_rpair_database)

#####
library(plyr)

temp <-
  rhea_rpair_database %>%
  plyr::dlply(.variables = .(from_compound_name, to_compound_name))


temp <-
  temp %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }

    x$reaction_Enzyme_EC_number <-
      x$reaction_Enzyme_EC_number[!is.na(x$reaction_Enzyme_EC_number)] %>%
      paste(collapse = "{}")

    x$reaction_ECOCYC_ID <-
      x$reaction_ECOCYC_ID[!is.na(x$reaction_ECOCYC_ID)] %>%
      paste(collapse = "{}")

    x$reaction_METACYC_ID <-
      x$reaction_METACYC_ID[!is.na(x$reaction_METACYC_ID)] %>%
      paste(collapse = "{}")

    x$reaction_KEGG_ID <-
      x$reaction_KEGG_ID[!is.na(x$reaction_KEGG_ID)] %>%
      paste(collapse = "{}")

    x$reaction_REACTOME_ID <-
      x$reaction_REACTOME_ID[!is.na(x$reaction_REACTOME_ID)] %>%
      paste(collapse = "{}")

    x %>%
      dplyr::distinct(from_compound_name,
                      to_compound_name,
                      .keep_all = TRUE)

  }) %>%
  dplyr::bind_rows()

temp <-
  temp %>%
  dplyr::filter(from_compound_name != to_compound_name)

rpair_id <-
  seq_len(nrow(temp)) %>%
  purrr::map(function(i) {
    c(temp$from_compound_name[i],
      temp$to_compound_name[i]) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

temp <-
  temp %>%
  dplyr::mutate(rpair_id) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE) %>%
  dplyr::select(-rpair_id)

rhea_rpair_database <-
  temp %>%
  tibble::as_tibble()

save(rhea_rpair_database, file = "rhea_rpair_database")
