no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("2_data/RECON3/reaction/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
library(R.matlab)

rm(list = ls())

load("recon3_reaction_data")
load("../compound/recon3_compound_data")

recon3_compound_data <-
  recon3_compound_data %>%
  dplyr::select(
    compound_RECON3_ID = abbreviation,
    compound_name = Compound.name,
    compound_KEGG_ID = KEGG_ID,
    compound_PUBCHEM_ID = PUBCHEM_ID,
    compound_CHEBI_ID = CHEBl_ID,
    compound_HMDB_ID = HMDB_ID,
    compound_FOODB_ID = FOODB_ID,
    compound_CHEMSPIDER_ID = CHEMSPIDER_ID,
    compound_BIGG_ID = BIGG_ID,
    compound_DRUGBANK_ID = DRUGBANK_ID,
    compound_METLIN_ID = METLIN_ID,
    compound_CAS_ID = CAS_ID,
    compound_INCHIKEY = INCHIKEY,
    compound_INCHISTRING = inchiString,
    compound_SMILE = SMILE
  ) %>%
  dplyr::filter(!is.na(compound_HMDB_ID) | !is.na(compound_KEGG_ID))

recon3_rpair_database <-
  vector(mode = "list", length = nrow(recon3_reaction_data))

for (i in 1:nrow(recon3_reaction_data)) {
  # cat(i, " ")
  # if (!is.na(recon3_reaction_data$reaction_KEGG_ID[i])) {
  #   recon3_rpair_database[[i]] <- NULL
  #   next()
  # }

  x <-
    recon3_reaction_data[i,]$formula

  x <-
    x %>%
    stringr::str_split("(\\-\\>)|(\\<\\=\\>)") %>%
    `[[`(1) %>%
    stringr::str_trim() %>%
    purrr::map(function(y) {
      y <-
        stringr::str_replace_all(y, "\\[[a-z]{1,2}\\]", "") %>%
        stringr::str_split("\\+") %>%
        `[[`(1) %>%
        stringr::str_trim()
      recon3_compound_data[match(y, recon3_compound_data$compound_RECON3_ID),] %>%
        dplyr::filter(!is.na(compound_RECON3_ID)) %>%
        dplyr::filter(!compound_RECON3_ID %in% c("h", "o2", "h2o"))
    })

  if (nrow(x[[1]]) == 0 | nrow(x[[2]]) == 0) {
    recon3_rpair_database[[i]] <- NULL
    next()
  }

  if ((i / 20) %in% seq(1, 100000, 1)) {
    cat(i, " ")
  }

  recon3_rpair_database[[i]] <-
    purrr::map(seq_len(nrow(x[[1]])), function(index1) {
      purrr::map(seq_len(nrow(x[[2]])), function(index2) {
        cbind(
          x[[1]][index1,] %>%
            dplyr::rename(
              from_compound_RECON3_ID = compound_RECON3_ID,
              from_compound_name = compound_name,
              from_compound_KEGG_ID = compound_KEGG_ID,
              from_compound_PUBCHEM_ID = compound_PUBCHEM_ID,
              from_compound_CHEBI_ID = compound_CHEBI_ID,
              from_compound_HMDB_ID = compound_HMDB_ID,
              from_compound_FOODB_ID = compound_FOODB_ID,
              from_compound_CHEMSPIDER_ID = compound_CHEMSPIDER_ID,
              from_compound_BIGG_ID = compound_BIGG_ID,
              from_compound_DRUGBANK_ID = compound_DRUGBANK_ID,
              from_compound_METLIN_ID = compound_METLIN_ID,
              from_compound_CAS_ID = compound_CAS_ID,
              from_compound_INCHIKEY = compound_INCHIKEY,
              from_compound_INCHISTRING = compound_INCHISTRING,
              from_compound_SMILE = compound_SMILE
            ),
          x[[2]][index2,] %>%
            dplyr::rename(
              to_compound_RECON3_ID = compound_RECON3_ID,
              to_compound_name = compound_name,
              to_compound_KEGG_ID = compound_KEGG_ID,
              to_compound_PUBCHEM_ID = compound_PUBCHEM_ID,
              to_compound_CHEBI_ID = compound_CHEBI_ID,
              to_compound_HMDB_ID = compound_HMDB_ID,
              to_compound_FOODB_ID = compound_FOODB_ID,
              to_compound_CHEMSPIDER_ID = compound_CHEMSPIDER_ID,
              to_compound_BIGG_ID = compound_BIGG_ID,
              to_compound_DRUGBANK_ID = compound_DRUGBANK_ID,
              to_compound_METLIN_ID = compound_METLIN_ID,
              to_compound_CAS_ID = compound_CAS_ID,
              to_compound_INCHIKEY = compound_INCHIKEY,
              to_compound_INCHISTRING = compound_INCHISTRING,
              to_compound_SMILE = compound_SMILE
            )
        )
      }) %>%
        dplyr::bind_rows()

    }) %>%
    dplyr::bind_rows()

  recon3_rpair_database[[i]] <-
    cbind(recon3_rpair_database[[i]],
          recon3_reaction_data[i, c(
            "abbreviation",
            "formula",
            "ecnumber",
            "reaction_KEGG_ID",
            "reaction_MODELSEED_ID",
            "reaction_METANETX_ID"
          )]) %>%
    dplyr::rename(
      reaction_RECON3_ID = abbreviation,
      reaction_RECON3_equation = formula,
      reaction_Enzyme_EC_number = ecnumber
    )

}


recon3_rpair_database <-
  recon3_rpair_database %>%
  dplyr::bind_rows()

dim(recon3_rpair_database)

#####
library(plyr)

temp <-
  recon3_rpair_database %>%
  plyr::dlply(.variables = .(from_compound_RECON3_ID, to_compound_RECON3_ID))

temp <-
temp %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }

    x$reaction_Enzyme_EC_number <-
      x$reaction_Enzyme_EC_number[!is.na(x$reaction_Enzyme_EC_number)] %>%
      paste(collapse = "{}")

    x$reaction_RECON3_ID <-
      x$reaction_RECON3_ID[!is.na(x$reaction_RECON3_ID)] %>%
      paste(collapse = "{}")

    x$reaction_RECON3_equation <-
      x$reaction_RECON3_equation[!is.na(x$reaction_RECON3_equation)] %>%
      paste(collapse = "{}")

    x$reaction_KEGG_ID <-
      x$reaction_KEGG_ID[!is.na(x$reaction_KEGG_ID)] %>%
      paste(collapse = "{}")

    x$reaction_MODELSEED_ID <-
      x$reaction_MODELSEED_ID[!is.na(x$reaction_MODELSEED_ID)] %>%
      paste(collapse = "{}")

    x$reaction_METANETX_ID <-
      x$reaction_METANETX_ID[!is.na(x$reaction_METANETX_ID)] %>%
      paste(collapse = "{}")

    x %>%
      dplyr::distinct(from_compound_RECON3_ID,
                      to_compound_RECON3_ID,
                      .keep_all = TRUE)

  }) %>%
  dplyr::bind_rows()


temp[temp == ''] <- NA


temp <-
  temp %>%
  dplyr::filter(from_compound_RECON3_ID != to_compound_RECON3_ID)

rpair_id <-
  seq_len(nrow(temp)) %>%
  purrr::map(function(i){
    c(temp$from_compound_RECON3_ID[i],
      temp$to_compound_RECON3_ID[i]) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

temp <-
temp %>%
  dplyr::mutate(rpair_id) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE) %>%
  dplyr::select(-rpair_id)

recon3_rpair_database <-
  temp %>%
  tibble::as_tibble()

dim(recon3_rpair_database)

save(recon3_rpair_database, file = "recon3_rpair_database")
