no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("2_data/MODELSEED/reaction/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)

rm(list = ls())


load("../compound/modelseed_compound")
load("modelseed_reaction")

modelseed_compound <-
  modelseed_compound %>%
  dplyr::filter(!is.na(inchikey) | !is.na(smiles) | is.na(KEGG_ID))

modelseed_compound$KEGG_ID <-
  modelseed_compound$KEGG_ID %>%
  stringr::str_split("\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

modelseed_rpair_database <-
  vector(mode = "list", length = nrow(modelseed_reaction))

for (i in 1:nrow(modelseed_compound)) {
  # if (!is.na(modelseed_reaction$reaction_KEGG_ID[i])) {
  #   modelseed_rpair_database[[i]] <- NULL
  #   next()
  # }

  x <-
    modelseed_reaction$code[i] %>%
    stringr::str_split(pattern = "\\<\\=\\>") %>%
    `[[`(1) %>%
    stringr::str_split("\\+") %>%
    purrr::map(function(y) {
      y <-
        y %>%
        stringr::str_trim() %>%
        stringr::str_extract_all("cpd[0-9]{4,7}") %>%
        unlist()

      modelseed_compound[match(y, modelseed_compound$compound_MODELSEED_ID),
                         c("compound_MODELSEED_ID",
                           "name",
                           "inchikey",
                           "smiles",
                           "KEGG_ID",
                           "BIGG_ID")] %>%
        dplyr::filter(!is.na(inchikey) |
                        !is.na(smiles) | !is.na(KEGG_ID))
    })

  if (nrow(x[[1]]) == 0 | nrow(x[[2]]) == 0) {
    modelseed_rpair_database[[i]] <- NULL
    next()
  }

  modelseed_rpair_database[[i]] <-
    purrr::map(seq_len(nrow(x[[1]])), function(index1) {
      purrr::map(seq_len(nrow(x[[2]])), function(index2) {
        cbind(
          x[[1]][index1,] %>%
            dplyr::rename(
              from_compound_MODELSEED_ID = compound_MODELSEED_ID,
              from_compound_name = name,
              from_compound_INCHIKEY = inchikey,
              from_compound_SMILES = smiles,
              from_compound_KEGG_ID = KEGG_ID,
              from_compound_BIGG_ID = BIGG_ID
            ),
          x[[2]][index2,] %>%
            dplyr::rename(
              to_compound_MODELSEED_ID = compound_MODELSEED_ID,
              to_compound_name = name,
              to_compound_INCHIKEY = inchikey,
              to_compound_SMILES = smiles,
              to_compound_KEGG_ID = KEGG_ID,
              to_compound_BIGG_ID = BIGG_ID
            )
        )
      }) %>%
        dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows()

  modelseed_rpair_database[[i]] <-
    data.frame(
      modelseed_rpair_database[[i]],
      modelseed_reaction[i, c(
        "reaction_MODELSEED_ID",
        "reaction_MODELSEED_name",
        "equation",
        "ec_numbers",
        "reaction_BIGG_ID",
        "reaction_KEGG_ID",
        "reaction_METACYC_ID"
      )] %>%
        dplyr::rename(
          reaction_MODELSEED_equation = equation,
          reaction_Enzyme_EC_numbers = ec_numbers
        )
    )
  if ((i / 100) %in% seq(1, 100000, 1)) {
    cat(i, " ")
  }
}

modelseed_rpair_database <-
  modelseed_rpair_database %>%
  dplyr::bind_rows()

dim(modelseed_rpair_database)


library(plyr)

####remove some compounds
rbind(
  modelseed_rpair_database[, ("from_compound_name")] %>%
    dplyr::rename(name = from_compound_name),
  modelseed_rpair_database[, ("to_compound_name")] %>%
    dplyr::rename(name = to_compound_name)
) %>%
  dplyr::distinct() %>%
  dplyr::mutate(n = nchar(name)) %>%
  dplyr::filter(n <= 3)

remove_compound <-
  c(
    "H2O",
    "CO2",
    "H2O2",
    "O2",
    "Fe+2",
    "S",
    "NH3",
    "O2-",
    "NO",
    "H2S",
    "HCN",
    "Hg",
    "I-",
    "F-",
    "Cl-",
    "CO",
    "Mg",
    "K+",
    "H2",
    "N2",
    "Ag",
    "I2",
    "Fe",
    "Na+"
  )

modelseed_rpair_database <-
  modelseed_rpair_database %>%
  dplyr::filter(!from_compound_name %in% remove_compound &
                  !to_compound_name %in% remove_compound) %>%
  tibble::as_tibble()


temp <-
  modelseed_rpair_database %>%
  plyr::dlply(.variables = .(from_compound_MODELSEED_ID, to_compound_MODELSEED_ID))

temp <-
  seq_along(temp) %>%
  purrr::map(function(i) {
    # cat(i, " ")
    x <-
      temp[[i]]
    if (nrow(x) == 1) {
      return(x)
    }

    x$reaction_MODELSEED_ID <-
      x$reaction_MODELSEED_ID[!is.na(x$reaction_MODELSEED_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_MODELSEED_name <-
      x$reaction_MODELSEED_name[!is.na(x$reaction_MODELSEED_name)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_MODELSEED_equation <-
      x$reaction_MODELSEED_equation[!is.na(x$reaction_MODELSEED_equation)] %>%
      unique() %>%
      paste(collapse = "{}")

    reaction_Enzyme_EC_numbers <-
      x$reaction_Enzyme_EC_numbers[!is.na(x$reaction_Enzyme_EC_numbers)] %>%
      unique()

    x$reaction_Enzyme_EC_numbers <-
      reaction_Enzyme_EC_numbers[reaction_Enzyme_EC_numbers != "null"] %>%
      paste(collapse = "{}")

    x$reaction_BIGG_ID <-
      x$reaction_BIGG_ID[!is.na(x$reaction_BIGG_ID)] %>%
      stringr::str_split(';') %>%
      unlist() %>%
      stringr::str_trim() %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_KEGG_ID <-
      x$reaction_KEGG_ID[!is.na(x$reaction_KEGG_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_METACYC_ID <-
      x$reaction_METACYC_ID[!is.na(x$reaction_METACYC_ID)] %>%
      unique() %>%
      paste(collapse = "{}")

    x <-
      x %>%
      dplyr::distinct(from_compound_MODELSEED_ID,
                      to_compound_MODELSEED_ID,
                      .keep_all = TRUE)

    x
  })


modelseed_rpair_database <-
  temp %>%
  dplyr::bind_rows() %>%
  tibble::as_tibble()

modelseed_rpair_database <-
  modelseed_rpair_database %>%
  dplyr::filter(from_compound_MODELSEED_ID != to_compound_MODELSEED_ID)

###remove some duplicated reaction pair
rpair_id <-
  seq_len(nrow(modelseed_rpair_database)) %>%
  purrr::map(function(i) {
    c(
      modelseed_rpair_database$from_compound_MODELSEED_ID[i],
      modelseed_rpair_database$to_compound_MODELSEED_ID[i]
    ) %>%
      sort() %>%
      paste(collapse = "{}")
  }) %>%
  unlist()

modelseed_rpair_database <-
  modelseed_rpair_database %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE) %>%
  dplyr::select(-rpair_id)

save(modelseed_rpair_database, file = "modelseed_rpair_database")
