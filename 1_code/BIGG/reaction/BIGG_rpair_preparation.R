## http://bigg.ucsd.edu/
##
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())
load("other_files/HMDB/MS1/hmdb_ms1.rda")
load("other_files/KEGG/kegg_ms1.rda")
setwd("other_files/BIGG/reaction/")

load("../bigg_ms1.rda")
load("../bigg_reaction_database_human")

bigg_rpair_database_human <-
  seq_len(nrow(bigg_reaction_database_human)) %>%
  purrr::map(function(i) {
    if ((i / 10) %in% seq(1, 100000, 1)) {
      cat(i, " ")
    }
    from_id <-
      stringr::str_split(bigg_reaction_database_human$from[i], "\\{\\}")[[1]]
    to_id <-
      stringr::str_split(bigg_reaction_database_human$to[i], "\\{\\}")[[1]]

    from_idx <-
      match(from_id, bigg_ms1@spectra.info$Lab.ID)
    to_idx <-
      match(to_id, bigg_ms1@spectra.info$Lab.ID)

    from_idx <- from_idx[!is.na(from_idx)]
    to_idx <- to_idx[!is.na(to_idx)]

    if (length(from_idx) == 0 | length(to_idx) == 0) {
      return(NULL)
    }

    from_metabolite <-
      bigg_ms1@spectra.info[from_idx, c("Lab.ID",
                                        "Compound.name",
                                        "HMDB.ID",
                                        "KEGG.ID",
                                        "Formula",
                                        "mz")] %>%
      dplyr::filter(
        !Lab.ID %in% c(
          "pi",
          "h",
          "CE5049",
          "oh1",
          "na1",
          "cyan",
          "co",
          "no",
          "fald",
          "mma",
          "o2",
          "o2s",
          "meoh",
          "HC00250",
          "h2o2",
          "nh4",
          "cl",
          "k",
          "ca2",
          "co2",
          "acald",
          "for",
          "etoh",
          "no2",
          "fe2",
          "fe3",
          "tcynt",
          "acetone",
          "M02035",
          "ac",
          "gcald",
          "urea",
          "etha",
          "CE5643",
          "zn2",
          "CE1944",
          "mthgxl",
          "CE0737"
        )
      )
    to_metabolite <-
      bigg_ms1@spectra.info[to_idx, c("Lab.ID",
                                      "Compound.name",
                                      "HMDB.ID",
                                      "KEGG.ID",
                                      "Formula",
                                      "mz")] %>%
      dplyr::filter(
        !Lab.ID %in% c(
          "pi",
          "h",
          "CE5049",
          "oh1",
          "na1",
          "cyan",
          "co",
          "no",
          "fald",
          "mma",
          "o2",
          "o2s",
          "meoh",
          "HC00250",
          "h2o2",
          "nh4",
          "cl",
          "k",
          "ca2",
          "co2",
          "acald",
          "for",
          "etoh",
          "no2",
          "fe2",
          "fe3",
          "tcynt",
          "acetone",
          "M02035",
          "ac",
          "gcald",
          "urea",
          "etha",
          "CE5643",
          "zn2",
          "CE1944",
          "mthgxl",
          "CE0737"
        )
      )

    if (nrow(from_metabolite) == 0 | nrow(to_metabolite) == 0) {
      return(NULL)
    }

    seq_len(nrow(from_metabolite)) %>%
      purrr::map(function(idx1) {
        seq_len(nrow(to_metabolite)) %>%
          purrr::map(function(idx2) {
            cbind(
              from_metabolite[idx1, , drop = FALSE] %>%
                dplyr::rename(
                  from_Lab.ID = Lab.ID,
                  from_Compound.name = Compound.name,
                  from_HMDB.ID = HMDB.ID,
                  from_KEGG.ID = KEGG.ID,
                  from_Formula = Formula,
                  from_mz = mz
                ),
              to_metabolite[idx2, , drop = FALSE] %>%
                dplyr::rename(
                  to_Lab.ID = Lab.ID,
                  to_Compound.name = Compound.name,
                  to_HMDB.ID = HMDB.ID,
                  to_KEGG.ID = KEGG.ID,
                  to_Formula = Formula,
                  to_mz = mz
                ),
              bigg_reaction_database_human[i, c(
                "reaction_BIGG_ID",
                "reaction_BIGG_name",
                "reaction_BIOCYC_ID",
                "reaction_KEGG_ID",
                "reaction_METANETX_ID",
                "reaction_REACTOME_ID",
                "reaction_RHEA_ID",
                "reaction_MODELSEED_ID",
                "reaction_Enzyme_EC_Number"
              )]
            )
          })
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

rownames(bigg_rpair_database_human) <- NULL

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    from_compound_name = from_Compound.name,
    to_compound_name = to_Compound.name,
    from_compound_HMDB_ID = from_HMDB.ID,
    to_compound_HMDB_ID = to_HMDB.ID,
    from_compound_KEGG_ID = from_KEGG.ID,
    to_compound_KEGG_ID = to_KEGG.ID,
    from_compound_formula = from_Formula,
    to_compound_formula = to_Formula,
    from_compound_mz = from_mz,
    to_compound_mz = to_mz
  )


save(bigg_rpair_database_human, file = "bigg_rpair_database_human")

####remove some compounds
remove_compound <-
  c(
    "H2O H2O",
    "Cyanate",
    "Hypochlorous acid",
    "Nitryl chloride",
    "Glycophosphatidylinositol (GPI) signal sequence (C-terminal peptide)",
    "I c",
    "Band membrane protein (universal, erythrocyte -> 2.1,3,4.1)",
    "Band membrane protein (methylated, universal, erythrocyte -> 2.1,3,4.1)",
    "Iodine c"
  )

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::filter(!from_compound_name %in% remove_compound &
                  !to_compound_name %in% remove_compound) %>%
  tibble::as_tibble()

temp <-
  rbind(
    data.frame(bigg_rpair_database_human[, c("from_compound_name",
                                             "from_compound_formula",
                                             "from_compound_mz")]) %>%
      dplyr::rename(name = from_compound_name,
                    formula = from_compound_formula,
                    mz = from_compound_mz),
    data.frame(bigg_rpair_database_human[, c("to_compound_name", "to_compound_formula", "to_compound_mz")]) %>%
      dplyr::rename(name = to_compound_name,
                    formula = to_compound_formula,
                    mz = to_compound_mz)
  )

temp <-
  temp %>%
  dplyr::filter(!is.na(name)) %>%
  dplyr::arrange(mz) %>%
  dplyr::distinct(name, .keep_all = TRUE)

temp[which(nchar(temp$formula) <= 5),]

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::rename(from_compound_BIGG_ID = from_Lab.ID,
                to_compound_BIGG_ID = to_Lab.ID)

###remove water

library(plyr)

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  plyr::dlply(.variables = .(from_compound_name, to_compound_name)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }

    reaction_BIGG_ID <-
      x$reaction_BIGG_ID[!is.na(x$reaction_BIGG_ID)]

    if (length(reaction_BIGG_ID) == 0) {
      reaction_BIGG_ID <- NA
    } else{
      reaction_BIGG_ID <- paste(reaction_BIGG_ID, collapse = "{}")
    }

    x$reaction_BIGG_ID <- reaction_BIGG_ID

    reaction_BIGG_name <-
      x$reaction_BIGG_name[!is.na(x$reaction_BIGG_name)]

    if (length(reaction_BIGG_name) == 0) {
      reaction_BIGG_name <- NA
    } else{
      reaction_BIGG_name <- paste(reaction_BIGG_name, collapse = "{}")
    }
    x$reaction_BIGG_name <- reaction_BIGG_name


    reaction_BIOCYC_ID <-
      x$reaction_BIOCYC_ID[!is.na(x$reaction_BIOCYC_ID)]

    if (length(reaction_BIOCYC_ID) == 0) {
      reaction_BIOCYC_ID <- NA
    } else{
      reaction_BIOCYC_ID <- paste(reaction_BIOCYC_ID, collapse = "{}")
    }

    x$reaction_BIOCYC_ID <- reaction_BIOCYC_ID

    reaction_KEGG_ID <-
      x$reaction_KEGG_ID[!is.na(x$reaction_KEGG_ID)]

    if (length(reaction_KEGG_ID) == 0) {
      reaction_KEGG_ID <- NA
    } else{
      reaction_KEGG_ID <- paste(reaction_KEGG_ID, collapse = "{}")
    }
    x$reaction_KEGG_ID <- reaction_KEGG_ID

    reaction_METANETX_ID <-
      x$reaction_METANETX_ID[!is.na(x$reaction_METANETX_ID)]

    if (length(reaction_METANETX_ID) == 0) {
      reaction_METANETX_ID <- NA
    } else{
      reaction_METANETX_ID <- paste(reaction_METANETX_ID, collapse = "{}")
    }
    x$reaction_METANETX_ID <- reaction_METANETX_ID

    reaction_REACTOME_ID <-
      x$reaction_REACTOME_ID[!is.na(x$reaction_REACTOME_ID)]

    if (length(reaction_REACTOME_ID) == 0) {
      reaction_REACTOME_ID <- NA
    } else{
      reaction_REACTOME_ID <- paste(reaction_REACTOME_ID, collapse = "{}")
    }
    x$reaction_REACTOME_ID <- reaction_REACTOME_ID

    reaction_RHEA_ID <-
      x$reaction_RHEA_ID[!is.na(x$reaction_RHEA_ID)]

    if (length(reaction_RHEA_ID) == 0) {
      reaction_RHEA_ID <- NA
    } else{
      reaction_RHEA_ID <- paste(reaction_RHEA_ID, collapse = "{}")
    }
    x$reaction_RHEA_ID <- reaction_RHEA_ID

    reaction_MODELSEED_ID <-
      x$reaction_MODELSEED_ID[!is.na(x$reaction_MODELSEED_ID)]

    if (length(reaction_MODELSEED_ID) == 0) {
      reaction_MODELSEED_ID <- NA
    } else{
      reaction_MODELSEED_ID <-
        paste(reaction_MODELSEED_ID, collapse = "{}")
    }
    x$reaction_MODELSEED_ID <- reaction_MODELSEED_ID

    reaction_Enzyme_EC_Number <-
      x$reaction_Enzyme_EC_Number[!is.na(x$reaction_Enzyme_EC_Number)]

    if (length(reaction_Enzyme_EC_Number) == 0) {
      reaction_Enzyme_EC_Number <- NA
    } else{
      reaction_Enzyme_EC_Number <-
        paste(reaction_Enzyme_EC_Number, collapse = "{}")
    }
    x$reaction_Enzyme_EC_Number <- reaction_Enzyme_EC_Number

    x %>%
      dplyr::distinct(from_compound_name,
                      to_compound_name,
                      .keep_all = TRUE)

  }) %>%
  dplyr::bind_rows()


dim(bigg_rpair_database_human)

####remove the reaction paris that compound without KEGG ID and HMDB_ID

###Add new HMDB_ID and KEGG_ID

name_id <-
  rbind(
    bigg_rpair_database_human[, c("from_compound_name", "from_compound_HMDB_ID")] %>%
      dplyr::rename(name = from_compound_name,
                    HMDB_ID = from_compound_HMDB_ID) %>%
      dplyr::filter(is.na(HMDB_ID)),
    bigg_rpair_database_human[, c("to_compound_name", "to_compound_HMDB_ID")] %>%
      dplyr::rename(name = to_compound_name,
                    HMDB_ID = to_compound_HMDB_ID) %>%
      dplyr::filter(is.na(HMDB_ID))
  ) %>%
  dplyr::distinct(name, .keep_all = TRUE)

new_hmdb_id <-
  seq_len(nrow(name_id)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    idx <-
      match(name_id$name[i],
            hmdb_ms1@spectra.info$Compound.name)
    if (!is.na(idx)) {
      return(hmdb_ms1@spectra.info$HMDB.ID[idx])
    }

    idx <-
      tryCatch(
        stringr::str_detect(hmdb_ms1@spectra.info$Synonyms,
                            name_id$name[i]) %>%
          which(),
        error = function(e)
          NULL
      )

    if (length(idx) == 0) {
      return(NA)
    }

    if (length(idx) == 1) {
      return(hmdb_ms1@spectra.info$HMDB.ID[idx])
    }

    temp_idx <-
      hmdb_ms1@spectra.info$Synonyms[idx] %>%
      stringr::str_split("\\{\\}") %>%
      lapply(function(z) {
        any(z == name_id$name[i])
      }) %>%
      unlist() %>%
      which()

    if (length(temp_idx) == 0) {
      return(NA)
    }
    return(hmdb_ms1@spectra.info$HMDB.ID[idx[temp_idx[1]]])
  })

new_hmdb_id <-
  new_hmdb_id %>%
  unlist()

name_id$HMDB_ID <-
  new_hmdb_id

name_id <-
  name_id %>%
  dplyr::filter(!is.na(HMDB_ID))

from_compound_HMDB_ID <-
  bigg_rpair_database_human[, c("from_compound_name", "from_compound_HMDB_ID")] %>%
  dplyr::left_join(name_id, by = c("from_compound_name" = "name")) %>%
  apply(1, function(x) {
    x <- as.character(x[2:3])
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(tail(x, 1))
  })

bigg_rpair_database_human$from_compound_HMDB_ID <-
  from_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  lapply(function(x) {
    x[1]
  }) %>%
  unlist()

to_compound_HMDB_ID <-
  bigg_rpair_database_human[, c("to_compound_name", "to_compound_HMDB_ID")] %>%
  dplyr::left_join(name_id, by = c("to_compound_name" = "name")) %>%
  apply(1, function(x) {
    x <- as.character(x[2:3])
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(tail(x, 1))
  })

bigg_rpair_database_human$to_compound_HMDB_ID <-
  to_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  lapply(function(x) {
    x[1]
  }) %>%
  unlist()


###Add new HMDB_ID and KEGG_ID
name_id <-
  rbind(
    bigg_rpair_database_human[, c("from_compound_name", "from_compound_KEGG_ID")] %>%
      dplyr::rename(name = from_compound_name,
                    KEGG_ID = from_compound_KEGG_ID) %>%
      dplyr::filter(is.na(KEGG_ID)),
    bigg_rpair_database_human[, c("to_compound_name", "to_compound_KEGG_ID")] %>%
      dplyr::rename(name = to_compound_name,
                    KEGG_ID = to_compound_KEGG_ID) %>%
      dplyr::filter(is.na(KEGG_ID))
  ) %>%
  dplyr::distinct(name, .keep_all = TRUE)

new_kegg_id <-
  seq_len(nrow(name_id)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    idx <-
      match(name_id$name[i],
            kegg_ms1@spectra.info$Compound.name)
    if (!is.na(idx)) {
      return(kegg_ms1@spectra.info$KEGG.ID[idx])
    }

    idx <-
      tryCatch(
        stringr::str_detect(kegg_ms1@spectra.info$Synonyms,
                            name_id$name[i]) %>%
          which(),
        error = function(e)
          NULL
      )

    if (length(idx) == 0) {
      return(NA)
    }

    if (length(idx) == 1) {
      return(kegg_ms1@spectra.info$KEGG.ID[idx])
    }

    temp_idx <-
      kegg_ms1@spectra.info$Synonyms[idx] %>%
      stringr::str_split("\\{\\}") %>%
      lapply(function(z) {
        any(z == name_id$name[i])
      }) %>%
      unlist() %>%
      which()

    if (length(temp_idx) == 0) {
      return(NA)
    }
    return(kegg_ms1@spectra.info$KEGG.ID[idx[temp_idx[1]]])
  })

new_kegg_id <-
  new_kegg_id %>%
  unlist()

name_id$KEGG_ID <-
  new_kegg_id

name_id <-
  name_id %>%
  dplyr::filter(!is.na(KEGG_ID))

from_compound_KEGG_ID <-
  bigg_rpair_database_human[, c("from_compound_name", "from_compound_KEGG_ID")] %>%
  dplyr::left_join(name_id, by = c("from_compound_name" = "name")) %>%
  apply(1, function(x) {
    x <- as.character(x[2:3])
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(tail(x, 1))
  })

bigg_rpair_database_human$from_compound_KEGG_ID <-
  from_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  lapply(function(x) {
    x[1]
  }) %>%
  unlist()

to_compound_KEGG_ID <-
  bigg_rpair_database_human[, c("to_compound_name", "to_compound_KEGG_ID")] %>%
  dplyr::left_join(name_id, by = c("to_compound_name" = "name")) %>%
  apply(1, function(x) {
    x <- as.character(x[2:3])
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(tail(x, 1))
  })

bigg_rpair_database_human$to_compound_KEGG_ID <-
  to_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  lapply(function(x) {
    x[1]
  }) %>%
  unlist()


bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::filter((
    !is.na(from_compound_HMDB_ID) | !is.na(from_compound_KEGG_ID)
  ) &
    (!is.na(to_compound_HMDB_ID) |
       !is.na(to_compound_KEGG_ID)))


# check <-
# rbind(
#   bigg_rpair_database_human[,c("from_compound_name", "from_compound_HMDB_ID")] %>%
#     dplyr::rename(compound_name = from_compound_name,
#                   compound_HMDB_ID = from_compound_HMDB_ID) %>%
#     dplyr::distinct(),
#   bigg_rpair_database_human[,c("to_compound_name", "to_compound_HMDB_ID")] %>%
#     dplyr::rename(compound_name = to_compound_name,
#                   compound_HMDB_ID = to_compound_HMDB_ID) %>%
#     dplyr::distinct()
# ) %>%
#   dplyr::distinct()
#
# idx <-
# check$compound_HMDB_ID %>%
#   stringr::str_split("\\{\\}") %>%
#   lapply(length) %>%
#   unlist() %>%
#   `>`(1) %>%
#   which()
#
# check <-
# check[idx,]
#
# write.csv(check, file = "check.csv", row.names = FALSE)

save(bigg_rpair_database_human, file = "bigg_rpair_database_human")
