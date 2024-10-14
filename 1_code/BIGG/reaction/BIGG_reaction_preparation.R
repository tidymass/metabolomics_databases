## http://bigg.ucsd.edu/
##
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())
# source("R/BIGG.R")
setwd("other_files/BIGG/")

load("bigg_ms1.rda")

library(jsonlite)
library(rjson)
library(massdatabase)
library(tidyverse)

####download universal reactions
# model_info <-
#   request_bigg_model_info()
#
# unique(model_info$organism)
#
# ######human information
# bigg_id <-
#   model_info %>%
#   dplyr::filter(organism == "Homo sapiens")
#
# reaction_info <-
#   purrr::map(bigg_id$bigg_id, function(x) {
#     request_bigg_reaction_info(model_bigg_id = x)
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame() %>%
#   dplyr::distinct(bigg_id, .keep_all = TRUE)
#
# bigg_reaction_human <-
#   vector(mode = "list", length = nrow(reaction_info))
#
# for (i in 1:length(bigg_reaction_human)) {
#   cat(i, " ")
#   # Sys.sleep(time = "5")
#   bigg_reaction_human[[i]] <-
#     request_bigg_universal_reaction(reaction_id = reaction_info$bigg_id[i],
#                                     return_form = "data.frame")
# }
#
# save(bigg_reaction_human, file = "bigg_reaction_human")

load("bigg_reaction_human")

lapply(bigg_reaction_human, function(x) {
  ncol(x)
}) %>%
  unlist()

colnames_unique <-
  lapply(bigg_reaction_human, function(x) {
    colnames(x)
  }) %>%
  unlist() %>%
  unique()

colnames_unique

bigg_reaction_human <-
  bigg_reaction_human %>%
  purrr::map(function(x) {
    diff_name <- setdiff(colnames_unique, colnames(x))
    if (length(diff_name) > 0) {
      new_info <-
        matrix(NA, nrow = nrow(x), ncol = length(diff_name)) %>%
        as.data.frame()
      colnames(new_info) <- diff_name
      x <-
        cbind(x, new_info)
    }
    x[, colnames_unique]
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

head(bigg_reaction_human)

bigg_reaction_human <-
  bigg_reaction_human %>%
  dplyr::select(-`}, bigg_id:`)

which(is.na(bigg_reaction_human$reaction_string))

bigg_reaction_human <-
  bigg_reaction_human %>%
  dplyr::filter(!is.na(reaction_string))

###remove the reaction that are from KEGG

# bigg_reaction_human <-
#   bigg_reaction_human %>%
#   dplyr::filter(is.na(KEGG))

bigg_reaction_human$reaction_string[1]
bigg_reaction_human$reaction_string[2]

from_to_bigg_id  <-
  bigg_reaction_human$reaction_string %>%
  purrr::map(function(x) {
    x <-
      stringr::str_split(x, "\\&\\#8652\\;")[[1]] %>%
      stringr::str_trim()

    x <-
      x %>%
      purrr::map(function(y) {
        id <-
          stringr::str_split(y, "\\+")[[1]] %>%
          stringr::str_trim() %>%
          stringr::str_replace("\\_c$", "") %>%
          stringr::str_replace("\\_m$", "") %>%
          stringr::str_replace("\\_x$", "") %>%
          stringr::str_replace("\\_e$", "") %>%
          stringr::str_replace("\\_r$", "") %>%
          stringr::str_replace("\\_l$", "") %>%
          stringr::str_replace("\\_g$", "") %>%
          stringr::str_replace("\\_n$", "") %>%
          stringr::str_replace("\\_i$", "") %>%
          stringr::str_replace("^[0-9]{1,4}.[0-9]{1,4} ", "") %>%
          stringr::str_trim() %>%
          paste0(collapse = "{}")
      }) %>%
      unlist()

    data.frame(from = x[1],
               to = x[2])
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

bigg_reaction_database_human <-
  data.frame(from_to_bigg_id,
             bigg_reaction_human[, c(
               "bigg_id",
               "name",
               "BioCyc",
               "KEGG",
               "MetaNetX",
               "Reactome",
               "RHEA",
               "SEED",
               "EC_Number"
             )]) %>%
  dplyr::rename(
    reaction_BIGG_ID = bigg_id,
    reaction_BIGG_name = name,
    reaction_BIOCYC_ID = BioCyc,
    reaction_KEGG_ID = KEGG,
    reaction_METANETX_ID = MetaNetX,
    reaction_REACTOME_ID = Reactome,
    reaction_RHEA_ID = RHEA,
    reaction_MODELSEED_ID = SEED,
    reaction_Enzyme_EC_Number = EC_Number
  )

na_idx <-
  seq_len(nrow(bigg_reaction_database_human)) %>%
  purrr::map(function(i) {
    from_id <-
      stringr::str_split(bigg_reaction_database_human$from[i], "\\{\\}")[[1]]
    to_id <-
      stringr::str_split(bigg_reaction_database_human$to[i], "\\{\\}")[[1]]

    from_idx <-
      match(from_id, bigg_ms1@spectra.info$Lab.ID)
    to_idx <-
      match(to_id, bigg_ms1@spectra.info$Lab.ID)

    if (any(is.na(from_idx)) | any(is.na(to_idx))) {
      return(i)
    } else{
      return(NULL)
    }

  }) %>%
  unlist()

na_idx

na_id <-
  seq_len(nrow(bigg_reaction_database_human)) %>%
  purrr::map(function(i) {
    from_id <-
      stringr::str_split(bigg_reaction_database_human$from[i], "\\{\\}")[[1]]
    to_id <-
      stringr::str_split(bigg_reaction_database_human$to[i], "\\{\\}")[[1]]

    from_id <-
      from_id[!from_id %in% bigg_ms1@spectra.info$Lab.ID]

    to_id <-
      to_id[!to_id %in% bigg_ms1@spectra.info$Lab.ID]

    if (length(from_id) > 0 | length(to_id) > 0) {
      return(unique(c(from_id, to_id)))
    } else{
      return(NULL)
    }
  }) %>%
  unlist() %>%
  unique()

all_metabolite_id <-
  unique(c(unique(unlist(
    stringr::str_split(bigg_reaction_database_human$from, "\\{\\}")
  )),
  unique(unlist(
    stringr::str_split(bigg_reaction_database_human$to, "\\{\\}")
  ))))

idx <-
  match(all_metabolite_id, bigg_ms1@spectra.info$Lab.ID)

idx <- idx[!is.na(idx)]

# bigg_ms1@spectra.info[idx, c("Lab.ID",
#                              "Compound.name",
#                              "HMDB.ID",
#                              "KEGG.ID",
#                              "Formula",
#                              "mz")]

save(bigg_reaction_database_human, file = "bigg_reaction_database_human")
