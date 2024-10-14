no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("other_files/KEGG")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)

rm(list = ls())

###to get KEGG database
library(KEGGgraph)
library(KEGGREST)
library(KEGGREST)
library(tidyverse)
library(massdatabase)

#####Reaction class database
# kegg_rclass_info <-
#   request_kegg_rclass_info()
#
# kegg_rpair_database <-
#   vector(mode = "list", length = nrow(kegg_rclass_info))
#
# for (i in seq_len(nrow(kegg_rclass_info))) {
#   cat(i, " ")
  # temp <-
  #   request_kegg_rclass(rclass_id = kegg_rclass_info$KEGG.ID[i])
#
#   rpair <-
#     temp$RPAIR %>%
#     stringr::str_split(" ") %>%
#     `[[`(1) %>%
#     stringr::str_trim()
#
#   rpair <- rpair[rpair != ""]
#   rpair <-
#     rpair %>%
#     purrr::map(function(x) {
#       stringr::str_split(x, pattern = "_")[[1]]
#     }) %>%
#     do.call(rbind, .) %>%
#     as.data.frame()
#
#   colnames(rpair) <- c("from_KEGG.ID", "to_KEGG.ID")
#
#   reaction_KEGG <-
#     stringr::str_split(temp$REACTION, " ") %>%
#     unlist() %>%
#     unique() %>%
#     paste0(collapse = "{}")
#   kegg_rpair_database[[i]] <-
#     data.frame(rpair, reaction_KEGG)
# }
#
# kegg_rpair_database <-
# kegg_rpair_database %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# save(kegg_rpair_database, file = "kegg_rpair_database")
load("kegg_rpair_database")


####Reaction database
kegg_reaction_info <-
  request_kegg_reaction_info()

kegg_reaction_database <-
  vector(mode = "list", length = nrow(kegg_reaction_info))

# for (i in 1:nrow(kegg_reaction_info)) {
#   cat(i, " ")
#   temp <-
#     request_kegg_reaction(kegg_reaction_info$KEGG.ID[i])
#
#   kegg_reaction_database[[i]] <-
#     data.frame(
#       KEGG.ID = ifelse(is.null(unname(temp$ENTRY)), NA, unname(temp$ENTRY)),
#       NAME = ifelse(is.null(unname(temp$NAME)), NA, unname(temp$NAME)),
#       EQUATION = ifelse(is.null(unname(temp$EQUATION)), NA, unname(temp$EQUATION)),
#       ENZYME = ifelse(is.null(unname(temp$ENZYME)), NA, unname(temp$ENZYME)),
#       DBLINKS = ifelse(is.null(unname(temp$DBLINKS)), NA, unname(temp$DBLINKS))
#     )
# }
#
# kegg_reaction_database <-
#   kegg_reaction_database %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# save(kegg_reaction_database, file = "kegg_reaction_database")
load("kegg_reaction_database")

kegg_rpair_database

# kegg_rpair_database <-
# seq_len(nrow(kegg_rpair_database)) %>%
#   purrr::map(function(i) {
#     temp_idx <-
#       match(unique(stringr::str_split(kegg_rpair_database$reaction_KEGG[i], "\\{\\}")[[1]]), kegg_reaction_database$KEGG.ID)
#     reaction_name <-
#       paste0(kegg_reaction_database$NAME[temp_idx], collapse = "{}")
#     reaction_EQUATION <-
#       paste0(kegg_reaction_database$EQUATION[temp_idx], collapse = "{}")
#
#     reaction_ENZYME <-
#       paste0(kegg_reaction_database$ENZYME[temp_idx], collapse = "{}")
#     reaction_RHEA <-
#       paste0(kegg_reaction_database$DBLINKS[temp_idx], collapse = "{}")
#
#     cbind(kegg_rpair_database[i,,drop = FALSE],
#           reaction_name,
#           reaction_EQUATION,
#           reaction_ENZYME,
#           reaction_RHEA)
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# kegg_reaction_universal_database <-
#   kegg_reaction_database
#
# save(kegg_reaction_universal_database, file = "kegg_reaction_universal_database")
load("kegg_reaction_universal_database")
