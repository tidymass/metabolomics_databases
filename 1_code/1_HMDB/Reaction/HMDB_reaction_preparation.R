setwd(masstools::get_project_wd())
rm(list = ls())

setwd("2_data/HMDB/Reaction/")

library(tidyverse)
library(xml2)
library(stringr)
library(massdatabase)

# hmdb_reaction_database <- vector(mode = "list", length = 18203)
#
# for (i in 1:18203) {
#   cat(i, " ")
#   hmdb_reaction_database[[i]] <-
#     request_hmdb_reaction(reaction_id = as.character(i))
# }
#
#
# save(hmdb_reaction_database, file = "hmdb_reaction_database")
load("hmdb_reaction_database")

# hmdb_reaction_database <-
#   hmdb_reaction_database %>%
#   dplyr::bind_rows()
#
# hmdb_reaction_universal_database <-
#   hmdb_reaction_database
#
# save(hmdb_reaction_universal_database, file = "hmdb_reaction_universal_database")

load("hmdb_reaction_universal_database")
