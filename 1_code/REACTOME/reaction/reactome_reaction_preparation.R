no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
source("R/read_sbml_data.R")
source("R/19_REACTOME.R")
setwd("2_data/REACTOME/reaction/")

library(metid)
library(tidyverse)
library(XML)
library(SBMLR)

reactome_reaction_human <-
  readr::read_delim("reactome_reaction_exporter_v81.txt") %>%
  dplyr::select(reaction_id, reaction_name) %>%
  dplyr::rename(reaction_REACTOME_ID = reaction_id,
                reaction_REACTOME_name = reaction_name) %>%
  dplyr::distinct(.keep_all = TRUE)

reactome_copund_data <-
  readr::read_delim(
    "https://reactome.org/download/current/ChEBI2Reactome_PE_Pathway.txt",
    col_names = FALSE
  ) %>%
  dplyr::select(-c(X4, X5, X6, X7)) %>%
  dplyr::rename(
    compound_CHEBI_ID = X1,
    compound_REACTOME_ID = X2,
    compound_REACTOME_name = X3,
    species = X8
  ) %>%
  dplyr::mutate(compound_CHEBI_ID =
                  paste0("CHEBI:", compound_CHEBI_ID))

reactome_copund_data <-
  reactome_copund_data %>%
  dplyr::filter(species == "Homo sapiens") %>%
  dplyr::select(-species) %>%
  dplyr::distinct()


###download the reaction data
reactome_reaction_human_database <-
  vector(mode = "list", length = nrow(reactome_reaction_human))

for (i in 1:nrow(reactome_reaction_human)) {
  if((i / 10) %in% seq(1, 100000, 1)){
    cat(i, " ")
  }

  reactome_reaction_human_database[[i]] <-
    request_reactome_reaction(reaction_id = reactome_reaction_human$reaction_REACTOME_ID[i])
}


####Add CHEBI ID to compounds
for (i in 11:nrow(reactome_reaction_human)) {
  cat(i, " ")
  if (!is.null(reactome_reaction_human_database[[i]])) {
    reactome_reaction_human_database[[i]]$reactants <-
      reactome_reaction_human_database[[i]]$reactants %>%
      dplyr::left_join(reactome_copund_data,
                       by = c("name" = "compound_REACTOME_name"))
  }

  if (!is.null(reactome_reaction_human_database[[i]]$products)) {
    reactome_reaction_human_database[[i]]$products <-
      reactome_reaction_human_database[[i]]$products %>%
      dplyr::left_join(reactome_copund_data,
                       by = c("name" = "compound_REACTOME_name"))
  }
}


save(reactome_reaction_human_database,
     file = "reactome_reaction_human_database")
