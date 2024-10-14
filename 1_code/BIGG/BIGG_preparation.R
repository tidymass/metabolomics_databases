setwd(masstools::get_project_wd())
rm(list = ls())
source("R/3_utils.R")
setwd("other_files/WIKIPATHWAYS/")

library(tidyverse)
library(xml2)
library(stringr)
library(rWikiPathways)
listOrganisms()
hs.pathways <- listPathways('Homo sapiens')
hs.pathways


file_name <-
  dir(path = "wikipathways-20220410-gpml-Homo_sapiens", full.names = TRUE)

wikipathway_reaction_database <-
  seq_along(file_name) %>%
  purrr::map(function(i) {
    cat(i, " ")
    read_gpml(file = file_name[i])$reaction
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

wikipathway_reaction_database <-
  wikipathway_reaction_database %>%
  dplyr::filter(
    !Database %in% c(
      "",
      "IntAct",
      "KEGG Pathway",
      "WikiPathways",
      "XMetDB",
      "Database",
      "Wikidata",
      "SPIKE",
      "ChEBI",
      "KEGG Compound",
      "PATO:0002220",
      "Uniprot-TrEMBL"
    )
  ) %>%
  dplyr::distinct()

unique(wikipathway_reaction_database$Database)

save(wikipathway_reaction_database, file = "wikipathway_reaction_database")

