no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("other_files/RHEA/reaction/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)

rm(list = ls())

rhea_reaction <-
  readr::read_tsv("Rhea.tsv")

rhea_reaction <-
  rhea_reaction %>%
  dplyr::rename(
    reaction_RHEA_ID = `Reaction identifier`,
    CHEBI_ID = `ChEBI identifier`,
    CHEBI_name = `ChEBI name`,
    Enzyme_class = `Enzyme class`,
    EC_number = `EC number`,
    reaction_ECOCYC_ID = `Cross-reference (EcoCyc)`,
    reaction_METACYC_ID = `Cross-reference (MetaCyc)`,
    reaction_KEGG = `Cross-reference (KEGG)`,
    reaction_REACTOME = `Cross-reference (Reactome)`
  ) %>%
  dplyr::mutate(reaction_RHEA_ID = stringr::str_replace_all(reaction_RHEA_ID, "RHEA:", ""))

rhea_reaction$CHEBI_name <-
  rhea_reaction$CHEBI_name %>%
  stringr::str_replace_all(";", "{}")

rhea_reaction$CHEBI_ID <-
  rhea_reaction$CHEBI_ID %>%
  stringr::str_replace_all(";", "{}")

rhea_reaction$EC_number <-
  rhea_reaction$EC_number %>%
  stringr::str_replace_all('EC:', "")

rhea_reaction$reaction_ECOCYC_ID <-
  rhea_reaction$reaction_ECOCYC_ID %>%
  stringr::str_replace_all('EcoCyc:', "")

rhea_reaction$reaction_METACYC_ID <-
  rhea_reaction$reaction_METACYC_ID %>%
  stringr::str_replace_all('MetaCyc:', "")

rhea_reaction$reaction_KEGG <-
  rhea_reaction$reaction_KEGG %>%
  stringr::str_replace_all('KEGG:', "")

rhea_reaction$reaction_KEGG <-
  rhea_reaction$reaction_KEGG %>%
  stringr::str_replace_all('KEGG:', "")

rhea_reaction$reaction_REACTOME <-
  rhea_reaction$reaction_REACTOME %>%
  stringr::str_replace_all('Reactome:', "")

save(rhea_reaction, file = "rhea_reaction")
