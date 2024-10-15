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

modelseed_reaction <-
  readr::read_tsv("reactions.tsv")

modelseed_reaction <-
  modelseed_reaction %>%
  dplyr::rename(reaction_MODELSEED_ID = id,
                reaction_MODELSEED_name = name)

other_id <-
  seq_len(nrow(modelseed_reaction)) %>%
  purrr::map(function(i) {
    if((i / 10) %in% seq(1,100000,1)){
      cat(i, " ")
    }

    x <- modelseed_reaction$aliases[i] %>%
      stringr::str_split(pattern = "\\|") %>%
      `[[`(1)

    reaction_ARACYC_ID <-
      grep("^AraCyc", x, value = TRUE) %>%
      stringr::str_replace("AraCyc: ", "") %>%
      stringr::str_trim()

    reaction_BIGG_ID <-
      grep("^BiGG", x, value = TRUE) %>%
      stringr::str_replace("BiGG: ", "") %>%
      stringr::str_trim()

    reaction_CHLAMYCYC_ID <-
      grep("^ChlamyCyc", x, value = TRUE) %>%
      stringr::str_replace("ChlamyCyc: ", "") %>%
      stringr::str_trim()

    reaction_BRACHYCYC_ID <-
      grep("^BrachyCyc", x, value = TRUE) %>%
      stringr::str_replace("BrachyCyc: ", "") %>%
      stringr::str_trim()

    reaction_KEGG_ID <-
      grep("^KEGG", x, value = TRUE) %>%
      stringr::str_replace("KEGG: ", "") %>%
      stringr::str_trim()

    reaction_METACYC_ID <-
      grep("^MetaCyc", x, value = TRUE) %>%
      stringr::str_replace("MetaCyc: ", "") %>%
      stringr::str_trim()

    reaction_ECOCYC_ID <-
      grep("^EcoCyc", x, value = TRUE) %>%
      stringr::str_replace("EcoCyc: ", "") %>%
      stringr::str_trim()

    reaction_RICECYC_ID <-
      grep("^RiceCyc", x, value = TRUE) %>%
      stringr::str_replace("RiceCyc: ", "") %>%
      stringr::str_trim()

    reaction_CORNCYC_ID <-
      grep("^CornCyc", x, value = TRUE) %>%
      stringr::str_replace("CornCyc: ", "") %>%
      stringr::str_trim()

    reaction_MAIZECYC_ID <-
      grep("^MaizeCyc", x, value = TRUE) %>%
      stringr::str_replace("MaizeCyc: ", "") %>%
      stringr::str_trim()

    reaction_POPLARCYC_ID <-
      grep("^PoplarCyc", x, value = TRUE) %>%
      stringr::str_replace("PoplarCyc: ", "") %>%
      stringr::str_trim()

    reaction_PLANTCYC_ID <-
      grep("^PlantCyc", x, value = TRUE) %>%
      stringr::str_replace("PlantCyc: ", "") %>%
      stringr::str_trim()

    reaction_SOYCYC_ID <-
      grep("^SoyCyc", x, value = TRUE) %>%
      stringr::str_replace("SoyCyc: ", "") %>%
      stringr::str_trim()

    reaction_SORGHUMCYC_ID <-
      grep("^SorghumCyc", x, value = TRUE) %>%
      stringr::str_replace("SorghumCyc: ", "") %>%
      stringr::str_trim()

    reaction_SORGHUMCYC_ID <-
      grep("^SorghumCyc", x, value = TRUE) %>%
      stringr::str_replace("SorghumCyc: ", "") %>%
      stringr::str_trim()

    data.frame(
      reaction_ARACYC_ID = ifelse(length(reaction_ARACYC_ID) == 0, NA, reaction_ARACYC_ID),
      reaction_BIGG_ID = ifelse(length(reaction_BIGG_ID) == 0, NA, reaction_BIGG_ID),
      reaction_CHLAMYCYC_ID = ifelse(
        length(reaction_CHLAMYCYC_ID) == 0,
        NA,
        reaction_CHLAMYCYC_ID
      ),
      reaction_BRACHYCYC_ID = ifelse(
        length(reaction_BRACHYCYC_ID) == 0,
        NA,
        reaction_BRACHYCYC_ID
      ),
      reaction_KEGG_ID = ifelse(length(reaction_KEGG_ID) == 0, NA, reaction_KEGG_ID),
      reaction_METACYC_ID = ifelse(length(reaction_METACYC_ID) == 0, NA, reaction_METACYC_ID),
      reaction_ECOCYC_ID = ifelse(length(reaction_ECOCYC_ID) == 0, NA, reaction_ECOCYC_ID),
      reaction_RICECYC_ID = ifelse(length(reaction_RICECYC_ID) == 0, NA, reaction_RICECYC_ID),
      reaction_CORNCYC_ID = ifelse(length(reaction_CORNCYC_ID) == 0, NA, reaction_CORNCYC_ID),
      reaction_MAIZECYC_ID = ifelse(length(reaction_MAIZECYC_ID) == 0, NA, reaction_MAIZECYC_ID),
      reaction_POPLARCYC_ID = ifelse(
        length(reaction_POPLARCYC_ID) == 0,
        NA,
        reaction_POPLARCYC_ID
      ),
      reaction_PLANTCYC_ID = ifelse(length(reaction_PLANTCYC_ID) == 0, NA, reaction_PLANTCYC_ID),
      reaction_SOYCYC_ID = ifelse(length(reaction_SOYCYC_ID) == 0, NA, reaction_SOYCYC_ID),
      reaction_SORGHUMCYC_ID = ifelse(
        length(reaction_SORGHUMCYC_ID) == 0,
        NA,
        reaction_SORGHUMCYC_ID
      )
    )
  })

length(other_id)

dim(modelseed_reaction)

other_id <-
  other_id %>%
  dplyr::bind_rows() %>%
  as.data.frame()

modelseed_reaction <-
  cbind(modelseed_reaction, other_id)

save(modelseed_reaction, file = "modelseed_reaction")
