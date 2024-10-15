no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("2_data/MODELSEED/compound/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)

rm(list = ls())

modelseed_compound <-
  readr::read_tsv("compounds.tsv")

modelseed_compound <-
  modelseed_compound %>%
  dplyr::rename(compound_MODELSEED_ID = id)

other_id <-
  seq_len(nrow(modelseed_compound)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <- modelseed_compound$aliases[i] %>%
      stringr::str_split(pattern = "\\|") %>%
      `[[`(1)

    KEGG_ID <-
      grep("^KEGG", x, value = TRUE) %>%
      stringr::str_replace("KEGG: ", "") %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("; ", "{}") %>%
      stringr::str_replace_all(";", "{}")

    BIGG_ID <-
      grep("^BiGG", x, value = TRUE) %>%
      stringr::str_replace("BiGG: ", "") %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("; ", "{}") %>%
      stringr::str_replace_all(";", "{}")

    Synonym <-
      grep("^Name", x, value = TRUE) %>%
      stringr::str_replace("Name: ", "") %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("; ", "{}") %>%
      stringr::str_replace_all(";", "{}")

    data.frame(
      KEGG_ID = ifelse(length(KEGG_ID) == 0, NA, KEGG_ID),
      BIGG_ID = ifelse(length(BIGG_ID) == 0, NA, BIGG_ID),
      Synonym = ifelse(length(Synonym) == 0, NA, Synonym)
    )
  })

# all_id <-
# modelseed_compound$aliases %>%
#   purrr::map(function(x){
#     stringr::str_split(x, "\\|")[[1]] %>%
#       stringr::str_split(":") %>%
#       lapply(function(y){
#         y[1]
#       }) %>%
#       unlist()
#   }) %>%
#   unlist() %>%
#   unique()


length(other_id)

dim(modelseed_compound)

other_id <-
  other_id %>%
  dplyr::bind_rows() %>%
  as.data.frame()

modelseed_compound <-
  cbind(modelseed_compound, other_id)

save(modelseed_compound, file = "modelseed_compound")
