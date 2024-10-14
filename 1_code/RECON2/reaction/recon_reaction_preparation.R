no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("other_files/RECON2/reaction/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)

rm(list = ls())

result <-
  XML::xmlParse(file = "MODEL1603150001_url.xml")

result <-
  result %>%
  XML::xmlToList()

recon2_compound_data <-
  seq_along(result$model$listOfSpecies) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <-
      result$model$listOfSpecies[[i]]
    c(FORMULA = ifelse(length(x$notes$body) < 2, NA, x$notes$body[[2]]),
      x$.attrs)[c("FORMULA", "id", "name", "metaid", "sboTerm")]
  })

recon2_compound_data %>%
  lapply(length) %>%
  unlist()

recon2_compound_data <-
  recon2_compound_data %>%
  dplyr::bind_rows() %>%
  as.data.frame()


recon2_reaction_data <-
  seq_along(result$model$listOfReactions) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <-
      result$model$listOfReactions[[i]]


  })
