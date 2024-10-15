no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
setwd("2_data/RECON3/reaction/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
library(R.matlab)

rm(list = ls())

result1 <-
  readr::read_tsv("recon-store-reactions-1.tsv")

result2 <-
  readr::read_tsv("recon-store-reactions-2.tsv")

result3 <-
  readr::read_tsv("recon-store-reactions-3.tsv")

result <-
  rbind(result1,
        result2,
        result3)

recon3_reaction_data <-
result %>%
  dplyr::rename(reaction_KEGG_ID = keggId,
                reaction_METANETX_ID = metanetx,
                reaction_MODELSEED_ID = seed)


save(recon3_reaction_data, file = "recon3_reaction_data")
