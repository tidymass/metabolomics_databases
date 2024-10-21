no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
## install.packages("MetaDBparse")
# library(MetaDBparse)

rm(list = ls())

### to get KEGG database
library(KEGGgraph)
library(KEGGREST)
library(tidyverse)

dir.create("2_data/KEGG/drug", showWarnings = FALSE)
setwd("2_data/KEGG/drug")

# drug_ID <-
#   keggList(database = "drug") %>%
#   names() %>%
#   unique() %>%
#   stringr::str_replace_all(., "cpd:", "")
#
# kegg_drug_database <-
#   pbapply::pblapply(drug_ID, function(x){
#     KEGGREST::keggGet(dbentries = x)[[1]]
#   })
#
# save(kegg_drug_database, file = "kegg_drug_database")

load("kegg_drug_database")
