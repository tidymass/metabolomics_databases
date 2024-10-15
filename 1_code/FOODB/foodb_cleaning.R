## http://bigg.ucsd.edu/
##
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())

source("1_code/3_utils.R")

load("2_data/source_system/source_system.rda")

setwd("2_data/FOODB/")

library(tidyverse)
library(tidyselect)
library(metid)

load("foodb_ms1.rda")

foodb_ms1 <-
  update_metid_database_source_system(
    database = foodb_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(foodb_ms1, file = "foodb_ms1.rda")
