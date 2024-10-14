## http://bigg.ucsd.edu/
##
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())

source("R/3_utils.R")

load("other_files/source_system/source_system.rda")

setwd("other_files/FOODB/")

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
