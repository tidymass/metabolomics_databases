## http://bigg.ucsd.edu/
##
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())

source("1_code/3_utils.R")

load("2_data/source_system/source_system.rda")

setwd("2_data/T3DB/")

library(tidyverse)
library(tidyselect)
library(metid)

load("t3db_ms1.rda")

t3db_ms1 <-
  update_metid_database_source_system(
    database = t3db_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(t3db_ms1, file = "t3db_ms1.rda")

