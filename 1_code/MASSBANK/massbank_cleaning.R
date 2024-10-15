## http://bigg.ucsd.edu/
##
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())

source("1_code/3_utils.R")

load("2_data/source_system/source_system.rda")

setwd("2_data/MASSBANK")

library(tidyverse)
library(tidyselect)
library(metid)

load("massbank_ms2.rda")

massbank_ms2 <-
  update_metid_database_source_system(
    database = massbank_ms2,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(massbank_ms2, file = "massbank_ms2.rda")

