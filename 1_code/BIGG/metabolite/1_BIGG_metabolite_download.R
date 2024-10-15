#####this is from the BIGG database
####http://bigg.ucsd.edu/

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

dir.create("2_data/BIGG/metabolite", showWarnings = FALSE)

##download data

# massdatabase::download_bigg_universal_metabolite(path = ".", 
# sleep = 1)

###read data
bigg_metabolites <-
read_bigg_universal_metabolite(path = ".")

###convert it to metID compound database
bigg_ms1 <-
convert_bigg_universal2metid(data = bigg_metabolites, 
path = ".", 
threads = 5)

setwd(get_project_wd())
setwd("3_data_analysis/BIGG/metabolite")
save(bigg_ms1, file = "bigg_ms1.rda")