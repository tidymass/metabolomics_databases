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
