#####this is from the BIGG database
####http://bigg.ucsd.edu/

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(massdatabase)

dir.create("3_data_analysis/BIGG/module", showWarnings = FALSE)
setwd("3_data_analysis/BIGG/module")

module_info <-
request_bigg_model_info(url = "http://bigg.ucsd.edu/api/v2/models")  

# for(i in 1:nrow(module_info)){
#     cat(i, " ")
#     massdatabase::read_bigg_model(path = "../../")
# }



