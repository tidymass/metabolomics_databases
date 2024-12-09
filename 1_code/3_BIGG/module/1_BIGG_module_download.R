#####this is from the BIGG database
####http://bigg.ucsd.edu/

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

####download the BIGG module databases
setwd("2_data/BIGG/module")

library(massdatabase)

module_info <-
request_bigg_model_info(url = "http://bigg.ucsd.edu/api/v2/models")  

class(module_info)

dim(module_info)

head(module_info)

###108 models in the database

for(i in 1:nrow(module_info)) {
    cat(i, " ")
    download_bigg_model(model_id = module_info$bigg_id[i], path = ".")
}
