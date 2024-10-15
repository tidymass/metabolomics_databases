#####this is from the BIGG database
####http://bigg.ucsd.edu/

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

###read data
library(massdatabase)

bigg_metabolites <-
read_bigg_universal_metabolite(path = "2_data/BIGG/metabolite")

dir.create("3_data_analysis/BIGG/metabolite", showWarnings = FALSE)
setwd("3_data_analysis/BIGG/metabolite")

###convert it to metID compound database
class(bigg_metabolites)

####source
####from_human TRUE FALSE
#####from_which_part blood, tissue, urine, cell
####from_bacteria TRUE FALSE
#####from_which_bacteria Escherichia coli str. K-12 substr NA
######from_fungi TRUE FALSE
####from_which_fungi Saccharomyces cerevisiae S288C
####from_eukaryote TRUE FALSE
#####from_which_eukaryote Saccharomyces cerevisiae S288C
####from_archaea TRUE FALSE
#####from_which_archaea Methanosarcina barkeri str. Fusaro
####from_plant TRUE FALSE
#####from_which_plant Arabidopsis thaliana
#####fropm_animal TRUE FALSE
######from_which_animal Mus musculus
#######from_environment TRUE FALSE
#######from_which_environment soil, water, air
#######from_virus TRUE FALSE
#######from_which_virus Human immunodeficiency virus 1
#######from_protist TRUE FALSE
#######from_which_protist Plasmodium falciparum 3D7
#######from_drug TRUE FALSE
#######from_which_drug Aspirin
#######from_food TRUE FALSE
#######from_which_food Tomato

convert_bigg_universal2metid(data = bigg_metabolites, 
path = ".", 
threads = 5)

load("bigg_ms1")

dim(bigg_ms1)

colnames(bigg_ms1@spectra.info)

module_info <-
readr::read_csv("../module/model_info_manual.csv")

spectra_info <-
bigg_ms1@spectra.info

spectra_info <-
spectra_info %>% 
dplyr::distinct(Lab.ID, .keep_all = TRUE)

bigg_ms1@spectra.info <-
spectra_info

source <-
data.frame(Lab.ID = spectra_info$Lab.ID,
from_human = rep(FALSE, nrow(spectra_info)),
from_which_part = rep(NA, nrow(spectra_info)),
from_bacteria = rep(FALSE, nrow(spectra_info)),
from_which_bacteria = rep(NA, nrow(spectra_info)),
from_fungi = rep(FALSE, nrow(spectra_info)),
from_which_fungi = rep(NA, nrow(spectra_info)),
from_eukaryote = rep(FALSE, nrow(spectra_info)),
from_which_eukaryote = rep(NA, nrow(spectra_info)),
from_archaea = rep(FALSE, nrow(spectra_info)),
from_which_archaea = rep(NA, nrow(spectra_info)),
from_plant = rep(FALSE, nrow(spectra_info)),
from_which_plant = rep(NA, nrow(spectra_info)),
from_animal = rep(FALSE, nrow(spectra_info)),
from_which_animal = rep(NA, nrow(spectra_info)),
from_environment = rep(FALSE, nrow(spectra_info)),
from_which_environment = rep(NA, nrow(spectra_info)),
from_virus = rep(FALSE, nrow(spectra_info)),
from_which_virus = rep(NA, nrow(spectra_info)),
from_drug = rep(FALSE, nrow(spectra_info)),
from_which_drug = rep(NA, nrow(spectra_info)),
from_food = rep(FALSE, nrow(spectra_info)),
from_which_food = rep(NA, nrow(spectra_info))
)

for(i in 1:nrow(source)){
    cat(i, " ")
    temp_organism <-
    spectra_info$organism[i] %>% 
    stringr::str_split("\\{\\}") %>% 
    `[[`(1) %>% 
    unique()

temp_module_info <-
module_info %>% 
dplyr::filter(organism %in% temp_organism)

##human
if("Human" %in% temp_module_info$class){
    source$from_human[i] <- TRUE
    source$from_which_part[i] <- NA
}

if("Bacteria" %in% temp_module_info$class){
    source$from_bacteria[i] <- TRUE
    source$from_which_bacteria[i] <- 
    temp_module_info$organism[temp_module_info$class == "Bacteria"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Archaea" %in% temp_module_info$class){
    source$from_archaea[i] <- TRUE
    source$from_which_archaea[i] <- 
    temp_module_info$organism[temp_module_info$class == "Archaea"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Fungi" %in% temp_module_info$class){
    source$from_fungi[i] <- TRUE
    source$from_which_fungi[i] <- 
    temp_module_info$organism[temp_module_info$class == "Fungi"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Eukaryote" %in% temp_module_info$class){
    source$from_eukaryote[i] <- TRUE
    source$from_which_eukaryote[i] <- 
    temp_module_info$organism[temp_module_info$class == "Eukaryote"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}


if("Plant" %in% temp_module_info$class){
    source$from_plant[i] <- TRUE
    source$from_which_plant[i] <- 
    temp_module_info$organism[temp_module_info$class == "Plant"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Animal" %in% temp_module_info$class){
    source$from_animal[i] <- TRUE
    source$from_which_animal[i] <- 
    temp_module_info$organism[temp_module_info$class == "Animal"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Environment" %in% temp_module_info$class){
    source$from_environment[i] <- TRUE
    source$from_which_environment[i] <- 
    temp_module_info$organism[temp_module_info$class == "Environment"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Virus" %in% temp_module_info$class){
    source$from_virus[i] <- TRUE
    source$from_which_virus[i] <- 
    temp_module_info$organism[temp_module_info$class == "Virus"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Mammal" %in% temp_module_info$class){
    source$from_animal[i] <- TRUE
    source$from_which_animal[i] <- 
    temp_module_info$organism[temp_module_info$class == "Mammal"] %>% 
    unique() %>% 
    paste(collapse = "{}")

}

if("Drug" %in% temp_module_info$class){
    source$from_drug[i] <- TRUE
    source$from_which_drug[i] <- 
    temp_module_info$organism[temp_module_info$class == "Drug"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

if("Food" %in% temp_module_info$class){
    source$from_food[i] <- TRUE
    source$from_which_food[i] <- 
    temp_module_info$organism[temp_module_info$class == "Food"] %>% 
    unique() %>% 
    paste(collapse = "{}")
}

}

spectra_info <-
spectra_info %>% 
dplyr::left_join(source, by = "Lab.ID")

bigg_ms1@spectra.info <-
spectra_info

save(bigg_ms1, file = "bigg_ms1.RData")


spectra_info %>% 
dplyr::count(Lab.ID) %>% 
dplyr::filter(n > 1)

spectra_info %>% 
dplyr::filter(Lab.ID == "Pald")
