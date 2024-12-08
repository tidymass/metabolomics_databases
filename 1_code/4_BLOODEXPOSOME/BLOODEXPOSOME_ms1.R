###no source
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())

load("3_data_analysis/KEGG/compound/kegg_compound_ms1.rda")
load("3_data_analysis/HMDB/MS1/hmdb_compound_ms1.rda")

source("1_code/3_utils.R")

data <- readxl::read_xlsx("2_data/BLOODEXPOSOME/BloodExpsomeDatabase_version_1.0.xlsx")

dir.create("3_data_analysis/BLOODEXPOSOME", showWarnings = FALSE)
setwd("3_data_analysis/BLOODEXPOSOME")

library(tidyverse)
library(tidymass)

data <-
data %>%
  dplyr::rename(
    PUBCHEM.ID = "PubChem CID",
    Compound.name = "Compound Name",
    KEGG.ID = "KEGG ID",
    HMDB.ID = "HMDB ID",
    Formula = "Molecular Formula",
    SMILES.ID = CanonicalSMILES,
    INCHIKEY.ID = InChIKey,
    mz = ExactMass
  ) %>%
  dplyr::select(PUBCHEM.ID:BloodPaperCount)

sum(is.na(data$KEGG.ID))
sum(is.na(data$HMDB.ID))

data$KEGG.ID <-
data$KEGG.ID %>%
  stringr::str_split(pattern = ";") %>%
  purrr::map(function(x){
    if(is.na(x)[1]){
      return(NA)
    }
    paste(x, collapse = "{}")
  }) %>%
  unlist()

data$HMDB.ID <-
  data$HMDB.ID %>%
  stringr::str_split(pattern = ";") %>%
  purrr::map(function(x){
    if(is.na(x)[1]){
      return(NA)
    }
    paste(x, collapse = "{}")
  }) %>%
  unlist()

data <-
data %>%
  dplyr::left_join(kegg_compound_ms1@spectra.info[,c("KEGG.ID", "from_human", "from_which_part","from_drug", "from_which_drug")],
                   by = "KEGG.ID",
                   na_matches = "never") %>%
  dplyr::left_join(hmdb_compound_ms1@spectra.info[,c("HMDB.ID", "from_human", "from_which_part")],
                   by = "HMDB.ID",
                   na_matches = "never")


# data$From_human <-
# data[,c("From_human.x", "From_human.y")] %>%
#   apply(1, function(x){
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if(length(x) == 0){
#       return(NA)
#     }
#     if(length(x) == 2){
#       return("YES")
#     }else{
#       return(x)
#     }
#   })


data$from_human <-
data[,c("from_human.x", "from_human.y")] %>%
  apply(1, function(x){
    x <- as.character(x)
    x <- x[!is.na(x)]
    if(length(x) == 0){
      return(NA)
    }
    if(length(x) == 2){
      return(any(as.character(x)))
    }else{
      return(x)
    }
  })


data <-
  data %>%
  dplyr::select(-c(From_human.x, From_human.y, XLogP))

openxlsx::write.xlsx(data, file = "bloodexposome.xlsx", asTable = TRUE)

bloodexposome_ms1 <-
  construct_database(
    path = ".",
    version = "2022-04-25",
    metabolite.info.name = "bloodexposome.xlsx",
    source = "BLOODEXPOSOME",
    link = "https://bloodexposome.org/#/dashboard",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

bloodexposome_ms1

load("../HMDB/MS1/hmdb_compound_ms1.rda")
load("../KEGG/kegg_compound_ms1.rda")

intersect(colnames(bloodexposome_ms1@spectra.info),
          colnames(hmdb_compound_ms1@spectra.info))

setdiff(colnames(hmdb_compound_ms1@spectra.info),
        colnames(bloodexposome_ms1@spectra.info))

bloodexposome_ms1 <-
  update_metid_database_info(
    database = bloodexposome_ms1,
    ref_database = hmdb_compound_ms1,
    by = c(
      "PUBCHEM.ID",
      "Compound.name",
      "KEGG.ID",
      "HMDB.ID",
      "Formula",
      "SMILES.ID",
      "INCHIKEY.ID"
    ),
    combine_columns = c(
      "PUBCHEM.ID",
      "Compound.name",
      "KEGG.ID",
      "HMDB.ID",
      "Formula",
      "SMILES.ID",
      "INCHIKEY.ID"
    ),
    new_columns = c(
      "CAS.ID",
      "mz.pos",
      "mz.neg",
      "Submitter",
      "version",
      "Create_date",
      "Updated_date",
      "status",
      "secondary_accessions",
      "Description",
      "Synonyms",
      "monisotopic_molecular_weight",
      "IUPAC_name",
      "Traditional_IUPAC_name",
      "INCHI.ID",
      "Kingdom",
      "Super_class",
      "Class",
      "Sub_class",
      "State",
      "Biospecimen_locations",
      "Cellular_locations",
      "Tissue_locations",
      "CHEMSPIDER.ID",
      "DRUGBANK.ID",
      "FOODB.ID",
      "CHEBI.ID",
      "BIOCYC.ID",
      "BIGG.ID",
      "WIKIPEDIA.ID",
      "METLIN.ID"
    )
  )

intersect(colnames(bloodexposome_ms1@spectra.info),
          colnames(kegg_compound_ms1@spectra.info))

setdiff(colnames(kegg_compound_ms1@spectra.info),
        colnames(bloodexposome_ms1@spectra.info))

source(here::here("1_code/3_utils.R"))

bloodexposome_ms1 <-
  update_metid_database_info(
    database = bloodexposome_ms1,
    ref_database = kegg_compound_ms1,
    by = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "BIGG.ID",
      "BIOCYC.ID",
      "CHEBI.ID",
      "INCHIKEY.ID"
    ),
    combine_columns = c(
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "BIGG.ID",
      "BIOCYC.ID",
      "CHEBI.ID",
      "INCHIKEY.ID"
    ),
    new_columns = c(
      "CHEMBL.ID",
      "LIPIDMAPS.ID",
      "LIPIDBANK.ID",
      "KEGG_DRUG.ID"
    )
  )

load(here::here("2_data/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

bloodexposome_ms1 <-
  update_metid_database_source_system(database = bloodexposome_ms1,
                                      source_system = source_system,
                                      by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
                                      prefer = "database")

save(bloodexposome_ms1, file = "bloodexposome_ms1.rda")




save(bloodexposome_ms1, file = "bloodexposome_ms1.rda", compress = "xz")
