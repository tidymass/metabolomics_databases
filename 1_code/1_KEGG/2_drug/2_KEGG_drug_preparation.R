no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

rm(list = ls())

### to get KEGG database
library(KEGGgraph)
library(KEGGREST)
library(tidyverse)

### load data
load("2_data/KEGG/drug/kegg_drug_database")

dir.create("3_data_analysis/KEGG/drug", showWarnings = FALSE)
setwd("3_data_analysis/KEGG/drug")

kegg_drug_database %>%
  lapply(function(x){
    stringr::str_split(x$DBLINK, ": ") %>%
      lapply(function(x){
        x[1]
      }) %>%
      unlist()
  }) %>%
  unlist() %>%
  unique()

kegg_drug <-
  kegg_drug_database %>%
  purrr::map(function(x) {
    cat(x$ENTRY, " ")
    KEGG_DRUG.ID = x$ENTRY
    x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
    Compound.name = paste(x$NAME, collapse = "{}")
    Formula = x$FORMULA
    if (is.null(x$FORMULA)) {
      Formula = NA
    }
    mz = as.numeric(x$EXACT_MASS)
    if (is.null(x$EXACT_MASS)) {
      mz = NA
    }

    CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
      stringr::str_trim(side = "both")

    PUBCHEM.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
      stringr::str_trim(side = "both")

    CHEBI.ID <-
      stringr::str_replace(grep("ChEBI", x$DBLINKS, value = TRUE), "ChEBI: ", "") %>%
      stringr::str_trim(side = "both")

    CHEMBL.ID <-
      stringr::str_replace(grep("ChEMBL", x$DBLINKS, value = TRUE), "ChEMBL: ", "") %>%
      stringr::str_trim(side = "both")

    LIPIDMAPS.ID <-
      stringr::str_replace(grep("LIPIDMAPS", x$DBLINKS, value = TRUE), "LIPIDMAPS: ", "") %>%
      stringr::str_trim(side = "both")

    LIPIDBANK.ID <-
      stringr::str_replace(grep("LipidBank", x$DBLINKS, value = TRUE), "LipidBank: ", "") %>%
      stringr::str_trim(side = "both")

    DRUGBANK.ID <-
      stringr::str_replace(grep("DrugBank", x$DBLINKS, value = TRUE), "DrugBank: ", "") %>%
      stringr::str_trim(side = "both")

    if (length(CAS.ID) == 0) {
      CAS.ID = NA
    }

    if (length(PUBCHEM.ID) == 0) {
      PUBCHEM.ID = NA
    }

    if (length(CHEBI.ID) == 0) {
      CHEBI.ID = NA
    }

    if (length(CHEMBL.ID) == 0) {
      CHEMBL.ID = NA
    }

    if (length(LIPIDMAPS.ID) == 0) {
      LIPIDMAPS.ID = NA
    }

    if (length(LIPIDBANK.ID) == 0) {
      LIPIDBANK.ID = NA
    }

    if (length(DRUGBANK.ID) == 0) {
      DRUGBANK.ID = NA
    }

    From_drug = "Yes"
    REMARK <- x$REMARK
    if(is.null(REMARK)){
      From_human <- "No"
      KEGG.ID <- NA
    }else{
      KEGG.ID <-
        stringr::str_extract_all(REMARK, "Same as: C[0-9]{5,6}")[[1]] %>%
        stringr::str_replace_all("Same as: ", "") %>%
        paste(collapse = "{}")

      if(KEGG.ID == ""){
        From_human <- "No"
        KEGG.ID <- NA
      }else{
        From_human <- "Yes"
        KEGG.ID <- KEGG.ID
      }
    }

    data.frame(
      Lab.ID = KEGG_DRUG.ID,
      Compound.name,
      Formula,
      mz,
      CAS.ID,
      HMDB.ID = NA,
      KEGG.ID,
      PUBCHEM.ID,
      CHEBI.ID,
      CHEMBL.ID,
      LIPIDMAPS.ID,
      LIPIDBANK.ID,
      DRUGBANK.ID = DRUGBANK.ID,
      From_human = From_human,
      From_drug = From_drug,
      KEGG_DRUG.ID = KEGG_DRUG.ID
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

kegg_drug <-
  kegg_drug %>%
  dplyr::filter(!is.na(mz)) %>%
  dplyr::mutate(Synonyms = Compound.name) %>%
  dplyr::mutate(
    RT = NA,
    mz.pos = NA,
    mz.neg = NA,
    Submitter = "KEGG"
  ) %>%
  dplyr::select(
    Lab.ID,
    Compound.name,
    mz,
    RT,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Formula,
    mz.pos,
    mz.neg,
    Submitter,
    everything()
  )

kegg_drug$Compound.name =
  kegg_drug$Compound.name %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist() %>%
  stringr::str_replace(pattern = ";", "")

kegg_drug[which(kegg_drug == "", arr.ind = TRUE)] <- NA


dim(kegg_metabolite)
# dim(kegg_drug)

kegg_metabolite[which(kegg_metabolite == "", arr.ind = TRUE)] <- NA
# kegg_drug[which(kegg_drug == "", arr.ind = TRUE)] <- NA
#
# kegg_drug.id <-
# kegg_metabolite %>%
#   dplyr::filter(!is.na(KEGG_DRUG.ID)) %>%
#   dplyr::pull(KEGG_DRUG.ID) %>%
#   stringr::str_split("\\{\\}") %>%
#   unlist() %>%
#   unique()
#
# kegg_drug <-
#   kegg_drug %>%
#   dplyr::filter(!KEGG_DRUG.ID %in% kegg_drug.id)
#
# kegg.id <-
# kegg_drug$KEGG.ID[!is.na(kegg_drug$KEGG.ID)]
#
# kegg_metabolite <-
# kegg_metabolite %>%
#   dplyr::filter(!KEGG.ID %in% kegg.id)

kegg_metabolite$From_human
kegg_metabolite$From_drug

kegg_metabolite <-
kegg_metabolite %>%
    dplyr::rename(
        from_human = From_human,
        from_drug = From_drug
    ) %>%
    dplyr::mutate(from_human = case_when(
        from_human == "Yes" ~ TRUE,
        from_human == "No" ~ FALSE,
        TRUE ~ NA
    )) %>%
    dplyr::mutate(from_drug = case_when(
        from_drug == "Yes" ~ TRUE,
        from_drug == "No" ~ FALSE,
        TRUE ~ NA
    )) %>%
    dplyr::mutate(
        from_which_part = NA,
        from_which_drug = NA
    )

kegg_metabolite$from_human
kegg_metabolite$from_drug

# sum(kegg_metabolite$From_drug == "Yes")
# sum(kegg_drug$From_human == "Yes")
# sum(kegg_drug$From_drug == "Yes")

colnames(kegg_metabolite)
# colnames(kegg_drug)

# kegg_ms1 <-
#   rbind(kegg_metabolite,
#         kegg_drug)

openxlsx::write.xlsx(kegg_metabolite, file = "kegg_metabolite.xlsx", asTable = TRUE, overwrite = TRUE)

library(metid)

kegg_metabolite_ms1 <-
    construct_database(
        path = ".",
        version = "2022-04-08",
        metabolite.info.name = "kegg_metabolite.xlsx",
        source = "KEGG",
        link = "https://www.genome.jp/kegg",
        creater = "Xiaotao Shen",
        email = "xiaotao.shen@outlook.com",
        rt = FALSE,
        threads = 3
    )

# load("../HMDB/MS1/hmdb_ms1.rda")

# intersect(
#     colnames(kegg_ms1@spectra.info),
#     colnames(hmdb_ms1@spectra.info)
# )
# setdiff(
#     colnames(hmdb_ms1@spectra.info),
#     colnames(kegg_ms1@spectra.info)
# )

# source(here::here("1_code/3_utils.R"))

# kegg_ms1 <-
#     update_metid_database_info(
#         database = kegg_ms1,
#         ref_database = hmdb_ms1,
#         by = c(
#             "Compound.name",
#             "CAS.ID",
#             "HMDB.ID",
#             "KEGG.ID",
#             "PUBCHEM.ID",
#             "CHEBI.ID",
#             "DRUGBANK.ID"
#         ),
#         combine_columns = c(
#             "CAS.ID",
#             "HMDB.ID",
#             "KEGG.ID",
#             "PUBCHEM.ID",
#             "CHEBI.ID",
#             "DRUGBANK.ID",
#             "From_human"
#         ),
#         new_columns = c(
#             "status",
#             "IUPAC_name",
#             "Traditional_IUPAC_name",
#             "SMILES.ID",
#             "INCHI.ID",
#             "INCHIKEY.ID",
#             "Kingdom",
#             "Super_class",
#             "Class",
#             "Sub_class",
#             "State",
#             "Biospecimen_locations",
#             "Cellular_locations",
#             "Tissue_locations",
#             "CHEMSPIDER.ID",
#             "FOODB.ID",
#             "BIOCYC.ID",
#             "BIGG.ID",
#             "WIKIPEDIA.ID",
#             "METLIN.ID"
#         )
#     )

# load(here::here("2_data/source_system/source_system.rda"))

# library(tidyverse)
# library(tidyselect)
# library(metid)

# kegg_ms1 <-
#     update_metid_database_source_system(
#         database = kegg_ms1,
#         source_system = source_system,
#         by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
#         prefer = "database"
#     )



# save(kegg_ms1, file = "kegg_ms1.rda")

# load("kegg_ms1.rda")

# idx <-
#     grep("C[0-9]{4,8}", kegg_ms1@spectra.info$Compound.name, value = FALSE)


# if (length(idx) > 0) {
#     kegg_ms1@spectra.info$Compound.name[idx] <-
#         kegg_ms1@spectra.info$Synonyms[idx] %>%
#         stringr::str_split("\\{\\}") %>%
#         purrr::map(function(x) {
#             x[2]
#         }) %>%
#         unlist()
# }

# save(kegg_ms1, file = "kegg_ms1.rda")
