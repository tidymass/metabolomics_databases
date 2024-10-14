###
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)
library(plyr)
library(tidyr)
source("R/6_CHEBI.R")
source("R/3_utils.R")

setwd("other_files/CHEBI/")

# # download_chebi_compound(path = ".")
# #
chebi_compound <-
  read_chebi_compound(path = ".")
# # save(chebi_compound, file = "chebi_compound.rda", compress = "xz")
#
# load("chebi_compound.rda")
#
# chebi_compound <-
#   chebi_compound %>%
#   dplyr::rename(
#     Lab.ID = COMPOUND_ID,
#     Formula = FORMULA,
#     mz = MONOISOTOPIC_MASS,
#     Charge = CHARGE,
#     INCHI.ID = InChI,
#     Species = SPECIES,
#     Species.ID = SPECIES_ACCESSION,
#     Status = STATUS,
#     Database_source = SOURCE,
#     Compound.name = NAME,
#     Definition = DEFINITION,
#     Star = STAR,
#     Updated_date = MODIFIED_ON
#   ) %>%
#   dplyr::mutate(
#     CHEBI.ID = Lab.ID,
#     RT = NA,
#     mz.pos = NA,
#     mz.neg = NA,
#     Submitter = "CHEBI"
#   ) %>%
#   dplyr::select(
#     Lab.ID,
#     Compound.name,
#     mz,
#     RT,
#     CAS.ID,
#     HMDB.ID,
#     KEGG.ID,
#     Formula,
#     mz.pos,
#     mz.neg,
#     Submitter,
#     everything()
#   )
#
# chebi_compound <-
#   as.data.frame(chebi_compound)
#
# chebi_compound[which(chebi_compound == "null", arr.ind = TRUE)] <-
#   NA
#
# chebi_compound$mz <- as.numeric(chebi_compound$mz)
#
# ###remove mz is NA
# chebi_compound <-
#   chebi_compound %>%
#   dplyr::filter(!is.na(mz))
#
# dim(chebi_compound)
#
# library(metid)
#
# # temp <-
# # chebi_compound %>%
# #   purrr::map(function(x){
# #     x = stringi::stri_enc_toutf8(x)
# #     x
# #   }) %>%
# #   dplyr::bind_cols()
# #
# # invalid_utf8 <- function(x){
# #   !is.na(x) & is.na(iconv(x, "UTF-8", "UTF-8"))
# # }
# #
# # which(invalid_utf8(chebi_compound$Lab.ID))
# # which(invalid_utf8(chebi_compound$Compound.name))
# # which(invalid_utf8(chebi_compound$mz))
# # which(invalid_utf8(chebi_compound$RT))
# # which(invalid_utf8(chebi_compound$CAS.ID))
#
# species <-
#   sort(unique(unlist(
#     stringr::str_split(chebi_compound$Species[!is.na(chebi_compound$Species)], "\\{\\}")
#   )))
#
# species.id <-
#   sort(unique(unlist(
#     stringr::str_split(chebi_compound$Species.ID[!is.na(chebi_compound$Species.ID)], "\\{\\}")
#   )))
#
# stringr::str_split(species.id, ":") %>% lapply(function(x)
#   x[1]) %>% unlist() %>% table()
#
# library("taxize")
# gnr_resolve(species[1])
#
# # eol_id <-
# #   get_eolid(species[1:5], ask = TRUE)
# #
# # eol_id <-
# #   get_eolid(species[1], ask = FALSE)
# #
# # eol_id %>%
# #   dplyr::filter(source == "wikipedia EN")
#
# # wiki_class <- vector(mode = "list", length = length(species))
# #
# # for(i in 1:length(wiki_class)){
# #   cat(i, " ")
# #   wiki_class[[i]] <-
# #     classification(species[i], db = 'wiki',
# #                    wiki_site = "species",
# #                    wiki = "en",
# #                    limit = 1)[[1]]
# # }
# #
# # save(wiki_class, file = "wiki_class")
# #
# # load("wiki_class")
# #
# # names(wiki_class) <- species
# #
# # wiki_class_info <-
# #   lapply(wiki_class, function(x) {
# #     if (all(is.na(x))) {
# #       return(NA)
# #     }
# #
# #     x <-
# #       x %>%
# #       dplyr::filter(!is.na(name))
# #
# #     if (any(x$name == "Bacteria")) {
# #       return("Bacteria")
# #     }
# #
# #     if (any(x$name == "Fungi")) {
# #       return("Fungi")
# #     }
# #
# #     if (any(x$name == "Archaea")) {
# #       return("Archaea")
# #     }
# #
# #     if (any(x$name == "Plantae")) {
# #       return("Plantae")
# #     }
# #
# #     if (any(x$name == "H. sapiens")) {
# #       return("Human")
# #     }
# #
# #     if (any(x$name == "Animalia")) {
# #       return("Animalia")
# #     }
# #
# #     return(x$name[1])
# #   }) %>%
# #   unlist()
# #
# # wiki_class_info
# #
# # wikipedia_class <- vector(mode = "list", length = length(species))
# # names(wikipedia_class) <- species
# #
# # for(i in 1:length(species)){
# #   cat(i, " ")
# #   wikipedia_class[[i]] <-
# #     request_wikipedia_scientific_classification(species_id = species[i])
# # }
# #
# # save(wikipedia_class, file = "wikipedia_class")
# # load("wikipedia_class")
# #
# # wikipedia_class_info <-
# #   lapply(wikipedia_class, function(x) {
# #     if (all(is.na(x))) {
# #       return(NA)
# #     }
# #
# #     if (any(x$value == "Bacteria")) {
# #       return("Bacteria")
# #     }
# #
# #     if (any(x$value == "Fungi")) {
# #       return("Fungi")
# #     }
# #
# #     if (any(x$value == "Archaea")) {
# #       return("Archaea")
# #     }
# #
# #     if (any(x$value == "Plantae")) {
# #       return("Plantae")
# #     }
# #
# #     if (any(x$value == "H. sapiens")) {
# #       return("Human")
# #     }
# #
# #     if (any(x$value == "Animalia")) {
# #       return("Animalia")
# #     }
# #
# #     return(x$value[1])
# #   }) %>%
# #   unlist()
# #
# # class_info <-
# #   data.frame(species = species,
# #              wikipedia = wikipedia_class_info,
# #              wiki = wiki_class_info)
# #
# # merge <-
# #   class_info[, c(2, 3)] %>%
# #   apply(1, function(x) {
# #     if (all(is.na(x))) {
# #       return(NA)
# #     }
# #     x <- x[!is.na(x)]
# #     x[1]
# #   })
# #
# # class_info$merge = merge
# #
# #
# # write.csv(class_info, "class_info.csv", row.names = FALSE)
#
# class_info <- readr::read_csv("class_info_manual.csv")
#
# match_table <-
#   class_info %>%
#   dplyr::select(species, merge)
#
# source <-
#   chebi_compound$Species %>%
#   purrr::map(function(x) {
#     convert_species2source(x, match_table)
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# dim(chebi_compound)
# dim(source)
#
# source$From_drug[which(!is.na(chebi_compound$KEGG_DRUG.ID))] <-
#   'Yes'
#
# chebi_compound <-
#   cbind(chebi_compound,
#         source)
#
# chebi_compound$From_drug[which(!is.na(chebi_compound$DRUGBANK.ID))] <-
#   'Yes'
#
# readr::write_csv(chebi_compound, file = "chebi_compound.csv")

chebi_ms1 <-
  construct_database(
    path = ".",
    version = "2022-04-11",
    metabolite.info.name = "chebi_compound.csv",
    source = "CHEBI",
    link = "https://www.ebi.ac.uk/chebi/init.do",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

chebi_ms1

load("../HMDB/MS1/hmdb_ms1.rda")
load("../KEGG/kegg_ms1.rda")

intersect(colnames(chebi_ms1@spectra.info),
          colnames(hmdb_ms1@spectra.info))

setdiff(colnames(hmdb_ms1@spectra.info),
        colnames(chebi_ms1@spectra.info))

chebi_ms1 <-
  update_metid_database_info(
    database = chebi_ms1,
    ref_database = hmdb_ms1,
    by = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "INCHI.ID",
      "FOODB.ID",
      "WIKIPEDIA.ID",
      "DRUGBANK.ID",
      "PUBCHEM.ID",
      "CHEBI.ID"
    ),
    combine_columns = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "INCHI.ID",
      "FOODB.ID",
      "WIKIPEDIA.ID",
      "DRUGBANK.ID",
      "PUBCHEM.ID",
      "CHEBI.ID"
    ),
    new_columns = c(
      "status",
      "IUPAC_name",
      "Description",
      "Traditional_IUPAC_name",
      "SMILES.ID",
      "INCHIKEY.ID",
      "Kingdom",
      "Super_class",
      "Class",
      "Sub_class",
      "State",
      "Biospecimen_locations",
      "Cellular_locations",
      "Tissue_locations",
      "CHEMSPIDER.ID",
      "BIOCYC.ID",
      "BIGG.ID",
      "METLIN.ID"
    )
  )

intersect(colnames(chebi_ms1@spectra.info),
          colnames(kegg_ms1@spectra.info))

setdiff(colnames(kegg_ms1@spectra.info),
        colnames(chebi_ms1@spectra.info))

source(here::here("R/3_utils.R"))

chebi_ms1 <-
  update_metid_database_info(
    database = chebi_ms1,
    ref_database = kegg_ms1,
    by = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "INCHI.ID",
      "FOODB.ID",
      "WIKIPEDIA.ID",
      "DRUGBANK.ID",
      "PUBCHEM.ID",
      "CHEBI.ID"
    ),
    combine_columns = c(
      "Compound.name",
      "CAS.ID",
      "HMDB.ID",
      "KEGG.ID",
      "INCHI.ID",
      "FOODB.ID",
      "WIKIPEDIA.ID",
      "DRUGBANK.ID",
      "PUBCHEM.ID",
      "CHEBI.ID"
    ),
    new_columns = c(
      "CHEMBL.ID",
      "LIPIDBANK.ID"
    )
  )

load(here::here("other_files/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

chebi_ms1 <-
  update_metid_database_source_system(database = chebi_ms1,
                                      source_system = source_system,
                                      by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
                                      prefer = "database")

save(chebi_ms1, file = "chebi_ms1.rda")

