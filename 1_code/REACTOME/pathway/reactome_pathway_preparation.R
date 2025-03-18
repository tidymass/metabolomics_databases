no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
# source("R/read_sbml_data.R")
dir.create("2_data/REACTOME/pathway/")
setwd("2_data/REACTOME/pathway/")

library(metid)
library(tidyverse)
library(XML)
library(SBMLR)
library(massdatabase)

pathway_info <-
  request_reactome_pathway_info(organism = "Homo sapiens")

dim(pathway_info)

pathway_info$pathway_id

# for(i in 1:length(pathway_info$pathway_id)){
#   cat(i, " ")
#   pathway <-
#     request_reactome_pathway(pathway_id = pathway_info$pathway_id[i])
#   save(pathway, file = paste0(pathway_info$pathway_id[i], ".rda"))
# }


reactome_hsa_pathway <-
  vector(mode = "list", length = length(pathway_info$pathway_id))

for (i in 1:length(pathway_info$pathway_id)) {
  cat(i, " ")
  load(paste0(pathway_info$pathway_id[i], ".rda"))
  pathway$component_info <-
    pathway$component_info %>%
    dplyr::filter(!is.na(node_type)) %>%
    dplyr::filter(node_type == "metabolite")
  
  reactome_hsa_pathway[[i]] <-
    pathway
}

description <-
  reactome_hsa_pathway %>%
  purrr::map(function(x) {
    x$pathway_description
  })

pathway_id <-
  pathway_info$pathway_id

pathway_name <-
  pathway_info$pathway_name

pathway_class <-
  vector(mode = "list", length = length(reactome_hsa_pathway))

database_info <-
  vector(mode = "list", length = 2)

names(database_info) <-
  c("source", "version")

database_info$source <- "reactome_homo_apiensa"
database_info$version <- "2025-2-27"

gene_list <-
  vector(mode = "list", length = length(reactome_hsa_pathway))

for (i in 1:length(gene_list)) {
  gene_list[[i]] <-
    as.data.frame(matrix(nrow = 0, ncol = 0))
}

protein_list <-
  vector(mode = "list", length = length(reactome_hsa_pathway))

for (i in 1:length(protein_list)) {
  protein_list[[i]] <-
    as.data.frame(matrix(nrow = 0, ncol = 0))
}

load("../../../3_data_analysis/HMDB/MS1/hmdb_ms1.rda")
load("../../../3_data_analysis/CHEBI/chebi_ms1.rda")
load("../../../3_data_analysis/KEGG/compound/kegg_compound_ms1.rda")


reactome_hsa_pathway %>%
  purrr::map(function(x) {
    x$component_info$node_id
  }) %>%
  unlist() %>%
  unique()

##metabolite list
metabolite_list <-
  1:length(reactome_hsa_pathway) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <- reactome_hsa_pathway[[i]]
    x <- x$component_info %>%
      dplyr::distinct(node_id, .keep_all = TRUE)
    
    x$node_id <-
      x$node_id %>%
      stringr::str_replace("CHEBI:", "")
    
    ##CAS, ChEBI, Chemspider, DrugBank, InChIKey, LIPID MAPS,
    ####PubChem-compound, Wikidata
    hmdb_kegg1 <-
      x %>%
      dplyr::mutate(row_id = row_number()) %>%
      dplyr::left_join(
        chebi_ms1@spectra.info %>%
          dplyr::mutate(Lab.ID = stringr::str_replace(Lab.ID, "CHEBI:", "")),
        by = c("node_id" = "Lab.ID"),
        na_matches = "never"
      ) %>%
      dplyr::distinct(row_id, .keep_all = TRUE) %>%
      dplyr::select(HMDB.ID, KEGG.ID)
    
    hmdb_kegg2 <-
      x %>%
      dplyr::mutate(row_id = row_number()) %>%
      dplyr::left_join(
        kegg_compound_ms1@spectra.info %>%
          dplyr::mutate(CHEBI.ID = stringr::str_replace(CHEBI.ID, "CHEBI:", "")),
        by = c("node_id" = "CHEBI.ID"),
        na_matches = "never"
      ) %>%
      dplyr::distinct(row_id, .keep_all = TRUE) %>%
      dplyr::select(HMDB.ID, KEGG.ID)
    
    hmdb_kegg3 <-
      x %>%
      dplyr::mutate(row_id = row_number()) %>%
      dplyr::left_join(
        hmdb_ms1@spectra.info %>%
          dplyr::mutate(CHEBI.ID = stringr::str_replace(CHEBI.ID, "CHEBI:", "")),
        by = c("node_id" = "CHEBI.ID"),
        na_matches = "never"
      ) %>%
      dplyr::distinct(row_id, .keep_all = TRUE) %>%
      dplyr::select(HMDB.ID, KEGG.ID)
    
    hmdb_id <-
      data.frame(
        HMDB.ID2 = hmdb_kegg1$HMDB.ID,
        HMDB.ID3 = hmdb_kegg2$HMDB.ID,
        HMDB.ID4 = hmdb_kegg3$HMDB.ID
      ) %>%
      apply(1, function(x) {
        as.character(x[which(!is.na(x))][1])
      })
    
    kegg_id <-
      data.frame(
        KEGG.ID2 = hmdb_kegg1$KEGG.ID,
        KEGG.ID3 = hmdb_kegg2$KEGG.ID,
        KEGG.ID4 = hmdb_kegg3$KEGG.ID
      ) %>%
      apply(1, function(x) {
        as.character(x[which(!is.na(x))][1])
      })
    
    x$HMDB_ID <- hmdb_id
    x$KEGG_ID <- kegg_id
    
    x <-
      x %>%
      dplyr::rename(HMDB.ID = HMDB_ID, KEGG.ID = KEGG_ID)
    colnames(x)[1] <- "CHEBI.ID"
    x
  })

library(metpath)

reactome_hsa_pathway =
  new(
    Class = "pathway_database",
    database_info = database_info,
    pathway_id = pathway_id,
    pathway_name = pathway_name,
    describtion = description,
    pathway_class = pathway_class,
    gene_list = gene_list,
    compound_list = metabolite_list,
    protein_list = list(),
    reference_list = list(),
    related_disease = list(),
    related_module = list()
  )

save(reactome_hsa_pathway, file = "../../../3_data_analysis/REACTOME/pathway/reactome_hsa_pathway.rda")
