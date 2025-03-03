setwd(masstools::get_project_wd())
rm(list = ls())
source("1_code/3_utils.R")
setwd("2_data/WIKIPATHWAYS/")

library(tidyverse)
library(xml2)
library(stringr)
library(rWikiPathways)
listOrganisms()

######Homo sapiens

#####pathways
hsa_pathways <- listPathways('Homo sapiens')
hsa_pathways


dir.create("Homo_sapiensa_pathways")

library(massdatabase)

###download database
# for (i in seq_len(nrow(hsa_pathways))) {
#   cat(i, " ")
#
#   pathway <-
#     request_wikipathway(pathway_id = hsa_pathways$id[i])
#   save(pathway, file = file.path("Homo_sapiensa_pathways", hsa_pathways$id[i]))
#
# }

hsa_pathway_data <-
  seq_len(nrow(hsa_pathways)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    load(file.path("Homo_sapiensa_pathways", hsa_pathways$id[i]))
    pathway$node_info <-
      pathway$node_info %>%
      dplyr::filter(Type == "Metabolite")
    pathway
  })

library(metpath)

description <-
  hsa_pathway_data %>%
  purrr::map(function(x) {
    x$pathway_description
  })

pathway_id <-
  hsa_pathways$id

pathway_name <-
  hsa_pathways$name

pathway_class <-
  vector(mode = "list", length = length(hsa_pathway_data))

database_info <-
  vector(mode = "list", length = 2)

names(database_info) <-
  c("source", "version")

database_info$source <- "wikipathways_homo_apiensa"
database_info$version <- "2025-2-27"

gene_list <-
  vector(mode = "list", length = length(hsa_pathway_data))

for (i in 1:length(gene_list)) {
  gene_list[[i]] <-
    as.data.frame(matrix(nrow = 0, ncol = 0))
}

protein_list <-
  vector(mode = "list", length = length(hsa_pathway_data))

for (i in 1:length(protein_list)) {
  protein_list[[i]] <-
    as.data.frame(matrix(nrow = 0, ncol = 0))
}

load("../../3_data_analysis/HMDB/MS1/hmdb_ms1.rda")
load("../../3_data_analysis/CHEBI/chebi_ms1.rda")
load("../../3_data_analysis/KEGG/compound/kegg_compound_ms1.rda")


hsa_pathway_data %>%
  purrr::map(function(x) {
    unique(x$node_info$Database)
  }) %>%
  unlist() %>%
  unique()

##metabolite list
metabolite_list <-
  1:length(hsa_pathway_data) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <- hsa_pathway_data[[i]]
    x <- x$node_info %>%
      dplyr::distinct(ID, .keep_all = TRUE)
    x$Database[x$Database == "ChEBI compound"] <- "ChEBI"
    x$ID <-
      x$ID %>%
      stringr::str_replace("CHEBI:", "")
    colnames(x) <-
      c("Compound.name", "GraphId", "Type", "Database", "ID")
    
    x <-
      x %>%
      dplyr::mutate(CAS_ID =
                      case_when(Database == "CAS" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(HMDB_ID =
                      case_when(Database == "HMDB" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(ChEBI_ID =
                      case_when(Database == "ChEBI" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(ChEMBL_ID =
                      case_when(Database == "ChEMBL compound" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(Chemspider_ID =
                      case_when(Database == "Chemspider" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(DrugBank =
                      case_when(Database == "DrugBank" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(InChIKey =
                      case_when(Database == "InChIKey" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(
        KEGG_ID =
          case_when(
            Database == "KEGG Compound" ~ ID,
            Database == "KEGG Drug" ~ ID,
            Database == "KEGG Glycan" ~ ID,
            TRUE ~ NA_character_
          )
      ) %>%
      dplyr::mutate(LIPIDMAPS_ID =
                      case_when(Database == "LIPID MAPS" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(LipidBank_ID =
                      case_when(Database == "LipidBank" ~ ID, TRUE ~ NA_character_)) %>%
      dplyr::mutate(
        PubChem_ID =
          case_when(
            Database == "PubChem Compound" ~ ID,
            Database == "PubChem Substance" ~ ID,
            TRUE ~ NA_character_
          )
      ) %>%
      dplyr::mutate(Wikidata_ID =
                      case_when(Database == "Wikidata" ~ ID, TRUE ~ NA_character_))
    
    ##CAS, ChEBI, Chemspider, DrugBank, InChIKey, LIPID MAPS,
    ####PubChem-compound, Wikidata
    hmdb_kegg1 <-
      x %>%
      dplyr::mutate(row_id = row_number()) %>%
      dplyr::left_join(
        chebi_ms1@spectra.info %>%
          dplyr::mutate(Lab.ID = stringr::str_replace(Lab.ID, "CHEBI:", "")),
        by = c("ChEBI_ID" = "Lab.ID"),
        na_matches = "never"
      ) %>%
      dplyr::distinct(row_id, .keep_all = TRUE) %>%
      dplyr::select(HMDB.ID, KEGG.ID)
    
    hmdb_kegg2 <-
      x %>%
      dplyr::mutate(row_id = row_number()) %>%
      dplyr::left_join(hmdb_ms1@spectra.info,
                       by = c("CAS_ID" = "CAS.ID"),
                       na_matches = "never") %>%
      dplyr::distinct(row_id, .keep_all = TRUE) %>%
      dplyr::select(HMDB.ID, KEGG.ID)
    
    hmdb_kegg3 <-
      x %>%
      dplyr::mutate(row_id = row_number()) %>%
      dplyr::left_join(
        hmdb_ms1@spectra.info,
        by = c("PubChem_ID" = "PUBCHEM.ID"),
        na_matches = "never"
      ) %>%
      dplyr::distinct(row_id, .keep_all = TRUE) %>%
      dplyr::select(HMDB.ID, KEGG.ID)
    
    hmdb_id <-
      data.frame(
        HMDB.ID1 = x$HMDB_ID,
        HMDB.ID2 = hmdb_kegg1$HMDB.ID,
        HMDB.ID3 = hmdb_kegg2$HMDB.ID,
        HMDB.ID4 = hmdb_kegg3$HMDB.ID
      ) %>%
      apply(1, function(x) {
        as.character(x[which(!is.na(x))][1])
      })
    
    kegg_id <-
      data.frame(
        KEGG.ID1 = x$KEGG_ID,
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
    x
  })

library(metpath)

wikipathway_hsa_pathway =
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


save(wikipathway_hsa_pathway, file = "../../3_data_analysis/WIKIPATHWAYS/wikipathway_hsa_pathway.rda")




#####reaction
file_name <-
  dir(path = "wikipathways-20220410-gpml-Homo_sapiens", full.names = TRUE)

wikipathway_reaction_database <-
  seq_along(file_name) %>%
  purrr::map(function(i) {
    cat(i, " ")
    read_gpml(file = file_name[i])$reaction
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

wikipathway_reaction_database <-
  wikipathway_reaction_database %>%
  dplyr::filter(
    !Database %in% c(
      "",
      "IntAct",
      "KEGG Pathway",
      "WikiPathways",
      "XMetDB",
      "Database",
      "Wikidata",
      "SPIKE",
      "ChEBI",
      "KEGG Compound",
      "PATO:0002220",
      "Uniprot-TrEMBL"
    )
  ) %>%
  dplyr::distinct()

unique(wikipathway_reaction_database$Database)

save(wikipathway_reaction_database, file = "wikipathway_reaction_database")
