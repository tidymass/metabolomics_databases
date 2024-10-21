setwd(masstools::get_project_wd())
rm(list = ls())

source("1_code/3_utils.R")
setwd("2_data/HMDB/protein/")

library(tidyverse)
library(xml2)
library(stringr)

hmdb <- read_xml("hmdb_proteins.xml")

hmdb <- as_list(hmdb)

x = hmdb$hmdb[[1]]

hmdb_protein <- vector(mode = "list", length = length(hmdb$hmdb))

for(i in 1:length(hmdb$hmdb)){
  if((i/10) %in% seq(1,100000, 1)){
    cat(i, " ")
  }

  x <- hmdb$hmdb[[i]]

  accession <-
    paste(stringr::str_trim(unlist(x$accession)),
          collapse = "{}")

  name <-
    paste(stringr::str_trim(unlist(x$name)),
          collapse = "{}")

  protein_type <-
    paste(stringr::str_trim(unlist(x$protein_type)),
          collapse = "{}")

  synonyms <-
    paste(stringr::str_trim(unlist(x$synonyms)),
          collapse = "{}")

  gene_name <-
    paste(stringr::str_trim(unlist(x$gene_name)),
          collapse = "{}")


  general_function <-
    paste(stringr::str_trim(unlist(x$general_function)),
          collapse = "{}")

  specific_function <-
    paste(stringr::str_trim(unlist(x$specific_function)),
          collapse = "{}")

  genbank_protein_id <-
    paste(stringr::str_trim(unlist(x$genbank_protein_id)),
          collapse = "{}")

  uniprot_id <-
    paste(stringr::str_trim(unlist(x$uniprot_id)),
          collapse = "{}")

  uniprot_name <-
    paste(stringr::str_trim(unlist(x$uniprot_name)),
          collapse = "{}")

  pdb_ids <-
    paste(stringr::str_trim(unlist(x$pdb_ids)),
          collapse = "{}")


  hmdb_protein[[i]] <-
    data.frame(version = ifelse(is.null(unlist(x$version)), NA, unlist(x$version)),
               creation_date = ifelse(is.null(unlist(x$creation_date)), NA, unlist(x$creation_date)),
               update_date = ifelse(is.null(unlist(x$update_date)), NA, unlist(x$update_date)),
               accession = ifelse(is.null(unlist(x$accession)), NA, unlist(x$accession)),
               name,
               protein_type,
               synonyms,
               gene_name,
               general_function,
               specific_function,
               genbank_protein_id,
               uniprot_id,
               uniprot_name,
               pdb_ids)
}

hmdb_protein <-
  hmdb_protein %>%
  dplyr::bind_rows() %>%
  as.data.frame()

hmdb_protein <-
hmdb_protein %>%
  tibble::as_tibble()

save(hmdb_protein, file = "hmdb_protein", compress = "xz")
