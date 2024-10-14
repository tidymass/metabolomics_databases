no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

load("other_files/HMDB/MS1/hmdb_ms1.rda")
load("other_files/KEGG/kegg_ms1.rda")

###KEGG rpair database
load("other_files/KEGG/reaction/kegg_rpair_database_human")

setwd("other_files/metabolic_network/")

idx1 <-
  match(
    kegg_rpair_database_human$from_compound_KEGG_ID,
    hmdb_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

idx2 <-
  match(
    kegg_rpair_database_human$from_compound_KEGG_ID,
    kegg_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

from_compound_HMDB_ID <-
  cbind(
    kegg_rpair_database_human$from_compound_HMDB_ID,
    hmdb_ms1@spectra.info$HMDB.ID[idx1],
    kegg_ms1@spectra.info$HMDB.ID[idx2]
  ) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    } else{
      x[1]
    }
  })

idx1 <-
  match(
    kegg_rpair_database_human$to_compound_KEGG_ID,
    hmdb_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

idx2 <-
  match(
    kegg_rpair_database_human$to_compound_KEGG_ID,
    kegg_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

to_compound_HMDB_ID <-
  cbind(
    kegg_rpair_database_human$to_compound_HMDB_ID,
    hmdb_ms1@spectra.info$HMDB.ID[idx1],
    kegg_ms1@spectra.info$HMDB.ID[idx2]
  ) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    } else{
      x[1]
    }
  })

kegg_rpair_database_human$from_compound_HMDB_ID <-
  from_compound_HMDB_ID

kegg_rpair_database_human$to_compound_HMDB_ID <-
  to_compound_HMDB_ID

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::filter((
    !is.na(from_compound_HMDB_ID) |
      !is.na(from_compound_KEGG_ID)
  ) &
    (!is.na(to_compound_HMDB_ID) |
       !is.na(to_compound_KEGG_ID)))

kegg_rpair_database_human$from_compound_HMDB_ID <-
  kegg_rpair_database_human$from_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

kegg_rpair_database_human$to_compound_HMDB_ID <-
  kegg_rpair_database_human$to_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

kegg_rpair_database_human$from_compound_KEGG_ID <-
  kegg_rpair_database_human$from_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

kegg_rpair_database_human$to_compound_KEGG_ID <-
  kegg_rpair_database_human$to_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

rpair_id <-
  seq_len(nrow(kegg_rpair_database_human)) %>%
  purrr::map(function(i) {
    from <-
      c(
        kegg_rpair_database_human$from_compound_HMDB_ID[i],
        kegg_rpair_database_human$from_compound_KEGG_ID[i]
      )
    from <- from[!is.na(from)][1]

    to <-
      c(
        kegg_rpair_database_human$to_compound_HMDB_ID[i],
        kegg_rpair_database_human$to_compound_KEGG_ID[i]
      )
    to <- to[!is.na(to)][1]

    c(from, to) %>%
      sort() %>%
      paste(collapse = "_")

  }) %>%
  unlist()

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::select(rpair_id, dplyr::everything())

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::filter(from_compound_KEGG_ID != to_compound_KEGG_ID) %>%
  dplyr::distinct(from_compound_KEGG_ID, to_compound_KEGG_ID, .keep_all = TRUE) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE)

dim(kegg_rpair_database_human)

length(unique(
  c(
    kegg_rpair_database_human$from_compound_KEGG_ID,
    kegg_rpair_database_human$to_compound_KEGG_ID
  )
))

###update compound name
#########change compound_name
idx1 <- match(
  kegg_rpair_database_human$from_compound_HMDB_ID,
  hmdb_ms1@spectra.info$HMDB.ID,
  incomparables = NA
)

idx2 <- match(
  kegg_rpair_database_human$from_compound_KEGG_ID,
  kegg_ms1@spectra.info$KEGG.ID,
  incomparables = NA
)

from_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    kegg_rpair_database_human$from_compound_name
  ) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    } else{
      x[1]
    }
  })


idx1 <- match(
  kegg_rpair_database_human$to_compound_HMDB_ID,
  hmdb_ms1@spectra.info$HMDB.ID,
  incomparables = NA
)

idx2 <- match(
  kegg_rpair_database_human$to_compound_KEGG_ID,
  kegg_ms1@spectra.info$KEGG.ID,
  incomparables = NA
)

to_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    kegg_rpair_database_human$to_compound_name
  ) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    } else{
      x[1]
    }
  })

sum(is.na(from_compound_name))
sum(is.na(to_compound_name))

kegg_rpair_database_human$from_compound_name <- from_compound_name
kegg_rpair_database_human$to_compound_name <- to_compound_name

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::filter(!is.na(from_compound_name) &
                  !is.na(to_compound_name))

dim(kegg_rpair_database_human)

kegg_rpair_database_human[kegg_rpair_database_human == ""] <- NA

library(clusterProfiler)
library(org.Hs.eg.db)

reaction_Enzyme_UNIPROT_ID <-
  seq_len(nrow(kegg_rpair_database_human)) %>%
  purrr::map(function(i){
    # cat(i, " ")
    if(is.na(kegg_rpair_database_human$reaction_Enzyme_EC_number[i])){
      return(NA)
    }

    reaction_Enzyme_EC_number <-
    kegg_rpair_database_human$reaction_Enzyme_EC_number[i] %>%
      stringr::str_split("\\{\\}") %>%
      `[[`(1) %>%
      unique()

    reaction_Enzyme_UNIPROT_ID <-
      tryCatch(clusterProfiler::bitr(geneID = reaction_Enzyme_EC_number,
                                     fromType = "ENZYME",
                                     toType = "UNIPROT",
                                     OrgDb = org.Hs.eg.db) %>%
                 pull(UNIPROT) %>%
                 unique(), error = function(e) NULL)

    if(is.null(reaction_Enzyme_UNIPROT_ID)){
      return(NA)
    }

    reaction_Enzyme_UNIPROT_ID <-
      reaction_Enzyme_UNIPROT_ID[!is.na(reaction_Enzyme_UNIPROT_ID)] %>%
      paste(collapse = "{}")
    reaction_Enzyme_UNIPROT_ID
  }) %>%
  unlist()

reaction_Enzyme_UNIPROT_ID[reaction_Enzyme_UNIPROT_ID == ""] <- NA

kegg_rpair_database_human$reaction_Enzyme_UNIPROT_ID <-
  reaction_Enzyme_UNIPROT_ID

save(kegg_rpair_database_human, file = "kegg_rpair_database_human")
