no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

load("2_data/HMDB/MS1/hmdb_ms1.rda")
load("2_data/KEGG/kegg_ms1.rda")

###recon3 rpair database
load("2_data/RECON3/reaction/recon3_rpair_database")

setwd("2_data/metabolic_network/")

recon3_rpair_database <-
  recon3_rpair_database %>%
  tibble::as_tibble() %>%
  dplyr::filter(from_compound_RECON3_ID != to_compound_RECON3_ID)

recon3_rpair_database$from_compound_HMDB_ID <-
  recon3_rpair_database$from_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

recon3_rpair_database$to_compound_HMDB_ID <-
  recon3_rpair_database$to_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

recon3_rpair_database$from_compound_KEGG_ID <-
  recon3_rpair_database$from_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

recon3_rpair_database$to_compound_KEGG_ID <-
  recon3_rpair_database$to_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()


sum(is.na(recon3_rpair_database$from_compound_HMDB_ID))
sum(is.na(recon3_rpair_database$to_compound_HMDB_ID))

sum(is.na(recon3_rpair_database$from_compound_KEGG_ID))
sum(is.na(recon3_rpair_database$to_compound_KEGG_ID))

###Add HMDB ID
idx1 <-
  match(recon3_rpair_database$from_compound_KEGG_ID,
        hmdb_ms1@spectra.info$KEGG.ID,
        incomparables = NA)

idx2 <-
  match(
    recon3_rpair_database$from_compound_PUBCHEM_ID,
    hmdb_ms1@spectra.info$PUBCHEM.ID,
    incomparables = NA
  )

idx3 <-
  match(recon3_rpair_database$from_compound_CHEBI_ID,
        hmdb_ms1@spectra.info$CHEBI.ID,
        incomparables = NA)

idx4 <-
  match(recon3_rpair_database$from_compound_CAS_ID,
        hmdb_ms1@spectra.info$CAS.ID,
        incomparables = NA)

idx5 <-
  match(recon3_rpair_database$from_compound_INCHIKEY,
        hmdb_ms1@spectra.info$INCHIKEY.ID,
        incomparables = NA)

idx6 <-
  match(recon3_rpair_database$from_compound_SMILE,
        hmdb_ms1@spectra.info$SMILES.ID,
        incomparables = NA)

from_compound_HMDB_ID <-
  cbind(
    recon3_rpair_database$from_compound_HMDB_ID,
    hmdb_ms1@spectra.info$HMDB.ID[idx1],
    hmdb_ms1@spectra.info$HMDB.ID[idx2],
    hmdb_ms1@spectra.info$HMDB.ID[idx3],
    hmdb_ms1@spectra.info$HMDB.ID[idx4],
    hmdb_ms1@spectra.info$HMDB.ID[idx5],
    hmdb_ms1@spectra.info$HMDB.ID[idx6]
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
  match(recon3_rpair_database$to_compound_KEGG_ID,
        hmdb_ms1@spectra.info$KEGG.ID,
        incomparables = NA)

idx2 <-
  match(recon3_rpair_database$to_compound_PUBCHEM_ID,
        hmdb_ms1@spectra.info$PUBCHEM.ID,
        incomparables = NA)

idx3 <-
  match(recon3_rpair_database$to_compound_CHEBI_ID,
        hmdb_ms1@spectra.info$CHEBI.ID,
        incomparables = NA)

idx4 <-
  match(recon3_rpair_database$to_compound_CAS_ID,
        hmdb_ms1@spectra.info$CAS.ID,
        incomparables = NA)

idx5 <-
  match(recon3_rpair_database$to_compound_INCHIKEY,
        hmdb_ms1@spectra.info$INCHIKEY.ID,
        incomparables = NA)

idx6 <-
  match(recon3_rpair_database$to_compound_SMILE,
        hmdb_ms1@spectra.info$SMILES.ID,
        incomparables = NA)

to_compound_HMDB_ID <-
  cbind(
    recon3_rpair_database$to_compound_HMDB_ID,
    hmdb_ms1@spectra.info$HMDB.ID[idx1],
    hmdb_ms1@spectra.info$HMDB.ID[idx2],
    hmdb_ms1@spectra.info$HMDB.ID[idx3],
    hmdb_ms1@spectra.info$HMDB.ID[idx4],
    hmdb_ms1@spectra.info$HMDB.ID[idx5],
    hmdb_ms1@spectra.info$HMDB.ID[idx6]
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

sum(is.na(from_compound_HMDB_ID))
sum(is.na(to_compound_HMDB_ID))

recon3_rpair_database$from_compound_HMDB_ID <-
  from_compound_HMDB_ID

recon3_rpair_database$to_compound_HMDB_ID <-
  to_compound_HMDB_ID

recon3_rpair_database <-
  recon3_rpair_database %>%
  dplyr::filter((
    !is.na(from_compound_HMDB_ID) |
      !is.na(from_compound_KEGG_ID)
  ) &
    (!is.na(to_compound_HMDB_ID) |
       !is.na(to_compound_KEGG_ID)))


dim(recon3_rpair_database)

rpair_id <-
  seq_len(nrow(recon3_rpair_database)) %>%
  purrr::map(function(i) {
    from <-
      c(
        recon3_rpair_database$from_compound_HMDB_ID[i],
        recon3_rpair_database$from_compound_KEGG_ID[i]
      )
    from <- from[!is.na(from)][1]

    to <-
      c(
        recon3_rpair_database$to_compound_HMDB_ID[i],
        recon3_rpair_database$to_compound_KEGG_ID[i]
      )
    to <- to[!is.na(to)][1]

    c(from, to) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

recon3_rpair_database <-
  recon3_rpair_database %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::select(rpair_id, dplyr::everything())

recon3_rpair_database <-
  recon3_rpair_database %>%
  dplyr::filter(from_compound_RECON3_ID != to_compound_RECON3_ID) %>%
  dplyr::distinct(from_compound_RECON3_ID,
                  to_compound_RECON3_ID,
                  .keep_all = TRUE) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE)

dim(recon3_rpair_database)

length(unique(
  c(
    recon3_rpair_database$from_compound_RECON3_ID,
    recon3_rpair_database$to_compound_RECON3_ID
  )
))

###update compound_name
###add compound_name
#########change compound_name
idx1 <- match(recon3_rpair_database$from_compound_HMDB_ID,
              hmdb_ms1@spectra.info$HMDB.ID,
              incomparables = NA)

idx2 <- match(recon3_rpair_database$from_compound_KEGG_ID,
              kegg_ms1@spectra.info$KEGG.ID,
              incomparables = NA)

from_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    recon3_rpair_database$from_compound_name
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


idx1 <- match(recon3_rpair_database$to_compound_HMDB_ID,
              hmdb_ms1@spectra.info$HMDB.ID,
              incomparables = NA)

idx2 <- match(recon3_rpair_database$to_compound_KEGG_ID,
              kegg_ms1@spectra.info$KEGG.ID,
              incomparables = NA)

to_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    recon3_rpair_database$to_compound_name
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

recon3_rpair_database$from_compound_name <- from_compound_name
recon3_rpair_database$to_compound_name <- to_compound_name

from_compound_mz1 <-
  hmdb_ms1@spectra.info$mz[match(recon3_rpair_database$from_compound_HMDB_ID,
                                 hmdb_ms1@spectra.info$HMDB.ID,
                                 incomparables = NA)]

from_compound_mz2 <-
  kegg_ms1@spectra.info$mz[match(recon3_rpair_database$from_compound_KEGG_ID,
                                 kegg_ms1@spectra.info$KEGG.ID,
                                 incomparables = NA)]

from_compound_mz <-
  cbind(from_compound_mz1, from_compound_mz2) %>%
  apply(1, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(as.numeric(x[1]))
  })

sum(is.na(from_compound_mz))

recon3_rpair_database$from_compound_mz <-
  from_compound_mz

to_compound_mz1 <-
  hmdb_ms1@spectra.info$mz[match(recon3_rpair_database$to_compound_HMDB_ID,
                                 hmdb_ms1@spectra.info$HMDB.ID,
                                 incomparables = NA)]

to_compound_mz2 <-
  kegg_ms1@spectra.info$mz[match(recon3_rpair_database$to_compound_KEGG_ID,
                                 kegg_ms1@spectra.info$KEGG.ID,
                                 incomparables = NA)]

to_compound_mz <-
  cbind(to_compound_mz1, 2) %>%
  apply(1, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(as.numeric(x[1]))
  })

sum(is.na(to_compound_mz))
which(is.na(to_compound_mz))

recon3_rpair_database$to_compound_mz <-
  to_compound_mz

recon3_rpair_database <-
  recon3_rpair_database %>%
  dplyr::filter(!is.na(from_compound_mz) &
                  !is.na(to_compound_mz))

remove_compund <-
  c(
    "Hydrogen Ion",
    "Hydrogen",
    "Helium",
    "Iron",
    "Oxygen",
    "Phosphate",
    "I(-)",
    "Pyrophosphate",
    "Superoxide",
    "Nitric oxide",
    "Calcium",
    "Lithium",
    "Beryllium",
    "Water",
    "Nitrite",
    "NH4OH",
    "Carbon dioxide",
    "Hydrogen sulfide",
    "Fluoride",
    "Thiocyanate",
    "Sulfide",
    "Hydrogen sulfide",
    "Cyanide",
    "Potassium",
    "Ammonium",
    "Chloride ion",
    "Cotinine methonium ion",
    "Glutathione episulfonium ion",
    "Calcium ionophore",
    "Magnesium ionophore IV",
    "Ammonia",
    "Mercury",
    "Magnesium",
    "trdox",
    "UWM6",
    "Silver",
    "Zinc",
    "Cobalt",
    "Nitrate",
    "Cadmium",
    "Copper",
    "Fe2+",
    "Sodium",
    "Lipid A",
    "Nitrogen",
    "Hydrogen peroxide",
    "Hydrogen carbonate",
    "Carbon monoxide",
    "Hydrogen cyanide",
    "Hydrogen selenide",
    "Hydrogen cyanide",
    "Bromide",
    "Formaldehyde",
    "Hydrochloric acid",
    "Arsenite",
    "Iodine"
  )

recon3_rpair_database <-
  recon3_rpair_database %>%
  dplyr::filter(!from_compound_name %in% remove_compund &
                  !to_compound_name %in% remove_compund)

dim(recon3_rpair_database)


library(clusterProfiler)
library(org.Hs.eg.db)

reaction_Enzyme_UNIPROT_ID <-
  seq_len(nrow(recon3_rpair_database)) %>%
  purrr::map(function(i){
    # cat(i, " ")
    if(is.na(recon3_rpair_database$reaction_Enzyme_EC_number[i])){
      return(NA)
    }

    reaction_Enzyme_EC_number <-
      recon3_rpair_database$reaction_Enzyme_EC_number[i] %>%
      stringr::str_split("\\{\\}") %>%
      `[[`(1) %>%
      stringr::str_replace_all("EC\\:", "") %>%
      stringr::str_trim() %>%
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

recon3_rpair_database$reaction_Enzyme_UNIPROT_ID <-
  reaction_Enzyme_UNIPROT_ID

save(recon3_rpair_database, file = "recon3_rpair_database")
