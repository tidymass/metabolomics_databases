no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

load("other_files/HMDB/MS1/hmdb_ms1.rda")
load("other_files/KEGG/kegg_ms1.rda")

###metanetx rpair database
load("other_files/METANETX/reaction/metanetx_rpair_database")

setwd("other_files/metabolic_network/")

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  tibble::as_tibble() %>%
  dplyr::filter(from_compound_METANETX_ID != to_compound_METANETX_ID)

metanetx_rpair_database$from_compound_HMDB_ID <-
  metanetx_rpair_database$from_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

metanetx_rpair_database$to_compound_HMDB_ID <-
  metanetx_rpair_database$to_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

metanetx_rpair_database$from_compound_KEGG_ID <-
  metanetx_rpair_database$from_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

metanetx_rpair_database$to_compound_KEGG_ID <-
  metanetx_rpair_database$to_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

sum(is.na(metanetx_rpair_database$from_compound_HMDB_ID))
sum(is.na(metanetx_rpair_database$to_compound_HMDB_ID))

sum(is.na(metanetx_rpair_database$from_compound_KEGG_ID))
sum(is.na(metanetx_rpair_database$to_compound_KEGG_ID))

idx1 <-
  match(metanetx_rpair_database$from_compound_KEGG_ID,
        hmdb_ms1@spectra.info$KEGG.ID,
        incomparables = NA)

idx2 <-
  match(metanetx_rpair_database$from_compound_KEGG_ID,
        kegg_ms1@spectra.info$KEGG.ID,
        incomparables = NA)

from_compound_HMDB_ID <-
  cbind(
    metanetx_rpair_database$from_compound_HMDB_ID,
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
  match(metanetx_rpair_database$to_compound_KEGG_ID,
        hmdb_ms1@spectra.info$KEGG.ID,
        incomparables = NA)

idx2 <-
  match(metanetx_rpair_database$to_compound_KEGG_ID,
        kegg_ms1@spectra.info$KEGG.ID,
        incomparables = NA)

to_compound_HMDB_ID <-
  cbind(
    metanetx_rpair_database$to_compound_HMDB_ID,
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

sum(is.na(from_compound_HMDB_ID))
sum(is.na(to_compound_HMDB_ID))

metanetx_rpair_database$from_compound_HMDB_ID <-
  from_compound_HMDB_ID

metanetx_rpair_database$to_compound_HMDB_ID <-
  to_compound_HMDB_ID

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::filter((
    !is.na(from_compound_HMDB_ID) |
      !is.na(from_compound_KEGG_ID)
  ) &
    (!is.na(to_compound_HMDB_ID) |
       !is.na(to_compound_KEGG_ID)))

dim(metanetx_rpair_database)

rpair_id <-
  seq_len(nrow(metanetx_rpair_database)) %>%
  purrr::map(function(i) {
    from <-
      c(
        metanetx_rpair_database$from_compound_HMDB_ID[i],
        metanetx_rpair_database$from_compound_KEGG_ID[i]
      )
    from <- from[!is.na(from)][1]

    to <-
      c(
        metanetx_rpair_database$to_compound_HMDB_ID[i],
        metanetx_rpair_database$to_compound_KEGG_ID[i]
      )
    to <- to[!is.na(to)][1]

    c(from, to) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::select(rpair_id, dplyr::everything())

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::filter(from_compound_METANETX_ID != to_compound_METANETX_ID) %>%
  dplyr::distinct(from_compound_METANETX_ID,
                  to_compound_METANETX_ID,
                  .keep_all = TRUE) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE)

dim(metanetx_rpair_database)


####Change compound name
from_compound_name1 <-
  hmdb_ms1@spectra.info$Compound.name[match(metanetx_rpair_database$from_compound_HMDB_ID,
                                            hmdb_ms1@spectra.info$HMDB.ID,
                                            incomparables = NA)]

from_compound_name2 <-
  kegg_ms1@spectra.info$Compound.name[match(metanetx_rpair_database$from_compound_KEGG_ID,
                                            kegg_ms1@spectra.info$KEGG.ID,
                                            incomparables = NA)]

from_compound_name <-
  cbind(from_compound_name1, from_compound_name2) %>%
  apply(1, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(x[1])
  })

sum(is.na(from_compound_name))
which(is.na(from_compound_name))

metanetx_rpair_database$from_compound_name <-
  from_compound_name

to_compound_name1 <-
  hmdb_ms1@spectra.info$Compound.name[match(metanetx_rpair_database$to_compound_HMDB_ID,
                                            hmdb_ms1@spectra.info$HMDB.ID,
                                            incomparables = NA)]

to_compound_name2 <-
  kegg_ms1@spectra.info$Compound.name[match(metanetx_rpair_database$to_compound_KEGG_ID,
                                            kegg_ms1@spectra.info$KEGG.ID,
                                            incomparables = NA)]

to_compound_name <-
  cbind(to_compound_name1, to_compound_name2) %>%
  apply(1, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(x[1])
  })

sum(is.na(to_compound_name))
which(is.na(to_compound_name))

metanetx_rpair_database$to_compound_name <-
  to_compound_name

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::filter(!is.na(from_compound_name) &
                  !is.na(to_compound_name))

from_compound_mz1 <-
  hmdb_ms1@spectra.info$mz[match(metanetx_rpair_database$from_compound_HMDB_ID,
                                 hmdb_ms1@spectra.info$HMDB.ID,
                                 incomparables = NA)]

from_compound_mz2 <-
  kegg_ms1@spectra.info$mz[match(metanetx_rpair_database$from_compound_KEGG_ID,
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

metanetx_rpair_database$from_compound_mz <-
  from_compound_mz

to_compound_mz1 <-
  hmdb_ms1@spectra.info$mz[match(metanetx_rpair_database$to_compound_HMDB_ID,
                                 hmdb_ms1@spectra.info$HMDB.ID,
                                 incomparables = NA)]

to_compound_mz2 <-
  kegg_ms1@spectra.info$mz[match(metanetx_rpair_database$to_compound_KEGG_ID,
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

metanetx_rpair_database$to_compound_mz <-
  to_compound_mz

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::filter(!is.na(from_compound_mz) &
                  !is.na(to_compound_mz))

dim(metanetx_rpair_database)

length(unique(
  c(
    metanetx_rpair_database$from_compound_METANETX_ID,
    metanetx_rpair_database$to_compound_METANETX_ID
  )
))

unique(
  c(
    metanetx_rpair_database$from_compound_name,
    metanetx_rpair_database$to_compound_name
  )
) %>%
  sort()


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

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::filter(!from_compound_name %in% remove_compund &
                  !to_compound_name %in% remove_compund)



library(clusterProfiler)
library(org.Hs.eg.db)

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::rename(reaction_Enzyme_EC_number = reaction_EC_number)

reaction_Enzyme_UNIPROT_ID <-
  seq_len(nrow(metanetx_rpair_database)) %>%
  purrr::map(function(i){
    # cat(i, " ")
    if(is.na(metanetx_rpair_database$reaction_Enzyme_EC_number[i])){
      return(NA)
    }

    reaction_Enzyme_EC_number <-
      metanetx_rpair_database$reaction_Enzyme_EC_number[i] %>%
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

metanetx_rpair_database$reaction_Enzyme_UNIPROT_ID <-
  reaction_Enzyme_UNIPROT_ID

dim(metanetx_rpair_database)

save(metanetx_rpair_database, file = "metanetx_rpair_database")
