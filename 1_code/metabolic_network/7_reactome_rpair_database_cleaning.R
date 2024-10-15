# no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

load("2_data/HMDB/MS1/hmdb_ms1.rda")
load("2_data/KEGG/kegg_ms1.rda")

###reactome rpair database
load("2_data/REACTOME/reaction/reactome_rpair_database_human")

setwd("2_data/metabolic_network/")

reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  tibble::as_tibble() %>%
  dplyr::filter(from_compound_REACTOME_ID != to_compound_REACTOME_ID)

reactome_rpair_database_human$from_compound_HMDB_ID <-
  reactome_rpair_database_human$from_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

reactome_rpair_database_human$to_compound_HMDB_ID <-
  reactome_rpair_database_human$to_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

reactome_rpair_database_human$from_compound_KEGG_ID <-
  reactome_rpair_database_human$from_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

reactome_rpair_database_human$to_compound_KEGG_ID <-
  reactome_rpair_database_human$to_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

sum(is.na(reactome_rpair_database_human$from_compound_HMDB_ID))
sum(is.na(reactome_rpair_database_human$to_compound_HMDB_ID))

sum(is.na(reactome_rpair_database_human$from_compound_KEGG_ID))
sum(is.na(reactome_rpair_database_human$to_compound_KEGG_ID))

dim(reactome_rpair_database_human)

###Add HMDB ID
idx1 <-
  match(
    reactome_rpair_database_human$from_compound_KEGG_ID,
    hmdb_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

idx2 <-
  match(
    stringr::str_replace(
      reactome_rpair_database_human$from_compound_CHEBI_ID,
      "CHEBI\\:",
      ""
    ),
    hmdb_ms1@spectra.info$CHEBI.ID,
    incomparables = NA
  )

from_compound_HMDB_ID <-
  cbind(
    reactome_rpair_database_human$from_compound_HMDB_ID,
    hmdb_ms1@spectra.info$HMDB.ID[idx1],
    hmdb_ms1@spectra.info$HMDB.ID[idx2]
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
    reactome_rpair_database_human$to_compound_KEGG_ID,
    hmdb_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

idx2 <-
  match(
    stringr::str_replace(
      reactome_rpair_database_human$to_compound_CHEBI_ID,
      "CHEBI\\:",
      ""
    ),
    hmdb_ms1@spectra.info$CHEBI.ID,
    incomparables = NA
  )

to_compound_HMDB_ID <-
  cbind(
    reactome_rpair_database_human$to_compound_HMDB_ID,
    hmdb_ms1@spectra.info$HMDB.ID[idx1],
    hmdb_ms1@spectra.info$HMDB.ID[idx2]
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

reactome_rpair_database_human$from_compound_HMDB_ID <-
  from_compound_HMDB_ID

reactome_rpair_database_human$to_compound_HMDB_ID <-
  to_compound_HMDB_ID

reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  dplyr::filter((
    !is.na(from_compound_HMDB_ID) |
      !is.na(from_compound_KEGG_ID)
  ) &
    (!is.na(to_compound_HMDB_ID) |
       !is.na(to_compound_KEGG_ID)))

dim(reactome_rpair_database_human)


rpair_id <-
  seq_len(nrow(reactome_rpair_database_human)) %>%
  purrr::map(function(i) {
    from <-
      c(
        reactome_rpair_database_human$from_compound_HMDB_ID[i],
        reactome_rpair_database_human$from_compound_KEGG_ID[i]
      )
    from <- from[!is.na(from)][1]

    to <-
      c(
        reactome_rpair_database_human$to_compound_HMDB_ID[i],
        reactome_rpair_database_human$to_compound_KEGG_ID[i]
      )
    to <- to[!is.na(to)][1]

    c(from, to) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::select(rpair_id, dplyr::everything())

reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  dplyr::filter(from_compound_REACTOME_ID != to_compound_REACTOME_ID) %>%
  dplyr::distinct(from_compound_REACTOME_ID,
                  to_compound_REACTOME_ID,
                  .keep_all = TRUE) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE)

dim(reactome_rpair_database_human)

length(unique(
  c(
    reactome_rpair_database_human$from_compound_REACTOME_ID,
    reactome_rpair_database_human$to_compound_REACTOME_ID
  )
))


###update compound_name
###add compound_name
#########change compound_name
idx1 <- match(
  reactome_rpair_database_human$from_compound_HMDB_ID,
  hmdb_ms1@spectra.info$HMDB.ID,
  incomparables = NA
)

idx2 <- match(
  reactome_rpair_database_human$from_compound_KEGG_ID,
  kegg_ms1@spectra.info$KEGG.ID,
  incomparables = NA
)

from_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    reactome_rpair_database_human$from_compound_name
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
  reactome_rpair_database_human$to_compound_HMDB_ID,
  hmdb_ms1@spectra.info$HMDB.ID,
  incomparables = NA
)

idx2 <- match(
  reactome_rpair_database_human$to_compound_KEGG_ID,
  kegg_ms1@spectra.info$KEGG.ID,
  incomparables = NA
)

to_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    reactome_rpair_database_human$to_compound_name
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

reactome_rpair_database_human$from_compound_name <-
  from_compound_name
reactome_rpair_database_human$to_compound_name <- to_compound_name

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

reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  dplyr::filter(!from_compound_name %in% remove_compund &
                  !to_compound_name %in% remove_compund)

dim(reactome_rpair_database_human)

save(reactome_rpair_database_human, file = "reactome_rpair_database_human")
