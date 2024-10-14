no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

load("other_files/HMDB/MS1/hmdb_ms1.rda")
load("other_files/KEGG/kegg_ms1.rda")

###BIGG rpair database
load("other_files/BIGG/reaction/bigg_rpair_database_human")

setwd("other_files/metabolic_network/")

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  tibble::as_tibble()

bigg_rpair_database_human$from_compound_HMDB_ID <-
  bigg_rpair_database_human$from_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

bigg_rpair_database_human$to_compound_HMDB_ID <-
  bigg_rpair_database_human$to_compound_HMDB_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

bigg_rpair_database_human$from_compound_KEGG_ID <-
  bigg_rpair_database_human$from_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

bigg_rpair_database_human$to_compound_KEGG_ID <-
  bigg_rpair_database_human$to_compound_KEGG_ID %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()


sum(is.na(bigg_rpair_database_human$from_compound_HMDB_ID))
sum(is.na(bigg_rpair_database_human$to_compound_HMDB_ID))

####update HMDB.ID
idx1 <-
  match(
    bigg_rpair_database_human$from_compound_KEGG_ID,
    hmdb_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

idx2 <-
  match(
    bigg_rpair_database_human$from_compound_KEGG_ID,
    kegg_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

from_compound_HMDB_ID <-
  cbind(
    bigg_rpair_database_human$from_compound_HMDB_ID,
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
    bigg_rpair_database_human$to_compound_KEGG_ID,
    hmdb_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

idx2 <-
  match(
    bigg_rpair_database_human$to_compound_KEGG_ID,
    kegg_ms1@spectra.info$KEGG.ID,
    incomparables = NA
  )

to_compound_HMDB_ID <-
  cbind(
    bigg_rpair_database_human$to_compound_HMDB_ID,
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

bigg_rpair_database_human$from_compound_HMDB_ID <-
  from_compound_HMDB_ID

bigg_rpair_database_human$to_compound_HMDB_ID <-
  to_compound_HMDB_ID

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::filter((
    !is.na(from_compound_HMDB_ID) |
      !is.na(from_compound_KEGG_ID)
  ) &
    (!is.na(to_compound_HMDB_ID) |
       !is.na(to_compound_KEGG_ID)))


bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::filter((
    !is.na(from_compound_HMDB_ID) |
      !is.na(from_compound_KEGG_ID)
  ) &
    (!is.na(to_compound_HMDB_ID) |
       !is.na(to_compound_KEGG_ID)))

dim(bigg_rpair_database_human)

#########change compound_name
idx1 <- match(
  bigg_rpair_database_human$from_compound_HMDB_ID,
  hmdb_ms1@spectra.info$HMDB.ID,
  incomparables = NA
)

idx2 <- match(
  bigg_rpair_database_human$from_compound_KEGG_ID,
  kegg_ms1@spectra.info$KEGG.ID,
  incomparables = NA
)

from_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    bigg_rpair_database_human$from_compound_name
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
  bigg_rpair_database_human$to_compound_HMDB_ID,
  hmdb_ms1@spectra.info$HMDB.ID,
  incomparables = NA
)

idx2 <- match(
  bigg_rpair_database_human$to_compound_KEGG_ID,
  kegg_ms1@spectra.info$KEGG.ID,
  incomparables = NA
)

to_compound_name <-
  cbind(
    hmdb_ms1@spectra.info$Compound.name[idx1],
    kegg_ms1@spectra.info$Compound.name[idx2],
    bigg_rpair_database_human$to_compound_name
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

bigg_rpair_database_human$from_compound_name <- from_compound_name
bigg_rpair_database_human$to_compound_name <- to_compound_name

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::filter(from_compound_BIGG_ID != to_compound_BIGG_ID)

###rpair ID
rpair_id <-
  seq_len(nrow(bigg_rpair_database_human)) %>%
  purrr::map(function(i) {
    from <-
      c(
        bigg_rpair_database_human$from_compound_HMDB_ID[i],
        bigg_rpair_database_human$from_compound_KEGG_ID[i]
      )
    from <- from[!is.na(from)][1]

    to <-
      c(
        bigg_rpair_database_human$to_compound_HMDB_ID[i],
        bigg_rpair_database_human$to_compound_KEGG_ID[i]
      )
    to <- to[!is.na(to)][1]

    c(from, to) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::select(rpair_id, dplyr::everything())

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::filter(from_compound_BIGG_ID != to_compound_BIGG_ID) %>%
  dplyr::distinct(from_compound_BIGG_ID, to_compound_BIGG_ID, .keep_all = TRUE) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE)

bigg_rpair_database_human$from_compound_BIGG_ID[1:10]

dim(bigg_rpair_database_human)

length(unique(
  c(
    bigg_rpair_database_human$from_compound_KEGG_ID,
    bigg_rpair_database_human$to_compound_KEGG_ID
  )
))

unique(
  c(
    bigg_rpair_database_human$from_compound_name,
    bigg_rpair_database_human$to_compound_name
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

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::filter(!from_compound_name %in% remove_compund &
                  !to_compound_name %in% remove_compund)

dim(bigg_rpair_database_human)


library(clusterProfiler)
library(org.Hs.eg.db)

bigg_rpair_database_human <-
bigg_rpair_database_human %>%
  dplyr::rename(reaction_Enzyme_EC_number = reaction_Enzyme_EC_Number)

reaction_Enzyme_UNIPROT_ID <-
  seq_len(nrow(bigg_rpair_database_human)) %>%
  purrr::map(function(i){
    # cat(i, " ")
    if(is.na(bigg_rpair_database_human$reaction_Enzyme_EC_number[i])){
      return(NA)
    }

    reaction_Enzyme_EC_number <-
      bigg_rpair_database_human$reaction_Enzyme_EC_number[i] %>%
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

bigg_rpair_database_human$reaction_Enzyme_UNIPROT_ID <-
  reaction_Enzyme_UNIPROT_ID

save(bigg_rpair_database_human, file = "bigg_rpair_database_human")
