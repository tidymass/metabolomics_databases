no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

load("2_data/HMDB/MS1/hmdb_ms1.rda")
load("2_data/KEGG/kegg_ms1.rda")

###HMDB rpair database
load("2_data/HMDB/Reaction/hmdb_rpair_database_human")

setwd("2_data/metabolic_network/")

load("kegg_rpair_database_human")


hmdb_rpair_database_human$reaction_Enzyme_EC_number[1:5]

rpair_id <-
  seq_len(nrow(hmdb_rpair_database_human)) %>%
  purrr::map(function(i) {
    from <-
      c(
        hmdb_rpair_database_human$from_compound_HMDB_ID[i],
        hmdb_rpair_database_human$from_compound_KEGG_ID[i]
      )
    from <- from[!is.na(from)][1]

    to <-
      c(
        hmdb_rpair_database_human$to_compound_HMDB_ID[i],
        hmdb_rpair_database_human$to_compound_KEGG_ID[i]
      )
    to <- to[!is.na(to)][1]

    c(from, to) %>%
      sort() %>%
      paste(collapse = "_")
  }) %>%
  unlist()

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  dplyr::mutate(rpair_id = rpair_id) %>%
  dplyr::select(rpair_id, dplyr::everything())

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  dplyr::filter(from_compound_HMDB_ID != to_compound_HMDB_ID) %>%
  dplyr::distinct(from_compound_HMDB_ID, to_compound_HMDB_ID, .keep_all = TRUE) %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE)

rownames(hmdb_rpair_database_human) <- NULL

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  tibble::as_tibble()

dim(hmdb_rpair_database_human)

length(unique(
  c(
    hmdb_rpair_database_human$from_compound_KEGG_ID,
    hmdb_rpair_database_human$to_compound_KEGG_ID
  )
))

dim(hmdb_rpair_database_human)

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

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  dplyr::filter(!from_compound_name %in% remove_compund &
                  !to_compound_name %in% remove_compund)

# temp <-
#   hmdb_rpair_database_human %>%
#   dplyr::select(from_compound_name, to_compound_name)
#
# temp <-
#   data.frame(Compound.name = unique(c(
#     temp$from_compound_name, temp$to_compound_name
#   ))) %>%
#   dplyr::left_join(hmdb_ms1@spectra.info[, c("Compound.name", "Formula", "mz")]) %>%
#   dplyr::arrange(Formula)
#
# temp[nchar(temp$Formula) < 7, ]


dim(hmdb_rpair_database_human)

save(hmdb_rpair_database_human, file = "hmdb_rpair_database_human")
