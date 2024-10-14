no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

setwd("other_files/HMDB/Reaction/")

library(tidyverse)
library(xml2)
library(stringr)
library(massdatabase)
library(tibble)

load("hmdb_reaction_universal_database")
load("../MS1/hmdb_ms1.rda")
load("../protein/hmdb_protein")

hmdb_reaction_database_human <-
  hmdb_reaction_universal_database %>%
  as_tibble() %>%
  dplyr::rename(reaction_equation_name = equation,
                reaction_KEGG_ID = External_Links) %>%
  dplyr::mutate(reaction_KEGG_ID = stringr::str_replace_all(reaction_KEGG_ID, "Kegg Reaction ID: ", ""))

hmdb_rpair_database_human <- vector(mode = "list",
                                    length = nrow(hmdb_reaction_database_human))

for (i in 1:nrow(hmdb_reaction_database_human)) {
  if ((i / 100) %in% seq(1, 10000, 1)) {
    cat(i, " ")
  }
  #
  #   if (!is.na(hmdb_reaction_database_human$reaction_KEGG_ID[i])) {
  #     hmdb_rpair_database_human[[i]] <- NULL
  #     next()
  #   }

  x <-
    stringr::str_split(hmdb_reaction_database_human$reaction_equation_name[i],
                       "\\=")[[1]] %>%
    stringr::str_split("\\+")

  x <-
    seq_along(x) %>%
    lapply(function(idx) {
      y <- x[[idx]]
      y <-
        y %>%
        stringr::str_trim() %>%
        stringr::str_replace_all("^A ", "") %>%
        stringr::str_replace_all("^a ", "") %>%
        stringr::str_replace_all("^An ", "") %>%
        stringr::str_replace_all("^an ", "")
      y <-
        hmdb_ms1@spectra.info[match(y, hmdb_ms1@spectra.info$Compound.name),
                              c("Compound.name", "mz", "KEGG.ID", "HMDB.ID")]
      if (idx == 1) {
        y <-
          y %>%
          dplyr::rename(
            from_compound_HMDB_ID = HMDB.ID,
            from_compound_KEGG_ID = KEGG.ID,
            from_compound_name = Compound.name,
            from_compound_mz = mz
          ) %>%
          dplyr::filter(!is.na(from_compound_HMDB_ID))
        # dplyr::filter(from_compound_mz > 18.010565)
      } else{
        y <-
          y %>%
          dplyr::rename(
            to_compound_HMDB_ID = HMDB.ID,
            to_compound_KEGG_ID = KEGG.ID,
            to_compound_name = Compound.name,
            to_compound_mz = mz
          ) %>%
          dplyr::filter(!is.na(to_compound_HMDB_ID))
        # dplyr::filter(to_compound_mz > 18.010565)
      }
    })

  if (any(nrow(x[[1]]) == 0 | nrow(x[[2]]) == 0)) {
    hmdb_rpair_database_human[[i]] <- NULL
    next()
  }

  hmdb_rpair_database_human[[i]] <-
    1:nrow(x[[1]]) %>%
    purrr::map(function(index1) {
      1:nrow(x[[2]]) %>%
        purrr::map(function(index2) {
          data.frame(x[[1]][index1, , drop = FALSE],
                     x[[2]][index2, , drop = FALSE],
                     reaction_KEGG_ID = hmdb_reaction_database_human$reaction_KEGG_ID[i])
        }) %>%
        dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows()

  hmdb_rpair_database_human[[i]]$reaction_Enzyme <-
    hmdb_reaction_database_human$Enzymes[i]
  hmdb_rpair_database_human[[i]]$reaction_Status <-
    hmdb_reaction_database_human$Status[i]
}

number <-
  hmdb_rpair_database_human %>%
  lapply(function(x) {
    if (is.null(x)) {
      return(0)
    } else{
      nrow(x)
    }
  }) %>%
  unlist()

plot(number)

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  dplyr::bind_rows()

dim(hmdb_rpair_database_human)

####remove some compounds
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
                  !to_compound_name %in% remove_compund) %>%
  tibble::as_tibble()

temp <-
  hmdb_rpair_database_human %>%
  # dplyr::filter((from_compound_mz > 150 & from_compound_mz < 200) | (to_compound_mz > 150 & to_compound_mz < 200)) %>%
  dplyr::select(from_compound_name, to_compound_name)

temp <-
  data.frame(Compound.name = unique(c(
    temp$from_compound_name, temp$to_compound_name
  ))) %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("Compound.name", "Formula", "mz")]) %>%
  dplyr::arrange(Formula)

temp[nchar(temp$Formula) < 7, ]

rownames(hmdb_rpair_database_human) <- NULL

####remove the duplicates
rpair_id <-
  seq_len(nrow(hmdb_rpair_database_human)) %>%
  purrr::map(function(i) {
    c(
      hmdb_rpair_database_human$from_compound_HMDB_ID[i],
      hmdb_rpair_database_human$to_compound_HMDB_ID[i]
    ) %>%
      sort() %>%
      paste(collapse = "{}")
  }) %>%
  unlist()

hmdb_rpair_database_human$rpair_id <-
  rpair_id

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  dplyr::distinct(rpair_id, .keep_all = TRUE)

index <-
  match(hmdb_rpair_database_human$reaction_Enzyme,
        hmdb_protein$name)

hmdb_rpair_database_human$reaction_Enzyme[which(is.na(index))]

hmdb_rpair_database_human$reaction_Enzyme_UNIPROT_ID <-
  hmdb_protein$uniprot_id[index]

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  dplyr::rename(reaction_Enzyme_name = reaction_Enzyme)

###Enzyme ID to EC number
library(clusterProfiler)
library(org.Hs.eg.db)

reaction_Enzyme_EC_number <-
  rep(NA, nrow(hmdb_rpair_database_human))

for (i in 3030:nrow(hmdb_rpair_database_human)) {
  if ((i / 10) %in% seq(1, 10000, 1)) {
    cat(i, " ")
  }

  if(is.na(hmdb_rpair_database_human$reaction_Enzyme_UNIPROT_ID[i])){
    reaction_Enzyme_EC_number[i] <- NA
    next()
    }

  reaction_Enzyme_EC_number[i] <-
    clusterProfiler::bitr(
      geneID = hmdb_rpair_database_human$reaction_Enzyme_UNIPROT_ID[i],
      fromType = "UNIPROT",
      toType = "ENZYME",
      OrgDb = org.Hs.eg.db
    ) %>%
    pull(ENZYME) %>%
    paste(collapse = "{}")
}

reaction_Enzyme_EC_number[reaction_Enzyme_EC_number == ""] <- NA

hmdb_rpair_database_human$reaction_Enzyme_EC_number <- NA

save(hmdb_rpair_database_human, file = "hmdb_rpair_database_human")
