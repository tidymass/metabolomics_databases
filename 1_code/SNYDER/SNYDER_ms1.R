###no source
no_source()
setwd(masstools::get_project_wd())
rm(list = ls())
load("2_data/KEGG/kegg_ms1.rda")
load("2_data/HMDB/MS1/hmdb_ms1.rda")

setwd("2_data/SNYDER/MS1/")

library(tidyverse)
library(tidymass)

load("snyder_database_hilic0.0.3.rda")
load("snyder_database_rplc0.0.3.rda")

# data_hilic <-
#   snyder_database_hilic0.0.3@spectra.info
#
# data_rplc <-
#   snyder_database_rplc0.0.3@spectra.info
#
# colnames(data_hilic)
# colnames(data_rplc)
#
# data <-
#   rbind(data_hilic,
#         data_rplc)
#
# HMDB.ID <- data$HMDB.ID
# KEGG.ID <- data$KEGG.ID
#
# HMDB.ID[grep("C", data$HMDB.ID)] <-
#   data$KEGG.ID[grep("C", data$HMDB.ID)]
#
# KEGG.ID[grep("HMDB", data$KEGG.ID)] <-
#   data$HMDB.ID[grep("HMDB", data$KEGG.ID)]
#
# grep("C", HMDB.ID)
# grep("HMDB", KEGG.ID)
#
# data$HMDB.ID <- HMDB.ID
# data$KEGG.ID <- KEGG.ID
#
# data[which(data == "", arr.ind = TRUE)]
#
# data <-
#   data %>%
#   dplyr::mutate(synonyms = Compound.name)
#
# grep("\\([0-9]{1,2}\\)", data$Compound.name, value = TRUE)
#
# data$Compound.name <-
#   data$Compound.name %>%
#   stringr::str_replace("\\([0-9]{1,2}\\)", "")
#
# data$synonyms <-
#   data$synonyms %>%
#   stringr::str_split(" \\(") %>%
#   purrr::map(function(x) {
#     x <-
#       x %>%
#       stringr::str_replace("from", "") %>%
#       stringr::str_replace("\\?", "") %>%
#       stringr::str_trim()
#     if (length(x) == 1) {
#       return(x)
#     } else{
#       x[-1] <-
#         x[-1] %>%
#         stringr::str_replace("\\)$", "")
#       paste(x, collapse = "{}")
#     }
#   }) %>%
#   unlist()
#
# data$Compound.name <-
#   data$synonyms %>%
#   stringr::str_split(pattern = "\\{\\}") %>%
#   purrr::map(function(x) {
#     x[1]
#   }) %>%
#   unlist()
#
# data$HMDB.ID <-
#   data$HMDB.ID %>%
#   stringr::str_replace("\\|$", "") %>%
#   stringr::str_replace("\\|", "{}")
#
# data$KEGG.ID <-
#   data$KEGG.ID %>%
#   stringr::str_replace("\\|$", "") %>%
#   stringr::str_replace("^\\|", "") %>%
#   stringr::str_replace("\\|\\|", "") %>%
#   stringr::str_replace("\\|", "{}")
#
# data[which(data ==  "0", arr.ind = TRUE)] <- NA
# data[which(data ==  "", arr.ind = TRUE)] <- NA
#
# data$HMDB.ID <-
#   data$HMDB.ID %>%
#   purrr::map(function(x) {
#     if (is.na(x)) {
#       return(x)
#     }
#
#     x <-
#       stringr::str_split(x, "\\{\\}")
#     x <- x[[1]]
#
#     x %>%
#       lapply(function(y) {
#         if (nchar(y) == 9) {
#           y = stringr::str_replace(y, "HMDB", "HMDB00")
#           y
#         }
#       }) %>%
#       unlist() %>%
#       paste(collapse = "{}")
#
#   }) %>%
#   unlist()
#
# data[which(data ==  "0", arr.ind = TRUE)] <- NA
# data[which(data ==  "", arr.ind = TRUE)] <- NA
#
# data$CAS.ID <-
#   data$CAS.ID %>%
#   purrr::map(function(x) {
#     if (is.na(x)) {
#       return(NA)
#     }
#     x %>%
#       stringr::str_replace("\\|\\|", "") %>%
#       stringr::str_replace("\\|", "\\{\\}")
#   }) %>%
#   unlist()
#
# database1 <-
#   trans_id_database()
#
# data <- readxl::read_xlsx("data_manual.xlsx")
#
# data <-
#   data %>%
#   dplyr::arrange(Compound.name)
#
# data$Compound.name <-
#   data$Compound.name %>%
#   stringr::str_to_sentence() %>%
#   stringr::str_replace("\\([0-9]{1}\\)$", "")
#
# library(plyr)
#
# data <-
#   data %>%
#   plyr::dlply(.variables = .(Compound.name)) %>%
#   purrr::map(function(x) {
#     cat(unique(x$Compound.name), " ")
#     if (nrow(x) == 1) {
#       return(x)
#     } else{
#       HMDB.ID <- unique(x$HMDB.ID)
#       HMDB.ID <- HMDB.ID[!is.na(HMDB.ID)]
#       if (length(HMDB.ID) == 0) {
#         return(x)
#       } else{
#         x$HMDB.ID <- HMDB.ID
#         return(x)
#       }
#     }
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# # HMDB.ID1 <-
# #   data$Compound.name %>%
# #   purrr::map(function(x) {
# #     cat(x, " ")
# #     masstools::trans_ID(query = x,
# #                         from = "Chemical Name",
# #                         to = "Human Metabolome Database")
# #   }) %>%
# #   do.call(rbind, .) %>%
# #   as.data.frame()
# #
# # HMDB.ID2 <-
# #   data$CAS.ID %>%
# #   purrr::map(function(x) {
# #     cat(x, " ")
# #     masstools::trans_ID(query = x,
# #                         from = "CAS",
# #                         to = "Human Metabolome Database")
# #   }) %>%
# #   do.call(rbind, .) %>%
# #   as.data.frame()
# #
# # HMDB.ID3 <-
# #   data$KEGG.ID %>%
# #   purrr::map(function(x) {
# #     cat(x, " ")
# #     masstools::trans_ID(query = x,
# #                         from = "KEGG",
# #                         to = "Human Metabolome Database")
# #   }) %>%
# #   do.call(rbind, .) %>%
# #   as.data.frame()
# #
# #
# # dim(HMDB.ID1)
# # dim(HMDB.ID2)
# # dim(HMDB.ID3)
# #
# # HMDB.ID <-
# #   cbind(HMDB.ID = data$HMDB.ID,
# #         HMDB.ID1 = HMDB.ID1$`Human Metabolome Database`,
# #         HMDB.ID2 = HMDB.ID2$`Human Metabolome Database`,
# #         HMDB.ID3 = HMDB.ID3$`Human Metabolome Database`)
# #
# # final_HMDB.ID <-
# # apply(HMDB.ID, 1, function(x){
# # x <- as.character(x)
# # if(!is.na(x[1])){
# # return(x[1])
# # }
# # x <- x[!is.na(x)]
# # if(length(x) == 0){
# #   return(NA)
# # }
# # return(x[1])
# # })
# #
# # data$HMDB.ID <- final_HMDB.ID
# #
# # openxlsx::write.xlsx(data,
# #                      file = "data.xlsx",
# #                      asTable = TRUE,
# #                      overwrite = TRUE)
#
#
# data <-
#   readxl::read_xlsx("data_manual2.xlsx") %>%
#   as.data.frame()
#
# data[994, ]
#
# data$mz <-
#   as.numeric(data$mz)
#
# # HMDB.ID2 <- rep(NA, nrow(data))
# # for (i in 1:nrow(data)) {
# #   cat(i, " ")
# #   result <-
# #     search_hmdb_database(
# #       name = data$Compound.name[i],
# #       mz = data$mz[i],
# #       hmdb_metabolite_database = hmdb_ms1,
# #       mz_error_cutoff = 1000,
# #       similarity_score_cutoff = 0.8
# #     )
# #
# #   if (nrow(result) == 0) {
# #     HMDB.ID2[i] <- NA
# #   } else{
# #     HMDB.ID2[i] <- result$HMDB.ID[1]
# #   }
# # }
# #
# # HMDB.ID <-
# #   data.frame(HMDB.ID1 = data$HMDB.ID,
# #              HMDB.ID2 = HMDB.ID2)
# #
# # hmdb_id <-
# #   apply(HMDB.ID, 1, function(x) {
# #     x <- as.character(x)
# #     x <- x[!is.na(x)]
# #     if (length(x) == 0) {
# #       return(NA)
# #     }
# #     x[1]
# #   })
# #
# #
# # data$HMDB.ID <- hmdb_id
# #
# # openxlsx::write.xlsx(data,
# #                      file = "data.xlsx",
# #                      asTable = TRUE,
# #                      overwrite = TRUE)
#
# data <-
#   readxl::read_xlsx("data_manual3.xlsx") %>%
#   as.data.frame()
# data[which(data == 0, arr.ind = TRUE)] <- NA
# data[which(data == "NA", arr.ind = TRUE)] <- NA
# idx <- which(is.na(data$HMDB.ID) & !is.na(data$SMILES.ID))
#
# match(data$SMILES.ID[idx], hmdb_ms1@spectra.info$SMILES.ID)
#
# hmdb_data <-
#   hmdb_ms1@spectra.info
#
# idx <-
#   match(data$HMDB.ID, hmdb_data$HMDB.ID)
#
# # Compound.name <-
# #   data.frame(
# #     Compound.name1 = data$Compound.name,
# #     Compound.name2 = hmdb_data$Compound.name[idx]
# #   )
#
# CAS.ID <-
#   data.frame(value1 = data$CAS.ID,
#              value2 = hmdb_data$CAS.ID[idx]) %>%
#   apply(1, function(x) {
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# data$CAS.ID <-
#   CAS.ID
#
# KEGG.ID <-
#   data.frame(value1 = data$KEGG.ID,
#              value2 = hmdb_data$KEGG.ID[idx]) %>%
#   apply(1, function(x) {
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# data$KEGG.ID <-
#   KEGG.ID
#
# SMILES.ID <-
#   data.frame(value1 = data$SMILES.ID,
#              value2 = hmdb_data$SMILES.ID[idx]) %>%
#   apply(1, function(x) {
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# data$SMILES.ID <-
#   SMILES.ID
#
# INCHIKEY.ID <-
#   data.frame(value1 = data$INCHIKEY.ID,
#              value2 = hmdb_data$INCHIKEY.ID[idx]) %>%
#   apply(1, function(x) {
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# data$INCHIKEY.ID <-
#   INCHIKEY.ID
#
# INCHI.ID <-
#   data.frame(value1 = data$INCHI.ID,
#              value2 = hmdb_data$INCHI.ID[idx]) %>%
#   apply(1, function(x) {
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     return(x[1])
#   })
#
# data$INCHI.ID <-
#   INCHI.ID
#
# synonyms <-
#   data.frame(value1 = data$synonyms,
#              value2 = hmdb_data$Synonyms[idx]) %>%
#   apply(1, function(x) {
#     x <- as.character(x)
#     x <- x[!is.na(x)]
#     if (length(x) == 0) {
#       return(NA)
#     }
#     stringr::str_split(x, "\\{\\}") %>%
#       unlist() %>%
#       unique() %>%
#       paste(collapse = "{}")
#   })
#
#
# data$synonyms <-
#   synonyms
#
# data <-
#   data %>%
#   dplyr::rename(Synonyms = synonyms)
#
# new_data <-
#   hmdb_data[idx,] %>%
#   dplyr::select(setdiff(colnames(hmdb_data), colnames(data))) %>%
#   dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions, ))
#
# data <-
#   cbind(data, new_data)
#
# data_rplc <-
#   data %>%
#   dplyr::filter(stringr::str_detect(data$Lab.ID, "RPLC"))
#
# data_hilic <-
#   data %>%
#   dplyr::filter(stringr::str_detect(data$Lab.ID, "HILIC"))
#
# openxlsx::write.xlsx(data_rplc, file = "data_rplc.xlsx", asTable = TRUE)
# openxlsx::write.xlsx(data_hilic, file = "data_hilic.xlsx", asTable = TRUE)

library(metid)

mpsnyder_rplc_ms1 <-
  construct_database(
    path = ".",
    version = "2022-04-23",
    metabolite.info.name = "data_rplc.xlsx",
    source = "Michael_Snyder_RPLC",
    link = "https://med.stanford.edu/snyderlab.html",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = TRUE,
    threads = 3
  )

mpsnyder_hilic_ms1 <-
  construct_database(
    path = ".",
    version = "2022-04-23",
    metabolite.info.name = "data_hilic.xlsx",
    source = "Michael_Snyder_HILIC",
    link = "https://med.stanford.edu/snyderlab.html",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = TRUE,
    threads = 3
  )

load(here::here("2_data/source_system/source_system.rda"))

source(here::here("1_code/3_utils.R"))

library(tidyverse)
library(tidyselect)
library(metid)

mpsnyder_rplc_ms1 <-
  update_metid_database_source_system(
    database = mpsnyder_rplc_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(mpsnyder_rplc_ms1, file = "mpsnyder_rplc_ms1.rda")


mpsnyder_hilic_ms1 <-
  update_metid_database_source_system(
    database = mpsnyder_hilic_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(mpsnyder_hilic_ms1, file = "mpsnyder_hilic_ms1.rda")

