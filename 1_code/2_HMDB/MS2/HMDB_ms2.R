###n
#####this should be run after hmdb_ms1 has been created
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
load("2_data/HMDB/MS1/hmdb_ms1.rda")
setwd("2_data/HMDB/MS2/")

library(tidyverse)
library(xml2)
library(stringr)

# file <-
# dir("hmdb_experimental_msms_spectra/")
#
# hmdb_ms2 <-
#   file %>%
# purrr::map(function(x){
#   cat(x, " ")
#   data <-
# read_xml(file.path("hmdb_experimental_msms_spectra/", x)) %>%
#     xml2::as_list()
#   Instrument_type <-
#   unlist(data$`ms-ms`$`instrument-type`)
#   Instrument_type <-
#     ifelse(is.null(Instrument_type), NA, Instrument_type)
#   Polarity <-
#     unlist(data$`ms-ms`$`ionization-mode`)
#   Polarity <-
#     ifelse(is.null(Polarity), NA, Polarity)
#   collision_energy_level <-
#   unlist(data$`ms-ms`$`collision-energy-level`)
#   collision_energy_level <-
#     ifelse(is.null(collision_energy_level), NA, collision_energy_level)
#   collision_energy_voltage <-
#     unlist(data$`ms-ms`$`collision-energy-voltage`)
#   collision_energy_voltage <-
#     ifelse(is.null(collision_energy_voltage), NA, collision_energy_voltage)
#
#   chromatography_type <-
#     unlist(data$`ms-ms`$`chromatography-type`)
#   chromatography_type <-
#     ifelse(is.null(chromatography_type), NA, chromatography_type)
#   analyzer_type <-
#     unlist(data$`ms-ms`$`analyzer-type`)
#   analyzer_type <-
#     ifelse(is.null(analyzer_type), NA, analyzer_type)
#   ionization_type <-
#     unlist(data$`ms-ms`$`ionization-type`)
#   ionization_type <-
#     ifelse(is.null(ionization_type), NA, ionization_type)
#   charge_type <-
#     unlist(data$`ms-ms`$`charge-type`)
#   charge_type <-
#     ifelse(is.null(charge_type), NA, charge_type)
#   adduct <-
#     unlist(data$`ms-ms`$adduct)
#   adduct <-
#     ifelse(is.null(adduct), NA, adduct)
#   adduct_type <-
#     unlist(data$`ms-ms`$`adduct-type`)
#   adduct_type <-
#     ifelse(is.null(adduct_type), NA, adduct_type)
#   adduct_mass <-
#     unlist(data$`ms-ms`$`adduct-mass`)
#   adduct_mass <-
#     ifelse(is.null(adduct_mass), NA, adduct_mass)
#   ms1_info <-
#     data.frame(HMDB.ID = stringr::str_extract(x, "HMDB[0-9]{7,9}"),
#                Instrument_type = Instrument_type,
#                Polarity = Polarity, collision_energy_level = collision_energy_level,
#                collision_energy_voltage = collision_energy_voltage,
#                chromatography_type = chromatography_type,
#                analyzer_type = analyzer_type,
#                ionization_type = ionization_type,
#                charge_type = charge_type,
#                adduct = adduct,
#                adduct_type = adduct_type,
#                adduct_mass = adduct_mass)
#
#   ms2 <-
#     data$`ms-ms`$`ms-ms-peaks`
#   if(is.null(ms2)){
#     ms2 <- data.frame()
#   }else{
#     ms2 <-
#       lapply(ms2, function(y){
#         unlist(y)
#       }) %>%
#       dplyr::bind_rows() %>%
#       as.data.frame() %>%
#       dplyr::select(`mass-charge`, intensity) %>%
#       dplyr::rename(mz = `mass-charge`)
#   }
#   list(ms1_info = ms1_info,
#        ms2 = ms2)
# })
#
# remove_idx <-
# hmdb_ms2 %>%
#   lapply(function(x) {
#     nrow(x$ms2)
#   }) %>%
#   unlist() %>%
#   `==`(0) %>%
#   which()
#
# hmdb_ms2[[16776]]$ms2
#
# hmdb_ms2 <-
#   hmdb_ms2[-16776]
#
#
# spectra_info <-
#   hmdb_ms2 %>%
#   purrr::map(function(x){
#     x$ms1_info
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# save(spectra_info, file = "spectra_info")
#
# spectra_data <-
#   hmdb_ms2 %>%
#   purrr::map(function(x){
#     x$ms2
#   })
#
# save(spectra_data, file = "spectra_data")

load("spectra_info")
load("spectra_data")

spectra_info$HMDB.ID

spectra_info[which(spectra_info == "NA", arr.ind = TRUE)] <- NA
spectra_info[which(spectra_info == "n/a", arr.ind = TRUE)] <- NA
spectra_info[which(spectra_info == "N/A", arr.ind = TRUE)] <- NA

spectra_info <-
spectra_info %>%
  dplyr::select(HMDB.ID, Instrument_type, Polarity, collision_energy_voltage, adduct)

spectra_info$Polarity

remove_idx <-
  which(is.na(spectra_info$Polarity))
remove_idx
dim(spectra_info)
length(spectra_data)

spectra_info <-
  spectra_info[-remove_idx,]

spectra_data <-
  spectra_data[-remove_idx]

spectra_info <-
  spectra_info %>%
  dplyr::mutate(Polarity = case_when(
    Polarity == "positive" ~ "Positive",
    Polarity == "negative" ~ "Negative",
    TRUE ~ Polarity
  ))

library(plyr)

spectra_info$Lab.ID <-
  masstools::name_duplicated(spectra_info$HMDB.ID) %>%
  paste("shen", sep = "_")

spectra_info2 <-
  spectra_info %>%
  plyr::dlply(.variables = .(HMDB.ID)) %>%
  purrr::map(function(y) {
    if (sum(is.na(y$collision_energy_voltage)) > 0) {
      y$collision_energy_voltage[is.na(y$collision_energy_voltage)] <-
        paste("Unknown", 1:length(y$collision_energy_voltage[is.na(y$collision_energy_voltage)]), sep = "_")
    }
    y
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

spectra_info2 <-
  spectra_info2[match(spectra_info$Lab.ID, spectra_info2$Lab.ID),]

sum(spectra_info2$Lab.ID == spectra_info$Lab.ID)

spectra_data2 <-
1:length(spectra_data) %>%
  purrr::map(function(i){
    x <- spectra_data[[i]]
    x <- list(x)
    names(x) <-
      spectra_info2$collision_energy_voltage[i]
    x
  })

names(spectra_data2) <- spectra_info2$Lab.ID

######positive mode
spectra_info2$Lab.ID == names(spectra_data2)

index_pos <- which(spectra_info2$Polarity == "Positive")
index_neg <- which(spectra_info2$Polarity == "Negative")

spectra_info_pos <- spectra_info2[index_pos,]
spectra_data_pos <- spectra_data2[index_pos]

spectra_info_neg <- spectra_info2[index_neg,]
spectra_data_neg <- spectra_data2[index_neg]

colnames(spectra_info2)
colnames(hmdb_ms1@spectra.info)

spectra_info2 <-
spectra_info2 %>%
  dplyr::rename(CE = "collision_energy_voltage")

spectra_info2 <-
spectra_info2 %>%
  dplyr::left_join(hmdb_ms1@spectra.info %>% dplyr::select(-Lab.ID), by = c("HMDB.ID"))

spectra_info2$mz

hmdb_ms2 <- hmdb_ms1

hmdb_ms2@spectra.info <- spectra_info2

hmdb_ms2@spectra.data$Spectra.positive <- spectra_data_pos
hmdb_ms2@spectra.data$Spectra.negative <- spectra_data_neg

hmdb_ms2

hmdb_ms2@spectra.info %>%
  dplyr::count(HMDB.ID, Polarity)

idx <-
hmdb_ms2@spectra.info %>%
dplyr::filter(HMDB.ID == "HMDB0000288") %>%
  pull(Lab.ID)

idx2 <-
which(names(hmdb_ms2@spectra.data$Spectra.positive) %in% idx)

masstools::ms2_plot(spectrum1 = hmdb_ms2@spectra.data$Spectra.positive[[idx2[1]]][[1]],
                    spectrum2 = hmdb_ms2@spectra.data$Spectra.positive[[idx2[5]]][[1]])


head(hmdb_ms2@spectra.info$mz)
head(hmdb_ms2@spectra.info$monisotopic_molecular_weight)

hmdb_ms2@spectra.data$Spectra.positive <-
  hmdb_ms2@spectra.data$Spectra.positive %>%
  purrr::map(function(x){
    x %>%
      lapply(function(y){
        y$mz <- as.numeric(y$mz)
        y$intensity <- as.numeric(y$intensity)
        y
      })
  })

hmdb_ms2@spectra.data$Spectra.negative <-
  hmdb_ms2@spectra.data$Spectra.negative %>%
  purrr::map(function(x){
    x %>%
      lapply(function(y){
        y$mz <- as.numeric(y$mz)
        y$intensity <- as.numeric(y$intensity)
        y
      })
  })

save(hmdb_ms2, file = "hmdb_ms2.rda")



