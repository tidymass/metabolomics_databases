# ###n
# #####this should be run after hmdb_compound_ms1 has been created
# no_source()

# setwd(masstools::get_project_wd())
# rm(list = ls())

# load("3_data_analysis/HMDB/MS1/hmdb_compound_ms1.rda")

# library(tidyverse)
# library(xml2)
# library(stringr)

# file <-
# dir("2_data/HMDB/MS2/hmdb_experimental_msms_spectra/")

# hmdb_compound_ms2 <-
#   file %>%
# purrr::map(function(x){
#   cat(x, " ")
#   data <-
# read_xml(file.path("2_data/HMDB/MS2/hmdb_experimental_msms_spectra/", x)) %>%
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


# # Create a directory for data analysis
# dir.create("3_data_analysis/HMDB/MS2", showWarnings = FALSE)
# # Set the working directory to the newly created directory
# setwd("3_data_analysis/HMDB/MS2")

# remove_idx <-
# hmdb_compound_ms2 %>%
#   lapply(function(x) {
#     nrow(x$ms2)
#   }) %>%
#   unlist() %>%
#   `==`(0) %>%
#   which()

# hmdb_compound_ms2[[16776]]$ms2

# hmdb_compound_ms2 <-
#   hmdb_compound_ms2[-16776]


# spectra_info <-
#   hmdb_compound_ms2 %>%
#   purrr::map(function(x){
#     x$ms1_info
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()

# save(spectra_info, file = "spectra_info")

# spectra_data <-
#   hmdb_compound_ms2 %>%
#   purrr::map(function(x){
#     x$ms2
#   })

# save(spectra_data, file = "spectra_data")
