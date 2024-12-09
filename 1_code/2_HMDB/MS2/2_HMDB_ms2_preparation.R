###n
#####this should be run after hmdb_compound_ms1 has been created
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
load("3_data_analysis/HMDB/MS1/hmdb_compound_ms1.rda")

setwd("3_data_analysis/HMDB/MS2/")

library(tidyverse)
library(xml2)
library(stringr)

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
colnames(hmdb_compound_ms1@spectra.info)

spectra_info2 <-
spectra_info2 %>%
  dplyr::rename(CE = "collision_energy_voltage")

spectra_info2 <-
spectra_info2 %>%
  dplyr::left_join(hmdb_compound_ms1@spectra.info %>% dplyr::select(-Lab.ID), by = c("HMDB.ID"))

spectra_info2$mz

hmdb_compound_ms2 <- hmdb_compound_ms1

hmdb_compound_ms2@spectra.info <- spectra_info2

hmdb_compound_ms2@spectra.data$Spectra.positive <- spectra_data_pos
hmdb_compound_ms2@spectra.data$Spectra.negative <- spectra_data_neg

hmdb_compound_ms2

hmdb_compound_ms2@spectra.info %>%
  dplyr::count(HMDB.ID, Polarity)

idx <-
hmdb_compound_ms2@spectra.info %>%
dplyr::filter(HMDB.ID == "HMDB0000288") %>%
  pull(Lab.ID)

idx2 <-
which(names(hmdb_compound_ms2@spectra.data$Spectra.positive) %in% idx)

masstools::ms2_plot(spectrum1 = hmdb_compound_ms2@spectra.data$Spectra.positive[[idx2[1]]][[1]],
                    spectrum2 = hmdb_compound_ms2@spectra.data$Spectra.positive[[idx2[5]]][[1]])


head(hmdb_compound_ms2@spectra.info$mz)
head(hmdb_compound_ms2@spectra.info$monisotopic_molecular_weight)

hmdb_compound_ms2@spectra.data$Spectra.positive <-
  hmdb_compound_ms2@spectra.data$Spectra.positive %>%
  purrr::map(function(x){
    x %>%
      lapply(function(y){
        y$mz <- as.numeric(y$mz)
        y$intensity <- as.numeric(y$intensity)
        y
      })
  })

hmdb_compound_ms2@spectra.data$Spectra.negative <-
  hmdb_compound_ms2@spectra.data$Spectra.negative %>%
  purrr::map(function(x){
    x %>%
      lapply(function(y){
        y$mz <- as.numeric(y$mz)
        y$intensity <- as.numeric(y$intensity)
        y
      })
  })

save(hmdb_compound_ms2, file = "hmdb_compound_ms2.rda")



