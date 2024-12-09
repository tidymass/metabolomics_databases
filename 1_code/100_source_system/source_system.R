no_source()

setwd(masstools::get_project_wd())
library(metid)

####load data
####BIGG
load("2_data/BIGG/bigg_ms1.rda")
bigg_ms1

###blood exposome
load("2_data/BLOODEXPOSOME/bloodexposome_ms1.rda")
bloodexposome_ms1

###CHEBI
load("2_data/CHEBI/chebi_ms1.rda")
chebi_ms1

###DrugBnak
load("2_data/DRUGBANK/drugbank_ms1.rda")

drugbank_ms1

###FOODB
load("2_data/FOODB/foodb_ms1.rda")

foodb_ms1

foodb_ms1@spectra.info$From_food <- "Yes"
foodb_ms1@spectra.info$From_food

###GNPS
load("2_data/GNPS/gnps_ms2.rda")

gnps_ms2

###HMDB
load("2_data/HMDB/MS1/hmdb_ms1.rda")

hmdb_ms1

load("2_data/HMDB/MS2/hmdb_ms2.rda")

hmdb_ms2

###KEGG
load("2_data/KEGG/kegg_ms1.rda")

kegg_ms1

###LIPIDBANKS
load("2_data/LIPIDBANK/lipidbank_ms1.rda")

lipidbank_ms1

###LIPIDMAPS
load("2_data/LIPIDMAPS/lipidmaps_ms1.rda")

lipidmaps_ms1

###MASSBANK
load("2_data/MASSBANK/massbank_ms2.rda")

massbank_ms2

###METLIN
load("2_data/METLIN/metlin_ms2.rda")

metlin_ms2

###MONA
load("2_data/MONA/mona_ms2.rda")

mona_ms2

###NIST
load("2_data/NIST/nist_ms2.rda")

nist_ms2

###T3DB
load("2_data/T3DB/t3db_ms1.rda")

t3db_ms1

###SNYDER
load("2_data/SNYDER/MS2/mpsnyder_rplc_ms2.rda")
mpsnyder_rplc_ms2

load("2_data/SNYDER/MS2/mpsnyder_hilic_ms2.rda")
mpsnyder_hilic_ms2


grep("From_", colnames(bigg_ms1@spectra.info), value = TRUE)
grep("From_", colnames(bloodexposome_ms1@spectra.info), value = TRUE)
grep("From_", colnames(chebi_ms1@spectra.info), value = TRUE)
grep("From_", colnames(drugbank_ms1@spectra.info), value = TRUE)
grep("From_", colnames(foodb_ms1@spectra.info), value = TRUE)
grep("From_", colnames(gnps_ms2@spectra.info), value = TRUE)
grep("From_", colnames(hmdb_ms1@spectra.info), value = TRUE)
grep("From_", colnames(kegg_ms1@spectra.info), value = TRUE)
grep("From_", colnames(lipidbank_ms1@spectra.info), value = TRUE)
grep("From_", colnames(lipidmaps_ms1@spectra.info), value = TRUE)
grep("From_", colnames(massbank_ms2@spectra.info), value = TRUE)
grep("From_", colnames(metlin_ms2@spectra.info), value = TRUE)
grep("From_", colnames(mona_ms2@spectra.info), value = TRUE)
grep("From_", colnames(nist_ms2@spectra.info), value = TRUE)
grep("From_", colnames(mpsnyder_rplc_ms2@spectra.info), value = TRUE)
grep("From_", colnames(mpsnyder_hilic_ms2@spectra.info), value = TRUE)

bigg_source_system <-
  bigg_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

bloodexposome_source_system <-
  bloodexposome_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

chebi_source_system <-
  chebi_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

drugbank_source_system <-
  drugbank_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

foodb_source_system <-
  foodb_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

gnps_source_system <-
  gnps_ms2@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

hmdb_source_system <-
  hmdb_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

kegg_source_system <-
  kegg_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

lipidbank_source_system <-
  lipidbank_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

lipidmaps_source_system <-
  lipidmaps_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

massbank_source_system <-
  massbank_ms2@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

metlin_source_system <-
  metlin_ms2@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

mona_source_system <-
  mona_ms2@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

nist_source_system <-
  nist_ms2@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

mpsnyder_rplc_source_system <-
  mpsnyder_rplc_ms2@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

mpsnyder_hilic_source_system <-
  mpsnyder_hilic_ms2@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

t3db_source_system <-
  t3db_ms1@spectra.info %>%
  dplyr::select(
    c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    tidyselect::starts_with("From")
  ) %>%
  dplyr::filter(!is.na(CAS.ID) & !is.na(HMDB.ID) & !is.na(KEGG.ID))

bigg_source_system <-
  bigg_source_system %>%
  dplyr::rename(From_plant = From_plan)

source_system <-
  bigg_source_system %>%
  dplyr::full_join(bloodexposome_source_system,
                   by = intersect(c(
                     grep("ID", colnames(bigg_source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(bloodexposome_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()

source_system <-
  source_system %>%
  dplyr::full_join(chebi_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(chebi_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()

source_system <-
  source_system %>%
  dplyr::full_join(drugbank_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(drugbank_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(foodb_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(foodb_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(gnps_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(gnps_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(hmdb_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(hmdb_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()

source_system <-
  source_system %>%
  dplyr::full_join(kegg_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(kegg_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(lipidbank_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(lipidbank_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()

source_system <-
  source_system %>%
  dplyr::full_join(lipidmaps_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(lipidmaps_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(massbank_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(massbank_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(metlin_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(metlin_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(mona_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(mona_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(nist_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(nist_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(mpsnyder_rplc_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(mpsnyder_rplc_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


source_system <-
  source_system %>%
  dplyr::full_join(mpsnyder_hilic_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(mpsnyder_hilic_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()




source_system <-
  source_system %>%
  dplyr::full_join(t3db_source_system,
                   by = intersect(c(
                     grep("ID", colnames(source_system), value = TRUE)
                   ),
                   c(
                     grep("ID", colnames(t3db_source_system), value = TRUE)
                   )),
                   na_matches = "never") %>%
  dplyr::distinct(CAS.ID, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  merge_same_source()


dim(source_system)

setwd(masstools::get_project_wd())
setwd("2_data/source_system/")
save(source_system, file = "source_system.rda")
