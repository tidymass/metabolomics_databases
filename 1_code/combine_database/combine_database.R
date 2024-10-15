no_source()

setwd(masstools::get_project_wd())

####load data
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

foodb_ms1@spectra.info$From_food

###HMDB
load("2_data/HMDB/MS1/hmdb_ms1.rda")

hmdb_ms1

load("2_data/HMDB/MS2/hmdb_ms2.rda")

hmdb_ms2

sum(!is.na(hmdb_ms1@spectra.info$HMDB.ID))
sum(!is.na(hmdb_ms1@spectra.info$KEGG.ID))

###KEGG
load("2_data/KEGG/kegg_ms1.rda")

kegg_ms1

sum(!is.na(kegg_ms1@spectra.info$HMDB.ID))
sum(!is.na(kegg_ms1@spectra.info$KEGG.ID))

idx <-
  match(kegg_ms1@spectra.info$KEGG.ID, hmdb_ms1@spectra.info$KEGG.ID, incomparables = NA)

kegg_ms1@spectra.info$HMDB.ID <- hmdb_ms1@spectra.info$HMDB.ID[idx]

kegg_ms1@spectra.info$INCHI.ID <- NA
kegg_ms1@spectra.info$INCHI.ID <- hmdb_ms1@spectra.info$INCHI.ID[idx]

kegg_ms1@spectra.info$INCHIKEY.ID <- NA
kegg_ms1@spectra.info$INCHIKEY.ID <- hmdb_ms1@spectra.info$INCHIKEY.ID[idx]

kegg_ms1@spectra.info$SMILES.ID <- NA
kegg_ms1@spectra.info$SMILES.ID <- hmdb_ms1@spectra.info$SMILES.ID[idx]

##



###LIPIDBANKS
load("2_data/LIPIDBANK/lipidbank_ms1.rda")

lipidbank_ms1

sum(!is.na(lipidbank_ms1@spectra.info$HMDB.ID))
sum(!is.na(lipidbank_ms1@spectra.info$KEGG.ID))


###LIPIDMAPS
load("2_data/LIPIDMAPS/lipidmaps_ms1.rda")

lipidmaps_ms1

###T3DB
load("2_data/T3DB/t3db_ms1.rda")

t3db_ms1

###MASSBANK
load("2_data/MASSBANK/massbank_ms2.rda")

massbank_ms2

###MoNA
load("2_data/MONA/mona_ms2.rda")

mona_ms2

###METLIN
load("2_data/METLIN/metlin_ms2.rda")

metlin_ms2

###NIST
load("2_data/NIST/nist_ms2.rda")

nist_ms2

###SNYDER
load("2_data/SNYDER/MS2/mpsnyder_rplc_ms2.rda")
mpsnyder_rplc_ms2

load("2_data/SNYDER/MS2/mpsnyder_hilic_ms2.rda")
mpsnyder_hilic_ms2














