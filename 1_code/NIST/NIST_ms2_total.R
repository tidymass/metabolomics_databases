###source
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
load("other_files/HMDB/MS1/hmdb_ms1.rda")
load("other_files/KEGG/kegg_ms1.rda")
source("R/read_msp_data.R")
source("R/3_utils.R")
source("R/12_NIST.R")
setwd("other_files/NIST")

library(tidyverse)
library(tidymass)

###
data1 <-
  read_msp_data(file = "190331_NIST17_MONA_MSMS_NEG.msp",
                source = "nist",
                threads = 5)

length(data1)

data2 <-
  read_msp_data(file = "MSMS-Pos-VS14.msp",
                source = "nist",
                threads = 5)

length(data1)
length(data2)

data <-
  c(data1, data2)

length(data)

nist_total_ms2_neg <-
  convert_nist2metid(data = data1,
                     path = "NIST_total_neg",
                     threads = 5)

nist_total_ms2_pos <-
  convert_nist2metid(data = data2,
                     path = "NIST_total_pos",
                     threads = 5)

save(nist_total_ms2_pos, file = "nist_total_ms2_pos.rda")
save(nist_total_ms2_neg, file = "nist_total_ms2_neg.rda")

colnames(nist_total_ms2_pos@spectra.info)
colnames(nist_total_ms2_neg@spectra.info)

nist_total_ms2_neg@spectra.info <-
nist_total_ms2_neg@spectra.info %>%
  dplyr::rename(Comments = comments,
                Ionization = ionization,
                Notes = notes
                ) %>%
  dplyr::mutate(Synonyms = synon,
                CCS = NA,
                Ontology = NA) %>%
  dplyr::select(-c("num peaks", collisiongas, mw, sampleinlet, spectrumtype,
                   synon))

nist_total_ms2_pos@spectra.info <-
nist_total_ms2_pos@spectra.info %>%
  dplyr::rename(CCS = ccs,
                Comments = comment,
                Ontology = ontology
  ) %>%
  dplyr::mutate(Ionization = NA,
                Notes = NA) %>%
  dplyr::select(-c("num peaks", retentiontime))


nist_total_ms2 <-
  nist_total_ms2_pos

nist_total_ms2@spectra.info <-
  rbind(nist_total_ms2@spectra.info,
        nist_total_ms2_neg@spectra.info)


nist_total_ms2@spectra.data$Spectra.positive <-
  c(nist_total_ms2@spectra.data$Spectra.positive,
    nist_total_ms2_neg@spectra.data$Spectra.positive)

nist_total_ms2@spectra.data$Spectra.negative <-
  c(nist_total_ms2@spectra.data$Spectra.negative,
    nist_total_ms2_neg@spectra.data$Spectra.negative)

idx1 <-
  match(
    nist_total_ms2@spectra.info$INCHI.ID,
    hmdb_ms1@spectra.info$INCHI.ID,
    incomparables = NA
  )

idx2 <-
  match(
    nist_total_ms2@spectra.info$INCHIKEY.ID,
    hmdb_ms1@spectra.info$INCHIKEY.ID,
    incomparables = NA
  )

idx3 <-
  match(
    nist_total_ms2@spectra.info$CAS.ID,
    hmdb_ms1@spectra.info$CAS.ID,
    incomparables = NA
  )

idx4 <-
  match(
    nist_total_ms2@spectra.info$SMILES.ID,
    hmdb_ms1@spectra.info$SMILES.ID,
    incomparables = NA
  )

idx <-
  data.frame(idx1, idx2, idx3, idx4) %>%
  apply(1, function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    return(x[1])
  })

cbind(
  x = nist_total_ms2@spectra.info$Compound.name,
  y = hmdb_ms1@spectra.info$Compound.name[idx]
) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(y))

colnames(nist_total_ms2@spectra.info)
colnames(hmdb_ms1@spectra.info)

spectra.info <-
  nist_total_ms2@spectra.info

spectra.info <-
  spectra.info %>%
  dplyr::select(-c(CAS.ID, HMDB.ID, KEGG.ID))

spectra.info <-
  data.frame(spectra.info,
             hmdb_ms1@spectra.info[idx, setdiff(colnames(hmdb_ms1@spectra.info), colnames(spectra.info))]) %>%
  dplyr::select(-c(version, Create_date, Updated_date, secondary_accessions))

spectra.info$Compound.name

nist_total_ms2@spectra.info <-
  spectra.info

save(nist_total_ms2, file = "nist_total_ms2.rda")

load("nist_ms2.rda")


match(nist_ms2@spectra.info$HMDB.ID,
      nist_total_ms2@spectra.info$HMDB.ID)

i = 1
j = 73655

which(names(nist_ms2@spectra.data$Spectra.positive) == nist_ms2@spectra.info$Lab.ID[i])
which(names(nist_total_ms2@spectra.data$Spectra.positive) == nist_total_ms2@spectra.info$Lab.ID[j])

which(names(nist_ms2@spectra.data$Spectra.negative) == nist_ms2@spectra.info$Lab.ID[i])
which(names(nist_total_ms2@spectra.data$Spectra.negative) == nist_total_ms2@spectra.info$Lab.ID[j])

masstools::ms2_plot(nist_ms2@spectra.data$Spectra.negative[[1]]$`35,15`,
                    nist_total_ms2@spectra.data$Spectra.negative[[1667]]$`35%`)



grep("Unknown", nist_total_ms2@spectra.info$Compound.name, value = TRUE)

nist_total_ms2 <-
nist_total_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "Unknown"))

nist_total_ms2@spectra.info$Compound.name

grep("^C[0-9]{1,2}H[0-9]{1,2}", nist_total_ms2@spectra.info$Compound.name, value = TRUE)

nist_total_ms2@spectra.info$Compound.name[which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "^C[0-9]{1,2}H[0-9]{1,2}"))]

nist_total_ms2 <-
  nist_total_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "^C[0-9]{1,2}H[0-9]{1,2}"))


nist_total_ms2@spectra.info$Compound.name[which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "base"))]

nist_total_ms2 <-
  nist_total_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "base"))

idx <-
match(nist_total_ms2@spectra.info$HMDB.ID,
      hmdb_ms1@spectra.info$HMDB.ID,
      incomparables = NA)

Compound.name2 <- hmdb_ms1@spectra.info$Compound.name[idx]


Compound.name <-
  data.frame(Compound.name1 = nist_total_ms2@spectra.info$Compound.name,
             Compound.name2) %>%
  apply(1, function(x){
    x <- as.character(x)
    if(!is.na(x[2])){
      return(x[2])
    }else{
      return(x[1])
    }
  })


nist_total_ms2@spectra.info$Compound.name <- Compound.name





nist_total_ms2@spectra.info$Compound.name[which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "CE[0-9]{1,2}"))]
idx <-
  which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "CE[0-9]{1,2}"))
nist_total_ms2@spectra.info$Compound.name[idx] <-
  nist_total_ms2@spectra.info$Compound.name[idx] %>%
  stringr::str_split(";") %>%
  purrr::map(function(x){
    stringr::str_trim(x[1])
  }) %>%
  unlist()



nist_total_ms2@spectra.info$Compound.name[which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "not validated"))]

nist_total_ms2 <-
  nist_total_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "not validated"))


nist_total_ms2@spectra.info$Compound.name[which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "PlaSMA"))]
idx <-
  which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "PlaSMA"))

nist_total_ms2@spectra.info$Compound.name[idx] <-
  nist_total_ms2@spectra.info$Compound.name[idx] %>%
  stringr::str_split(";") %>%
  purrr::map(function(x){
    stringr::str_trim(x[1])
  }) %>%
  unlist()



nist_total_ms2@spectra.info$Compound.name[which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "\\!"))]
idx <-
  which(stringr::str_detect(nist_total_ms2@spectra.info$Compound.name, "\\!"))

nist_total_ms2@spectra.info$Compound.name[idx] <-
  nist_total_ms2@spectra.info$Compound.name[idx] %>%
  stringr::str_split("\\!") %>%
  purrr::map(function(x){
    stringr::str_trim(x[2])
  }) %>%
  unlist()

nist_total_ms2@spectra.info[which(nist_total_ms2@spectra.info == "", arr.ind = TRUE)] <- NA

nist_total_ms2 <-
  nist_total_ms2 %>%
  dplyr::filter(!is.na(Compound.name))


save(nist_total_ms2, file = "nist_total_ms2")


#####remove the plant
grep("Unknown", nist_total_ms2@spectra.info$Compound.name, value = TRUE)

nist_clean_ms2 <-
  nist_total_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "Unknown"))

grep("PlaSMA", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "PlaSMA"))

grep(" CE", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)

nist_clean_ms2@spectra.info$Compound.name[grep(" CE", nist_clean_ms2@spectra.info$Compound.name)] <-
  nist_clean_ms2@spectra.info$Compound.name[grep(" CE", nist_clean_ms2@spectra.info$Compound.name)] %>%
  stringr::str_split(";", n = 2) %>%
  purrr::map(function(x){
  x[1]
  }) %>%
  unlist() %>%
  stringr::str_trim()



grep("MMV", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "MMV"))



grep("NCGC", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)
grep("NCGC", nist_clean_ms2@spectra.info$Compound.name, value = FALSE)

nist_clean_ms2@spectra.info$origin <-
nist_clean_ms2@spectra.info$Comments %>%
  purrr::map(function(x) {
    stringr::str_extract(x, "origin=.+")
  }) %>%
  unlist()

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(is.na(origin))


grep("MS2", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)

nist_clean_ms2@spectra.info$Compound.name[grep("MS2", nist_clean_ms2@spectra.info$Compound.name)] <-
  nist_clean_ms2@spectra.info$Compound.name[grep("MS2", nist_clean_ms2@spectra.info$Compound.name)] %>%
  stringr::str_split(";", n = 2) %>%
  purrr::map(function(x){
    x[1]
  }) %>%
  unlist() %>%
  stringr::str_trim()


###remove compounds in massbank
nist_clean_ms2@spectra.info$origin <-
nist_clean_ms2@spectra.info$Comments %>%
  purrr::map(.f = function(x){
    stringr::str_extract(x, "registered in MassBank")
  }) %>%
  unlist()

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(is.na(origin))

###remove compounds in registered in RIKEN PlaSMA
nist_clean_ms2@spectra.info$origin <-
  nist_clean_ms2@spectra.info$Comments %>%
  purrr::map(.f = function(x){
    stringr::str_extract(x, "registered in RIKEN PlaSMA")
  }) %>%
  unlist()

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(is.na(origin))

###remove compounds in registered in RIKEN PlaSMA
nist_clean_ms2@spectra.info$origin <-
  nist_clean_ms2@spectra.info$Comments %>%
  purrr::map(.f = function(x){
    stringr::str_extract(x, "registered in Respect")
  }) %>%
  unlist()

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(is.na(origin))

###remove compounds in Natural Product
nist_clean_ms2@spectra.info$origin <-
  nist_clean_ms2@spectra.info$Comments %>%
  purrr::map(.f = function(x){
    stringr::str_extract(x, "Natural Product")
  }) %>%
  unlist()

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(is.na(origin))


idx <-
  which(!is.na(nist_clean_ms2@spectra.info$Comments))

new_info <-
idx %>%
  purrr::map(function(i){
    cat(i, " ")
    x <- nist_clean_ms2@spectra.info$Comments[i]
    if(is.na(x)){
      data.frame(
        SMILES.ID = NA,
        CAS.ID = NA,
        CHEBI.ID = NA,
        HMDB.ID = NA,
        KEGG.ID = NA,
        LIPIDMAPS.ID = NA,
        PUBCHEM.ID = NA,
        INCHI.ID = NA,
        Author = NA,
        mz = NA,
        Instrument = NA,
        Instrument_type = NA,
        CE =NA,
        Adduct = NA,
        Precursor_mz= NA,
        Submitter_team = NA,
        accession = NA
      )
    }

    x <-
    x %>%
      stringr::str_split('\"') %>%
      `[[`(1) %>%
      stringr::str_trim()

    x <- x[x!=""]

    SMILES.ID = x[grep("^SMILES=", x)] %>% stringr::str_replace("^SMILES=", "")
    SMILES.ID <-
      ifelse(length(SMILES.ID) == 0, NA, SMILES.ID)
    CAS.ID = x[grep("^cas=", x)] %>% stringr::str_replace("^cas=", "")
    CAS.ID <-
      ifelse(length(CAS.ID) == 0, NA, CAS.ID)
    CHEBI.ID = x[grep("^chebi=", x)] %>% stringr::str_replace("^chebi=", "")
    CHEBI.ID <-
      ifelse(length(CHEBI.ID) == 0, NA, CHEBI.ID)
    HMDB.ID = x[grep("^hmdb=", x)] %>% stringr::str_replace("^hmdb=", "")
    HMDB.ID <-
      ifelse(length(HMDB.ID) == 0, NA, HMDB.ID)
    KEGG.ID = x[grep("^kegg=", x)] %>% stringr::str_replace("^kegg=", "")
    KEGG.ID <-
      ifelse(length(KEGG.ID) == 0, NA, KEGG.ID)
    LIPIDMAPS.ID = x[grep("^lipidmaps=", x)] %>% stringr::str_replace("^lipidmaps=", "")
    LIPIDMAPS.ID <-
      ifelse(length(LIPIDMAPS.ID) == 0, NA, LIPIDMAPS.ID)
    PUBCHEM.ID = x[grep("^pubchem sid=", x)] %>% stringr::str_replace("^pubchem sid=", "")
    PUBCHEM.ID <-
      ifelse(length(PUBCHEM.ID) == 0, NA, PUBCHEM.ID)
    INCHI.ID = x[grep("^InChI=", x)] %>% stringr::str_replace("^InChI=", "")
    INCHI.ID <-
      ifelse(length(INCHI.ID) == 0, NA, INCHI.ID)
    Author = x[grep("^author=", x)] %>% stringr::str_replace("^author=", "")
    Author <-
      ifelse(length(Author) == 0, NA, Author)
    mz = x[grep("^exact mass=", x)] %>% stringr::str_replace("^exact mass=", "")
    mz <-
      ifelse(length(mz) == 0, NA, mz)
    Instrument = x[grep("^instrument=", x)] %>% stringr::str_replace("^instrument=", "")
    Instrument <-
      ifelse(length(Instrument) == 0, NA, Instrument)
    Instrument_type = x[grep("^instrument type=", x)] %>% stringr::str_replace("^instrument type=", "")
    Instrument_type <-
      ifelse(length(Instrument_type) == 0, NA, Instrument_type)
    CE = x[grep("^collision energy=", x)] %>% stringr::str_replace("^collision energy=", "")
    CE <-
      ifelse(length(CE) == 0, NA, CE)
    Adduct = x[grep("^precursor type=", x)] %>% stringr::str_replace("^precursor type=", "")
    Adduct <-
      ifelse(length(Adduct) == 0, NA, Adduct)
    Precursor_mz= x[grep("^precursor m/z=", x)] %>% stringr::str_replace("^precursor m/z=", "")
    Precursor_mz <-
      ifelse(length(Precursor_mz) == 0, NA, Precursor_mz)
    Submitter_team = x[grep("^submitter=", x)] %>% stringr::str_replace("^submitter=", "")
    Submitter_team <-
      ifelse(length(Submitter_team) == 0, NA, Submitter_team)
    accession = x[grep("^accession=", x)] %>% stringr::str_replace("^accession=", "")
    accession <-
      ifelse(length(accession) == 0, NA, accession)

    data.frame(
      SMILES.ID = SMILES.ID,
      CAS.ID = CAS.ID,
      CHEBI.ID = CHEBI.ID,
      HMDB.ID = HMDB.ID,
      KEGG.ID = KEGG.ID,
      LIPIDMAPS.ID = LIPIDMAPS.ID,
      PUBCHEM.ID = PUBCHEM.ID,
      INCHI.ID = INCHI.ID,
      Author = Author,
      mz = mz,
      Instrument = Instrument,
      Instrument_type = Instrument_type,
      CE = CE,
      Adduct = Adduct,
      Precursor_mz= Precursor_mz,
      Submitter_team = Submitter_team,
      accession = accession
    )

  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

new_info <-
data.frame(idx,new_info)

new_info$HMDB.ID <-
new_info$HMDB.ID %>%
  stringr::str_replace("HMDB", "HMDB00")


nist_clean_ms2@spectra.info$LIPIDMAPS.ID <- NA
nist_clean_ms2@spectra.info$Author <- NA
nist_clean_ms2@spectra.info$Submitter_team <- NA
nist_clean_ms2@spectra.info$accession <- NA

nist_clean_ms2@spectra.info$LIPIDMAPS.ID[idx] <- new_info$LIPIDMAPS.ID
nist_clean_ms2@spectra.info$Author[idx] <- new_info$Author
nist_clean_ms2@spectra.info$Submitter_team[idx] <- new_info$Submitter_team
nist_clean_ms2@spectra.info$accession[idx] <- new_info$accession

data.frame(v1 = nist_clean_ms2@spectra.info$HMDB.ID[idx],
      v2 = new_info$HMDB.ID) %>%
  dplyr::filter(!is.na(v2) & is.na(v1))

for (i in 1:length(new_info$idx)) {
  cat(i, " ")
  if (is.na(nist_clean_ms2@spectra.info$SMILES.ID[idx[i]])) {
    nist_clean_ms2@spectra.info$SMILES.ID[idx[i]] <-
      new_info$SMILES.ID[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$CAS.ID[idx[i]])) {
    nist_clean_ms2@spectra.info$CAS.ID[idx[i]] <-
      new_info$CAS.ID[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$CHEBI.ID[idx[i]])) {
    nist_clean_ms2@spectra.info$CHEBI.ID[idx[i]] <-
      new_info$CHEBI.ID[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$HMDB.ID[idx[i]])) {
    nist_clean_ms2@spectra.info$HMDB.ID[idx[i]] <-
      new_info$HMDB.ID[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$KEGG.ID[idx[i]])) {
    nist_clean_ms2@spectra.info$KEGG.ID[idx[i]] <-
        new_info$KEGG.ID[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$PUBCHEM.ID[idx[i]])) {
    nist_clean_ms2@spectra.info$PUBCHEM.ID[idx[i]] <-
      new_info$PUBCHEM.ID[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$INCHI.ID[idx[i]])) {
    nist_clean_ms2@spectra.info$INCHI.ID[idx[i]] <-
      new_info$INCHI.ID[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$mz[idx[i]])) {
    nist_clean_ms2@spectra.info$mz[idx[i]] <-
      new_info$mz[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$Instrument[idx[i]])) {
    nist_clean_ms2@spectra.info$Instrument[idx[i]] <-
      new_info$Instrument[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$Instrument_type[idx[i]])) {
    nist_clean_ms2@spectra.info$Instrument_type[idx[i]] <-
      new_info$Instrument_type[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$CE[idx[i]])) {
    nist_clean_ms2@spectra.info$CE[idx[i]] <-
      new_info$CE[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$Adduct[idx[i]])) {
    nist_clean_ms2@spectra.info$Adduct[idx[i]] <-
      new_info$Adduct[i]
  }

  if (is.na(nist_clean_ms2@spectra.info$Precursor_mz[idx[i]])) {
    nist_clean_ms2@spectra.info$Precursor_mz[idx[i]] <-
      new_info$Precursor_mz[i]
  }

}


nist_clean_ms2@spectra.info$mz <- as.numeric(nist_clean_ms2@spectra.info$mz)
nist_clean_ms2@spectra.info$Precursor_mz <- as.numeric(nist_clean_ms2@spectra.info$Precursor_mz)

nist_clean_ms2@spectra.info$Compound.name


nist_clean_ms2@spectra.info %>%
  dplyr:::filter(is.na(HMDB.ID) & is.na(KEGG.ID)) %>%
  pull(Compound.name) %>%
  unique()




grep("\\#", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)
grep("\\?", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)

nist_clean_ms2 <-
nist_clean_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "\\#")) %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "\\?"))


nist_clean_ms2@spectra.info$Formula


nist_clean_ms2@spectra.info %>%
  dplyr:::filter(is.na(HMDB.ID) & is.na(KEGG.ID)) %>%
  pull(Compound.name) %>%
  unique()

grep("^Man[0-9]{1,3}", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)
grep("^Man[0-9]{1,3}", nist_clean_ms2@spectra.info$Compound.name, value = FALSE)

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(!stringr::str_detect(Compound.name, "^Man[0-9]{1,3}"))



grep("^G[0-9]{1,3}", nist_clean_ms2@spectra.info$Compound.name, value = TRUE)
grep("^G[0-9]{1,3}", nist_clean_ms2@spectra.info$Compound.name, value = FALSE)

nist_clean_ms2@spectra.info$Notes[grep("^G[0-9]{1,3}", nist_clean_ms2@spectra.info$Compound.name, value = FALSE)]

grep("Glycan", nist_clean_ms2@spectra.info$Notes, value = TRUE)

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(!stringr::str_detect(Notes, "Glycan"))

grep("Gly", nist_clean_ms2@spectra.info$Notes, value = TRUE)

nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(!stringr::str_detect(Notes, "Gly"))

grep("gly", nist_clean_ms2@spectra.info$Notes, value = TRUE)


nist_clean_ms2 <-
  nist_clean_ms2 %>%
  dplyr::filter(!is.na(HMDB.ID) | !is.na(KEGG.ID))

save(nist_clean_ms2, file = "nist_clean_ms2.rda")

