##
no_source()
setwd(masstools::get_project_wd())
library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
rm(list = ls())

library(RJSONIO)
library(readr)

library(xml2)

source("R/14_FOODB.R")

load("2_data/FOODB/MS1/foodb_ms1.rda")

setwd('2_data/FOODB/MS2/')

dir.create("ms2_data")

223:length(foodb_ms1@spectra.info$Lab.ID) %>%
  purrr::walk(
    .f = function(i) {
      cat(i, " ")
      x <- foodb_ms1@spectra.info$Lab.ID[i]
      ms2 <-
        tryCatch(request_foodb_compound_ms2(compound_id = x),
                 error = function(e) NULL)
      if(!is.null(ms2)){
        save(ms2, file = file.path("ms2_data", x))
      }

    }
  )



compound <-
  readr::read_csv("foodb_2020_04_07_csv/Compound.csv")

public_id <-
  compound$public_id

# for(i in 50001:length(public_id)){
#   cat(i, " ")
#   result <-
#     tryCatch(request_foodb_compound(compound_id = public_id[i]), error = function(e) NULL)
#   if(is.null(result)){
#     next()
#   }
#   save(result, file = file.path("compound", public_id[i]), compress = "xz")
# }

file <- dir("compound/")

foodb <-
  1:length(file) %>%
  purrr::map(function(i) {
    cat(i, " ")
    load(file.path("compound", file[i]))
    result
  })

foodb_result <-
  foodb %>%
  purrr::map(function(x) {
    cat(x$accession, " ")
    Kingdom <-
      ifelse(is.null(x$taxonomy$kingdom), NA, x$taxonomy$kingdom)
    Super_class <-
      ifelse(is.null(x$taxonomy$super_class),
             NA,
             x$taxonomy$super_class)
    Class <-
      ifelse(is.null(x$taxonomy$class), NA, x$taxonomy$class)
    Sub_class <-
      ifelse(is.null(x$taxonomy$sub_class), NA, x$taxonomy$sub_class)
    State <- ifelse(is.null(x$state), NA, x$state)
    HMDB.ID <- ifelse(is.null(x$hmdb_id), NA, x$hmdb_id)
    PUBCHEM.ID <-
      ifelse(is.null(x$pubchem_compound_id), NA, x$pubchem_compound_id)
    CHEMSPIDER.ID <-
      ifelse(is.null(x$chemspider_id), NA, x$chemspider_id)
    KEGG.ID <- ifelse(is.null(x$kegg_id), NA, x$kegg_id)
    CHEBI.ID <- ifelse(is.null(x$chebi_id), NA, x$chebi_id)
    BIOCYC.ID <- ifelse(is.null(x$biocyc_id), NA, x$biocyc_id)
    HET.ID <- ifelse(is.null(x$het_id), NA, x$het_id)
    WIKIPEDIA.ID <- ifelse(is.null(x$wikipidia), NA, x$wikipidia)
    VMH.ID <- ifelse(is.null(x$vmh_id), NA, x$vmh_id)
    Synonyms = ifelse(is.null(x$synonyms), NA, paste(unname(unlist(x$synonyms)), collapse = "{}"))
    foods <- x$foods
    if(is.null(foods)){
      Food_name <- NA
      Food_type <- NA
      Food_category <- NA
      Food_scientific_name <- NA
    }else{
      if(class(foods)[1] == "list"){
        foods <-
          foods %>%
          lapply(function(y){
            data.frame(name = ifelse(is.null(y$name), NA, y$name),
                       food_type = ifelse(is.null(y$food_type), NA, y$food_type),
                       category = ifelse(is.null(y$category), NA, y$category),
                       name_scientific = ifelse(is.null(y$name_scientific), NA, y$name_scientific))
          }) %>%
          dplyr::bind_rows() %>%
          as.data.frame()

        Food_name <- as.character(foods[,"name"]) %>%
          paste(collapse = "{}")
        Food_type <- as.character(foods[,"food_type"]) %>%
          paste(collapse = "{}")
        Food_category <- as.character(foods[,"category"]) %>%
          paste(collapse = "{}")
        Food_scientific_name <-
          as.character(foods[,"name_scientific"]) %>%
          paste(collapse = "{}")

      }else{
        colnames(foods) <- paste("V", 1:ncol(foods), sep = "")

        foods <-
          foods[1:4,,drop = FALSE] %>%
          as.data.frame() %>%
          apply(2, function(y){
            y <-
              y %>%
              lapply(function(z){
                if(is.null(z)) {
                  return(NA)
                } else{
                  return(z)
                }
              }) %>%
              unlist()
            y
          })
        Food_name <- as.character(foods["name", ]) %>%
          paste(collapse = "{}")
        Food_type <- as.character(foods["food_type", ]) %>%
          paste(collapse = "{}")
        Food_category <- as.character(foods["category", ]) %>%
          paste(collapse = "{}")
        Food_scientific_name <-
          as.character(foods["name_scientific", ]) %>%
          paste(collapse = "{}")
      }
    }

    data.frame(
      Lab.ID = x$accession,
      Create_date = ifelse(is.null(x$creation_date), NA, x$creation_date),
      Updated_date = ifelse(is.null(x$update_date), NA, x$update_date),
      Compound.name = ifelse(is.null(x$name), NA, x$name),
      Description = ifelse(is.null(x$description), NA, x$description),
      Formula = ifelse(is.null(x$chemical_formula), NA, x$chemical_formula),
      Synonyms = Synonyms,
      Average.mass =  ifelse(
        is.null(x$average_molecular_weight),
        NA,
        x$average_molecular_weight
      ),
      mz = ifelse(
        is.null(x$monisotopic_moleculate_weight),
        NA,
        x$monisotopic_moleculate_weight
      ),
      IUPAC_name = ifelse(is.null(x$iupac_name), NA, x$iupac_name),
      Traditional_IUPAC_name = ifelse(is.null(x$traditional_iupac), NA, x$traditional_iupac),
      CAS.ID = ifelse(is.null(x$cas_registry_number), NA, x$cas_registry_number),
      SMILES.ID = ifelse(is.null(x$smiles), NA, x$smiles),
      INCHI.ID = ifelse(is.null(x$inchi), NA, x$inchi),
      INCHIKEY.ID = ifelse(is.null(x$inchikey), NA, x$inchikey),
      Kingdom = Kingdom,
      Super_class = Super_class,
      Class = Class,
      Sub_class = Sub_class,
      State = State,
      FOODB.ID = ifelse(is.null(x$accession), NA, x$accession),
      HMDB.ID = HMDB.ID,
      PUBCHEM.ID = PUBCHEM.ID,
      CHEMSPIDER.ID = CHEMSPIDER.ID,
      KEGG.ID = KEGG.ID,
      CHEBI.ID = CHEBI.ID,
      BIOCYC.ID = BIOCYC.ID,
      HET.ID = HET.ID,
      WIKIPEDIA.ID = WIKIPEDIA.ID,
      VMH.ID = VMH.ID,
      Food_name,
      Food_type,
      Food_category,
      Food_scientific_name
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()


foodb_result$mz <-
  as.numeric(foodb_result$mz)

foodb_result$RT <- NA
foodb_result$From_food <- TRUE
foodb_result$mz.pos = NA
foodb_result$mz.neg = NA
foodb_result$Submitter = "FOODB"

foodb_result <-
foodb_result %>%
  dplyr::filter(!is.na(mz) & !is.na(Formula))


openxlsx::write.xlsx(foodb_result, file = "foodb_result.xlsx", asTable = TRUE)

library(metid)

foodb_ms1 <-
  construct_database(
    path = ".",
    version = "2022-04-21",
    metabolite.info.name = "foodb_result.xlsx",
    source = "FOODB",
    link = "https://foodb.ca/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

save(foodb_ms1, file = "foodb_ms1.rda", compress = "xz")
