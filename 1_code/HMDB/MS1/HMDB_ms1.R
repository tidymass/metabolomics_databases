setwd(masstools::get_project_wd())
rm(list = ls())

source("R/3_utils.R")
setwd("other_files/HMDB/MS1/")

library(tidyverse)
library(xml2)
library(stringr)

# hmdb <- read_xml("hmdb_metabolites.xml")
#
# hmdb <- as_list(hmdb)
#
# x = hmdb$hmdb[[1]]
#
# hmdb_metabolite <- vector(mode = "list", length = length(hmdb$hmdb))
#
# for(i in 1:length(hmdb$hmdb)){
#   cat(i, " ")
#   x <- hmdb$hmdb[[i]]
#   biospecimen_locations <-
#     paste(stringr::str_trim(unname(unlist(x$biological_properties$biospecimen_locations))),
#           collapse = "{}")
#   cellular_locations <-
#     paste(stringr::str_trim(unlist(x$biological_properties$cellular_locations)),
#           collapse = "{}")
#   tissue_locations <-
#     paste(stringr::str_trim(unlist(x$biological_properties$tissue_locations)),
#           collapse = "{}")
#   # pathway <-
#   # lapply(x$biological_properties$pathways, function(y){
#   #   pathway_name <- unlist(unlist(y$name))
#   #   pathway_smpdb_id <- unlist(unlist(y$smpdb_id))
#   #   pathway_kegg_id <- unlist(unlist(y$kegg_map_id))
#   #   pathway_name <-
#   #     ifelse(is.null(pathway_name), NA, pathway_name)
#   #   pathway_smpdb_id <-
#   #     ifelse(is.null(pathway_smpdb_id), NA, pathway_smpdb_id)
#   #   pathway_kegg_id <-
#   #     ifelse(is.null(pathway_kegg_id), NA, pathway_kegg_id)
#   #   c(pathway_name = pathway_name,
#   #     pathway_smpdb_id = pathway_smpdb_id,
#   #     pathway_kegg_id = pathway_kegg_id)
#   # }) %>%
#   #   do.call(rbind, .) %>%
#   #   as.data.frame()
#   #
#   # pathway$pathway_name <-
#   #   paste(pathway$pathway_name, collapse = "{}")
#   # pathway$pathway_smpdb_id <-
#   #   paste(pathway$pathway_smpdb_id[!is.na(pathway$pathway_smpdb_id)], collapse = "{}")
#   # pathway$pathway_kegg_id <-
#   #   paste(pathway$pathway_kegg_id[!is.na(pathway$pathway_kegg_id)], collapse = "{}")
#   #
#   # pathway <-
#   # pathway %>%
#   #   dplyr::distinct(pathway_name, .keep_all = TRUE)
#   #
#   # source <- x$ontology[[2]]
#
#   average_molecular_weight = ifelse(is.null(unlist(x$average_molecular_weight)),
#                                     NA,
#                                     as.numeric(unlist(x$average_molecular_weight)))
#   monisotopic_molecular_weight = ifelse(is.null(unlist(x$monisotopic_molecular_weight)),
#                                         NA,
#                                         as.numeric(unlist(x$monisotopic_molecular_weight)))
#
#   secondary_accessions = paste(stringr::str_trim(unlist(x$secondary_accessions)),
#                                collapse = "{}")
#
#   synonyms = paste(stringr::str_trim(unname(unlist(x$synonyms))), collapse = "{}")
#
#
#   hmdb_metabolite[[i]] <-
#     data.frame(version = ifelse(is.null(unlist(x$version)), NA, unlist(x$version)),
#                creation_date = ifelse(is.null(unlist(x$creation_date)), NA, unlist(x$creation_date)),
#                update_date = ifelse(is.null(unlist(x$update_date)), NA, unlist(x$update_date)),
#                accession = ifelse(is.null(unlist(x$accession)), NA, unlist(x$accession)),
#                status = ifelse(is.null(unlist(x$status)), NA, unlist(x$status)),
#                secondary_accessions = secondary_accessions,
#                name = unlist(x$name),
#                description = ifelse(is.null(unlist(x$description)), NA, unlist(x$description)),
#                synonyms = synonyms,
#                chemical_formula = ifelse(is.null(unlist(x$chemical_formula)), NA, unlist(x$chemical_formula)),
#                average_molecular_weight = average_molecular_weight,
#                monisotopic_molecular_weight = monisotopic_molecular_weight,
#                iupac_name = ifelse(is.null(unlist(x$iupac_name)), NA, unlist(x$iupac_name)),
#                traditional_iupac = ifelse(is.null(unlist(x$traditional_iupac)), NA, unlist(x$traditional_iupac)),
#                cas_registry_number = ifelse(is.null(unlist(x$cas_registry_number)), NA, unlist(x$cas_registry_number)),
#                smiles = ifelse(is.null(unlist(x$smiles)), NA, unlist(x$smiles)),
#                inchi = ifelse(is.null(unlist(x$inchi)), NA, unlist(x$inchi)),
#                inchikey = ifelse(is.null(unlist(x$inchikey)), NA, unlist(x$inchikey)),
#                kingdom = ifelse(is.null(unlist(x$taxonomy$kingdom)),
#                                 NA, unlist(x$taxonomy$kingdom)),
#                super_class = ifelse(is.null(unlist(x$taxonomy$super_class)),
#                                     NA, unlist(x$taxonomy$super_class)),
#                class = ifelse(is.null(unlist(x$taxonomy$class)),
#                               NA, unlist(x$taxonomy$class)),
#                sub_class = ifelse(is.null(unlist(x$taxonomy$sub_class)),
#                                   NA, unlist(x$taxonomy$sub_class)),
#                state = ifelse(is.null(unlist(x$state)), NA, unlist(x$state)),
#                biospecimen_locations = biospecimen_locations,
#                cellular_locations = cellular_locations,
#                tissue_locations = tissue_locations,
#                chemspider_id = ifelse(is.null(unlist(x$chemspider_id)), NA, unlist(x$chemspider_id)),
#                drugbank_id = ifelse(is.null(unlist(x$drugbank_id)), NA, unlist(x$drugbank_id)),
#                foodb_id = ifelse(is.null(unlist(x$foodb_id)), NA, unlist(x$foodb_id)),
#                pubchem_compound_id = ifelse(is.null(unlist(x$pubchem_compound_id)), NA, unlist(x$pubchem_compound_id)),
#                chebi_id = ifelse(is.null(unlist(x$chebi_id)), NA, unlist(x$chebi_id)),
#                kegg_id = ifelse(is.null(unlist(x$kegg_id)), NA, unlist(x$kegg_id)),
#                biocyc_id = ifelse(is.null(unlist(x$biocyc_id)), NA, unlist(x$biocyc_id)),
#                bigg_id = ifelse(is.null(unlist(x$bigg_id)), NA, unlist(x$bigg_id)),
#                wikipedia_id = ifelse(is.null(unlist(x$wikipedia_id)), NA, unlist(x$wikipedia_id)),
#                metlin_id = ifelse(is.null(unlist(x$metlin_id)), NA, unlist(x$metlin_id))
#     )
# }
#
# hmdb_metabolite <-
#   hmdb_metabolite %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# save(hmdb_metabolite, file = "hmdb_metabolite", compress = "xz")
# load("hmdb_metabolite")
#
# dim(hmdb_metabolite)
#
# hmdb_metabolite <-
#   hmdb_metabolite %>%
#   dplyr::filter(!is.na(monisotopic_molecular_weight))
#
# hmdb_metabolite[which(hmdb_metabolite == "", arr.ind = TRUE)] <- NA
#
# hmdb_metabolite <-
#   hmdb_metabolite %>%
#   dplyr::rename(
#     mz = average_molecular_weight,
#     Create_date = creation_date,
#     Updated_date = update_date,
#     Lab.ID = accession,
#     Compound.name = name,
#     Description = description,
#     Synonyms = synonyms,
#     Formula = chemical_formula,
#     IUPAC_name = iupac_name,
#     Traditional_IUPAC_name = traditional_iupac,
#     CAS.ID = cas_registry_number,
#     SMILES.ID = smiles,
#     INCHI.ID = inchi,
#     INCHIKEY.ID = inchikey,
#     Kingdom = kingdom,
#     Super_class = super_class,
#     Class = class,
#     Sub_class = sub_class,
#     State = state,
#     Biospecimen_locations = biospecimen_locations,
#     Cellular_locations = cellular_locations,
#     Tissue_locations = tissue_locations,
#     CHEMSPIDER.ID = chemspider_id,
#     DRUGBANK.ID = drugbank_id,
#     FOODB.ID = foodb_id,
#     PUBCHEM.ID = pubchem_compound_id,
#     CHEBI.ID = chebi_id,
#     KEGG.ID = kegg_id,
#     BIOCYC.ID = biocyc_id,
#     BIGG.ID = bigg_id,
#     WIKIPEDIA.ID = wikipedia_id,
#     METLIN.ID = metlin_id
#   ) %>%
#   dplyr::mutate(
#     HMDB.ID = Lab.ID,
#     RT = NA,
#     mz.pos = NA,
#     mz.neg = NA,
#     Submitter = "HMDB",
#     From_human = "Yes"
#   ) %>%
#   dplyr::select(
#     Lab.ID,
#     Compound.name,
#     mz,
#     RT,
#     CAS.ID,
#     HMDB.ID,
#     KEGG.ID,
#     Formula,
#     mz.pos,
#     mz.neg,
#     Submitter,
#     everything()
#   )
#
# load("../../KEGG/kegg_ms1.rda")
#
# idx1 <-
#   match(hmdb_metabolite$CAS.ID,
#         kegg_ms1@spectra.info$CAS.ID,
#         incomparables = NA)
#
# idx2 <-
#   match(hmdb_metabolite$KEGG.ID,
#         kegg_ms1@spectra.info$KEGG.ID,
#         incomparables = NA)
#
# idx <-
#   data.frame(idx1, idx2) %>%
#   keep_one_from_multiple()
#
# KEGG.ID <-
#   data.frame(hmdb_metabolite$KEGG.ID, kegg_ms1@spectra.info$KEGG.ID[idx]) %>%
#   keep_one_from_multiple()
# hmdb_metabolite$KEGG.ID <- KEGG.ID
#
# CAS.ID <-
#   data.frame(hmdb_metabolite$CAS.ID, kegg_ms1@spectra.info$CAS.ID[idx]) %>%
#   keep_one_from_multiple()
# hmdb_metabolite$CAS.ID <- CAS.ID
#
# openxlsx::write.xlsx(hmdb_metabolite,
#                      file = "hmdb_metabolite.xlsx",
#                      asTable = TRUE, overwrite = TRUE)
#
library(metid)

hmdb_ms1 <-
  construct_database(
    path = ".",
    version = "2022-04-11",
    metabolite.info.name = "hmdb_metabolite.xlsx",
    source = "HMDB",
    link = "https://hmdb.ca/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

hmdb_ms1@spectra.info$mz <-
hmdb_ms1@spectra.info$monisotopic_molecular_weight

load(here::here("other_files/source_system/source_system.rda"))

library(tidyverse)
library(tidyselect)
library(metid)

hmdb_ms1 <-
  update_metid_database_source_system(
    database = hmdb_ms1,
    source_system = source_system,
    by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
    prefer = "database"
  )

save(hmdb_ms1, file = "hmdb_ms1.rda")
