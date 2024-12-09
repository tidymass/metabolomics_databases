no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)

rm(list = ls())

### to get KEGG database
library(KEGGgraph)
library(KEGGREST)
library(tidyverse)

### load data
kegg_metabolite <- readxl::read_xlsx("3_data_analysis/KEGG/metabolite/kegg_metabolite.xlsx")
kegg_drug <- readxl::read_xlsx("3_data_analysis/KEGG/drug/kegg_drug.xlsx")

dim(kegg_metabolite)

dim(kegg_drug)

###16409 metabolites and 9218 drugs

dir.create("3_data_analysis/KEGG/compound", showWarnings = FALSE)

setwd("3_data_analysis/KEGG/compound")

colnames(kegg_metabolite)
colnames(kegg_drug)

colnames(kegg_metabolite) == colnames(kegg_drug)

intersect_kegg_id <-
intersect(kegg_metabolite$KEGG.ID,
kegg_drug$KEGG.ID)

intersect_kegg_id[1]   

length(intersect_kegg_id)

###1625 overlapped compounds
kegg_metabolite %>% 
dplyr::filter(KEGG.ID == intersect_kegg_id[1])

kegg_drug %>% 
dplyr::filter(KEGG.ID == intersect_kegg_id[1])

rbind(
    dplyr::filter(kegg_metabolite, KEGG.ID %in% intersect_kegg_id[1]),
dplyr::filter(kegg_drug, KEGG.ID %in% intersect_kegg_id[1])
)

####for the overlapped compounds, we will combine the information from metabolite and drug
overlapped_kegg_info <-
purrr::map(intersect_kegg_id, function(i){
    cat(i, " ")
    kegg_metabolite_info <- dplyr::filter(kegg_metabolite, KEGG.ID == i)
    kegg_drug_info <- dplyr::filter(kegg_drug, KEGG.ID == i)
    kegg_metabolite_info$Synonyms <-
    paste(kegg_metabolite_info$Synonyms,
    kegg_drug_info$Synonyms, sep = "{}")
    kegg_metabolite_info
}) %>% 
do.call(rbind, .) %>% 
as.data.frame()

kegg_metabolite <-
kegg_metabolite %>%
dplyr::filter(!KEGG.ID %in% intersect_kegg_id)

kegg_drug <-
kegg_drug %>%
dplyr::filter(!KEGG.ID %in% intersect_kegg_id)

kegg_compound <- rbind(kegg_metabolite, kegg_drug, overlapped_kegg_info)

write.csv(kegg_compound, "kegg_compound.csv", row.names = FALSE)

library(metid)

kegg_compound_ms1 <-
    construct_database(
        path = ".",
        version = "2024-12-03",
        metabolite.info.name = "kegg_compound.csv",
        source = "KEGG",
        link = "https://www.genome.jp/kegg",
        creater = "Xiaotao Shen",
        email = "xiaotao.shen@outlook.com",
        rt = FALSE,
        threads = 3
    )

save(kegg_compound_ms1, file = "kegg_compound_ms1.rda")

# load("../HMDB/MS1/hmdb_ms1.rda")

# intersect(
#     colnames(kegg_ms1@spectra.info),
#     colnames(hmdb_ms1@spectra.info)
# )
# setdiff(
#     colnames(hmdb_ms1@spectra.info),
#     colnames(kegg_ms1@spectra.info)
# )

# source(here::here("1_code/3_utils.R"))

# kegg_ms1 <-
#     update_metid_database_info(
#         database = kegg_ms1,
#         ref_database = hmdb_ms1,
#         by = c(
#             "Compound.name",
#             "CAS.ID",
#             "HMDB.ID",
#             "KEGG.ID",
#             "PUBCHEM.ID",
#             "CHEBI.ID",
#             "DRUGBANK.ID"
#         ),
#         combine_columns = c(
#             "CAS.ID",
#             "HMDB.ID",
#             "KEGG.ID",
#             "PUBCHEM.ID",
#             "CHEBI.ID",
#             "DRUGBANK.ID",
#             "From_human"
#         ),
#         new_columns = c(
#             "status",
#             "IUPAC_name",
#             "Traditional_IUPAC_name",
#             "SMILES.ID",
#             "INCHI.ID",
#             "INCHIKEY.ID",
#             "Kingdom",
#             "Super_class",
#             "Class",
#             "Sub_class",
#             "State",
#             "Biospecimen_locations",
#             "Cellular_locations",
#             "Tissue_locations",
#             "CHEMSPIDER.ID",
#             "FOODB.ID",
#             "BIOCYC.ID",
#             "BIGG.ID",
#             "WIKIPEDIA.ID",
#             "METLIN.ID"
#         )
#     )

# load(here::here("2_data/source_system/source_system.rda"))

# library(tidyverse)
# library(tidyselect)
# library(metid)

# kegg_ms1 <-
#     update_metid_database_source_system(
#         database = kegg_ms1,
#         source_system = source_system,
#         by = c("CAS.ID", "HMDB.ID", "KEGG.ID"),
#         prefer = "database"
#     )



# save(kegg_ms1, file = "kegg_ms1.rda")

# load("kegg_ms1.rda")

# idx <-
#     grep("C[0-9]{4,8}", kegg_ms1@spectra.info$Compound.name, value = FALSE)


# if (length(idx) > 0) {
#     kegg_ms1@spectra.info$Compound.name[idx] <-
#         kegg_ms1@spectra.info$Synonyms[idx] %>%
#         stringr::str_split("\\{\\}") %>%
#         purrr::map(function(x) {
#             x[2]
#         }) %>%
#         unlist()
# }

# save(kegg_ms1, file = "kegg_ms1.rda")
