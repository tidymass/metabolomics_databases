

no_source()
setwd(masstools::get_project_wd())
rm(list = ls())
setwd("2_data/UNIPROT_RAT/")
library(jsonlite)
library(rjson)
library(tidyverse)

# rat_data <-
#   fromJSON(file = "uniprot-compressed_true_download_true_format_json_query__28reviewed_-2023.05.19-06.10.12.43.json")
#
# rat_uniprot <-
#   purrr::map(1:length(rat_data[[1]]), function(idx) {
#     cat(idx, " ")
#     x <- rat_data[[1]][[idx]]
#     data.frame(UNIPROT = x$primaryAccession,
#                uniProtkbId = x$uniProtkbId,
#                protein_name = x$proteinDescription$recommendedName$fullName$value,
#                gene_name = ifelse(is.null(x$genes[[1]]$geneName$value), NA, x$genes[[1]]$geneName$value))
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# rat_data2 <-
#   fromJSON(file = "uniprot-compressed_true_download_true_format_json_query__28reviewed_-2023.05.19-06.26.19.13.json")
#
# rat_uniprot2 <-
#   purrr::map(1:length(rat_data2[[1]]), function(idx) {
#     cat(idx, " ")
#     x <- rat_data2[[1]][[idx]]
#     data.frame(UNIPROT = x$primaryAccession,
#                uniProtkbId = x$uniProtkbId,
#                protein_name = ifelse(is.null(x$proteinDescription$recommendedName$fullName$value), NA, x$proteinDescription$recommendedName$fullName$value),
#                gene_name = ifelse(is.null(x$genes[[1]]$geneName$value), NA, x$genes[[1]]$geneName$value))
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# rat_uniprot <-
#   rbind(rat_uniprot,
#         rat_uniprot2)
#
# save(rat_uniprot, file = "rat_uniprot")

load("rat_uniprot")

