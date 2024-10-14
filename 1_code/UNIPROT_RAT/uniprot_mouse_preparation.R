no_source()
setwd(masstools::get_project_wd())
rm(list = ls())
setwd("other_files/UNIPROT_MOUSE/")
library(jsonlite)
library(rjson)
library(tidyverse)

mouse_data <-
  fromJSON(file = "uniprot-compressed_true_download_true_format_json_query__28reviewed_-2023.05.19-05.57.03.54.json")

mouse_uniprot <-
  purrr::map(1:length(mouse_data[[1]]), function(idx) {
    cat(idx, " ")
    x <- mouse_data[[1]][[idx]]
    data.frame(
      UNIPROT = x$primaryAccession,
      uniProtkbId = x$uniProtkbId,
      protein_name = ifelse(
        is.null(x$proteinDescription$recommendedName$fullName$value),
        NA,
        x$proteinDescription$recommendedName$fullName$value
      ),
      gene_name = ifelse(
        is.null(x$genes[[1]]$geneName$value),
        NA,
        x$genes[[1]]$geneName$value
      )
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

save(mouse_uniprot, file = "mouse_uniprot")



