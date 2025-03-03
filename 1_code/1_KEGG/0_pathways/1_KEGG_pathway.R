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

dir.create("2_data/KEGG/pathway", showWarnings = FALSE)
setwd("2_data/KEGG/pathway")

###reference pathways
# reference_pathway_id <-
#   keggList(database = "pathway") %>%
#   names() %>%
#   unique()

# kegg_reference_pathway_database <-
#       pbapply::pblapply(reference_pathway_id, function(x) {
#         Sys.sleep(0.5)
#         KEGGREST::keggGet(dbentries = x)[[1]]
#     })

# pathway_id <-
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         unname(x$ENTRY)
#       }) %>%
#       unlist()

#  pathway_name <-
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         unname(x$PATHWAY_MAP)
#       }) %>%
#       unlist()
    
#     description =
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         unname(x$DESCRIPTION)
#       })
    
#     pathway_class =
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         unname(x$CLASS)
#       })
    
#     gene_list =
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         gene = x$GENE
#         if (is.null(gene)) {
#           return(data.frame())
#         }
#         data.frame(
#           KEGG.ID = gene[seq(1, length(gene) - 1, by = 2)],
#           Gene.name = gene[seq(2, length(gene), by = 2)],
#           stringsAsFactors = FALSE
#         )
#       })
    
#     compound_list =
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         data.frame(
#           KEGG.ID = names(x$COMPOUND),
#           Compound.name = x$COMPOUND,
#           stringsAsFactors = FALSE
#         )
#       })
    
#     reference_list =
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         purrr::map(
#           x$REFERENCE,
#           .f = function(y) {
#             y = lapply(y, function(z) {
#               if (length(z) > 1) {
#                 paste(z, collapse = "{}")
#               } else{
#                 z
#               }
#             })
#             y = unlist(y)
#             if (any(names(y) == "JOURNAL")) {
#               names(y)[names(y) == "JOURNAL"] = "JOURNAL1"
#               c(y, JOURNAL2 = "")
#             }
#           }
#         ) %>%
#           do.call(rbind, .) %>%
#           as.data.frame()
#       })
    
#     related_disease =
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         data.frame(
#           Disease.ID = names(x$DISEASE),
#           Disease.name = x$DISEASE,
#           stringsAsFactors = FALSE
#         )
#       })
    
    
#     related_module =
#       kegg_reference_pathway_database %>%
#       purrr::map(function(x) {
#         data.frame(
#           Module.ID = names(x$MODULE),
#           Module.name = x$MODULE,
#           stringsAsFactors = FALSE
#         )
#       })
    
# library(metpath)

#     pathway =
#       new(
#         Class = "pathway_database",
#         database_info = list(source = paste("KEGG", i, sep = "_"),
#                              version = as.character(Sys.Date())),
#         pathway_id = pathway_id,
#         pathway_name = pathway_name,
#         describtion = description,
#         pathway_class = pathway_class,
#         gene_list = gene_list,
#         compound_list = compound_list,
#         protein_list = list(),
#         reference_list = reference_list,
#         related_disease = related_disease,
#         related_module = related_module
#       )
    
#     if (length(pathway@gene_list) == 0) {
#       pathway@gene_list = vector(mode = "list",
#                                  length = length(pathway@pathway_id)) %>%
#         purrr::map(function(x) {
#           x = data.frame()
#           x
#         })
#     }
    
#     if (length(pathway@compound_list) == 0) {
#       pathway@compound_list = vector(mode = "list",
#                                      length = length(pathway@pathway_id)) %>%
#         purrr::map(function(x) {
#           x = data.frame()
#           x
#         })
#     }
    
#     if (length(pathway@protein_list) == 0) {
#       pathway@protein_list = vector(mode = "list",
#                                     length = length(pathway@pathway_id)) %>%
#         purrr::map(function(x) {
#           x = data.frame()
#           x
#         })
#     }

# kegg_reference_pathway <-
# pathway


# save(kegg_reference_pathway, file = "kegg_reference_pathway.rda")

load("kegg_reference_pathway.rda")



####specific organism pathways
###get the organism list
organism_list <-
  keggList(database = "organism") %>%
  as.data.frame()

retry_request <- function(dbentries, retries = 3) {
  for (i in seq_len(retries)) {
    tryCatch({
      return(KEGGREST::keggGet(dbentries)[[1]])
    }, error = function(e) {
      if (i == retries) {
        stop(e)  
      }
      message("Error encountered: ", conditionMessage(e), " - Retrying in 1 minute...")
      Sys.sleep(60)  
    })
  }
}

dir.create("organism_pathways")

for (i in organism_list$organism[8001:9000]) {

  print(paste0(i, " in progress"))

  if(any(dir() == i)) next
  
  pathway_id <- keggList(database = "pathway", organism = i) %>%
    names() %>%
    unique()
  

  kegg_organism_pathway_database <- pbapply::pblapply(pathway_id, function(x) {
    if (match(x, pathway_id) %% 10 == 0) {
      Sys.sleep(5)  
    }
    retry_request(x) 
  })
  
  pathway_id <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      unname(x$ENTRY)
    }) %>%
    unlist()
  
  pathway_name <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      unname(x$PATHWAY_MAP)
    }) %>%
    unlist()
  
  description <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      unname(x$DESCRIPTION)
    })
  
  pathway_class <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      unname(x$CLASS)
    })
  
  gene_list <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      gene = x$GENE
      if (is.null(gene)) {
        return(data.frame())
      }
      data.frame(
        KEGG.ID = gene[seq(1, length(gene) - 1, by = 2)],
        Gene.name = gene[seq(2, length(gene), by = 2)],
        stringsAsFactors = FALSE
      )
    })
  
  compound_list <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      data.frame(
        KEGG.ID = names(x$COMPOUND),
        Compound.name = x$COMPOUND,
        stringsAsFactors = FALSE
      )
    })
  
  reference_list <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      purrr::map(
        x$REFERENCE,
        .f = function(y) {
          y = lapply(y, function(z) {
            if (length(z) > 1) {
              paste(z, collapse = "{}")
            } else {
              z
            }
          })
          y = unlist(y)
          if (any(names(y) == "JOURNAL")) {
            names(y)[names(y) == "JOURNAL"] = "JOURNAL1"
            c(y, JOURNAL2 = "")
          }
        }
      ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    })
  
  related_disease <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      data.frame(
        Disease.ID = names(x$DISEASE),
        Disease.name = x$DISEASE,
        stringsAsFactors = FALSE
      )
    })
  
  related_module <- kegg_organism_pathway_database %>%
    purrr::map(function(x) {
      data.frame(
        Module.ID = names(x$MODULE),
        Module.name = x$MODULE,
        stringsAsFactors = FALSE
      )
    })
  
  
  pathway <- new(
    Class = "pathway_database",
    database_info = list(source = paste("KEGG", i, sep = "_"),
                         version = as.character(Sys.Date())),
    pathway_id = pathway_id,
    pathway_name = pathway_name,
    describtion = description,
    pathway_class = pathway_class,
    gene_list = gene_list,
    compound_list = compound_list,
    protein_list = list(),
    reference_list = reference_list,
    related_disease = related_disease,
    related_module = related_module
  )  

  if (length(pathway@gene_list) == 0) {
    pathway@gene_list <- vector(mode = "list", length = length(pathway@pathway_id)) %>%
      purrr::map(function(x) {
        data.frame()
      })
  }
  
  if (length(pathway@compound_list) == 0) {
    pathway@compound_list <- vector(mode = "list", length = length(pathway@pathway_id)) %>%
      purrr::map(function(x) {
        data.frame()
      })
  }
  
  if (length(pathway@protein_list) == 0) {
    pathway@protein_list <- vector(mode = "list", length = length(pathway@pathway_id)) %>%
      purrr::map(function(x) {
        data.frame()
      })
  }
  
  dir.create(i, showWarnings = FALSE)
  save(pathway, file = file.path(i, "pathway.rda"))
  
}



