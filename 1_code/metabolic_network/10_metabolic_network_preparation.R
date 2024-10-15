no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

library(tidyverse)

load("2_data/HMDB/MS1/hmdb_ms1.rda")
load("2_data/KEGG/kegg_ms1.rda")

setwd("2_data/metabolic_network/")

load("kegg_rpair_database_human")
load("hmdb_rpair_database_human")
load("bigg_rpair_database_human")
load("metanetx_rpair_database")
load("modelseed_rpair_database")
load("recon3_rpair_database")
load("reactome_rpair_database_human")
load("rhea_rpair_database")

dim(kegg_rpair_database_human)
dim(hmdb_rpair_database_human)
dim(bigg_rpair_database_human)
dim(metanetx_rpair_database)
dim(modelseed_rpair_database)
dim(recon3_rpair_database)
dim(reactome_rpair_database_human)
dim(rhea_rpair_database)

kegg_rpair_database_human <-
  kegg_rpair_database_human %>%
  dplyr::mutate(reaction_database = "KEGG")

hmdb_rpair_database_human <-
  hmdb_rpair_database_human %>%
  dplyr::mutate(reaction_database = "HMDB")

bigg_rpair_database_human <-
  bigg_rpair_database_human %>%
  dplyr::mutate(reaction_database = "BIGG")

metanetx_rpair_database <-
  metanetx_rpair_database %>%
  dplyr::mutate(reaction_database = "METANETX")

modelseed_rpair_database <-
  modelseed_rpair_database %>%
  dplyr::mutate(reaction_database = "MODELSEED")

recon3_rpair_database <-
  recon3_rpair_database %>%
  dplyr::mutate(reaction_database = "RECON3")

reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  dplyr::mutate(reaction_database = "REACTOME")

rhea_rpair_database <-
  rhea_rpair_database %>%
  dplyr::mutate(reaction_database = "RHEA")

colnames(kegg_rpair_database_human)

hmdb_rpair_database_human$reaction_RHEA_ID <- NA
modelseed_rpair_database$reaction_RHEA_ID <- NA
recon3_rpair_database$reaction_RHEA_ID <- NA

reactome_rpair_database_human <-
  reactome_rpair_database_human %>%
  dplyr::mutate(
    reaction_KEGG_ID = NA,
    reaction_RHEA_ID = NA,
    reaction_Enzyme_EC_number = NA,
    reaction_Enzyme_UNIPROT_ID = NA
  )

rhea_rpair_database <-
  rhea_rpair_database %>%
  dplyr::mutate(reaction_RHEA_ID = NA)

rhea_rpair_database[, colnames(kegg_rpair_database_human)]

metabolic_network_edge <-
  rbind(
    kegg_rpair_database_human,
    hmdb_rpair_database_human[, colnames(kegg_rpair_database_human)],
    bigg_rpair_database_human[, colnames(kegg_rpair_database_human)],
    metanetx_rpair_database[, colnames(kegg_rpair_database_human)],
    modelseed_rpair_database[, colnames(kegg_rpair_database_human)],
    recon3_rpair_database[, colnames(kegg_rpair_database_human)],
    reactome_rpair_database_human[, colnames(kegg_rpair_database_human)],
    rhea_rpair_database[, colnames(kegg_rpair_database_human)]
  )

metabolic_network_edge[metabolic_network_edge == ""] <- NA
metabolic_network_edge[metabolic_network_edge == "null"] <- NA

sum(duplicated(metabolic_network_edge$rpair_id))

sum(is.na(metabolic_network_edge$from_compound_name))
sum(is.na(metabolic_network_edge$to_compound_name))

library(plyr)

temp <-
  metabolic_network_edge %>%
  plyr::dlply(.variables = .(rpair_id))

temp <-
  temp %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }

    x$reaction_KEGG_ID <-
      x$reaction_KEGG_ID[!is.na(x$reaction_KEGG_ID)] %>%
      stringr::str_split("\\{\\}") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_RHEA_ID <-
      x$reaction_RHEA_ID[!is.na(x$reaction_RHEA_ID)] %>%
      stringr::str_split("\\{\\}") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_Enzyme_EC_number <-
      x$reaction_Enzyme_EC_number[!is.na(x$reaction_Enzyme_EC_number)] %>%
      stringr::str_split("\\{\\}") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_Enzyme_UNIPROT_ID <-
      x$reaction_Enzyme_UNIPROT_ID[!is.na(x$reaction_Enzyme_UNIPROT_ID)] %>%
      stringr::str_split("\\{\\}") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = "{}")

    x$reaction_database <-
      x$reaction_database[!is.na(x$reaction_database)] %>%
      stringr::str_split("\\{\\}") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = "{}")

    x %>%
      dplyr::distinct(rpair_id,
                      .keep_all = TRUE)
  }) %>%
  dplyr::bind_rows()

metabolic_network_edge <-
  temp

LAB_ID <-
  seq_len(nrow(metabolic_network_edge)) %>%
  purrr::map(function(i) {
    if ((i / 10) %in% seq(1, 100000)) {
      cat(i, " ")
    }
    if (!is.na(metabolic_network_edge$from_compound_HMDB_ID[i])) {
      from_compound_LAB_ID <-
        metabolic_network_edge$from_compound_HMDB_ID[i]
    } else{
      from_compound_LAB_ID <-
        metabolic_network_edge$from_compound_KEGG_ID[i]
    }

    if (!is.na(metabolic_network_edge$to_compound_HMDB_ID[i])) {
      to_compound_LAB_ID <-
        metabolic_network_edge$to_compound_HMDB_ID[i]
    } else{
      to_compound_LAB_ID <-
        metabolic_network_edge$to_compound_KEGG_ID[i]
    }

    data.frame(from_compound_LAB_ID,
               to_compound_LAB_ID)

  }) %>%
  dplyr::bind_rows()

metabolic_network_edge <-
  data.frame(metabolic_network_edge,
             LAB_ID) %>%
  tibble::as_tibble() %>%
  dplyr::select(
    rpair_id,
    from_compound_LAB_ID,
    from_compound_HMDB_ID,
    from_compound_KEGG_ID,
    from_compound_name,
    from_compound_mz,
    from_compound_LAB_ID,
    to_compound_LAB_ID,
    to_compound_HMDB_ID,
    to_compound_KEGG_ID,
    to_compound_name,
    to_compound_mz,
    dplyr::everything()
  )

sum(is.na(metabolic_network_edge$from_compound_name))
sum(is.na(metabolic_network_edge$to_compound_name))

temp1 <-
  metabolic_network_edge[, c(
    "rpair_id",
    "from_compound_name",
    "from_compound_HMDB_ID",
    "from_compound_KEGG_ID"
  )] %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Formula", "mz")],
                   by = c("from_compound_HMDB_ID" = "HMDB.ID")) %>%
  dplyr::left_join(kegg_ms1@spectra.info[, c("KEGG.ID", "Formula", "mz")],
                   by = c("from_compound_KEGG_ID" = "KEGG.ID")) %>%
  dplyr::select(from_compound_name, Formula.x, Formula.y, mz.x, mz.y) %>%
  dplyr::rename(name = from_compound_name)

temp2 <-
  metabolic_network_edge[, c("rpair_id",
                             "to_compound_name",
                             "to_compound_HMDB_ID",
                             "to_compound_KEGG_ID")] %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Formula", "mz")],
                   by = c("to_compound_HMDB_ID" = "HMDB.ID")) %>%
  dplyr::left_join(kegg_ms1@spectra.info[, c("KEGG.ID", "Formula", "mz")],
                   by = c("to_compound_KEGG_ID" = "KEGG.ID")) %>%
  dplyr::select(to_compound_name, Formula.x, Formula.y, mz.x, mz.y) %>%
  dplyr::rename(name = to_compound_name)


temp <-
  rbind(temp1,
        temp2) %>%
  dplyr::distinct(name, .keep_all = TRUE)

temp <-
  temp %>%
  dplyr::mutate(nchar1 = nchar(Formula.x),
                nchar2 = nchar(Formula.y)) %>%
  dplyr::arrange(nchar1, nchar2) %>%
  as.data.frame()


remove_compund <-
  c(
    "Hydrogen Ion",
    "Hydrogen",
    "Helium",
    "Iron",
    "Oxygen",
    "Phosphate",
    "I(-)",
    "Pyrophosphate",
    "Superoxide",
    "Nitric oxide",
    "Calcium",
    "Lithium",
    "Beryllium",
    "Water",
    "Nitrite",
    "NH4OH",
    "Carbon dioxide",
    "Hydrogen sulfide",
    "Fluoride",
    "Thiocyanate",
    "Sulfide",
    "Cyanide",
    "Potassium",
    "Ammonium",
    "Chloride ion",
    "Cotinine methonium ion",
    "Glutathione episulfonium ion",
    "Calcium ionophore",
    "Magnesium ionophore IV",
    "Ammonia",
    "Mercury",
    "Magnesium",
    "trdox",
    "UWM6",
    "Silver",
    "Zinc",
    "Cobalt",
    "Nitrate",
    "Cadmium",
    "Copper",
    "Fe2+",
    "Sodium",
    "Lipid A",
    "Nitrogen",
    "Hydrogen peroxide",
    "Hydrogen carbonate",
    "Carbon monoxide",
    "Hydrogen cyanide",
    "Hydrogen selenide",
    "Hydrogen cyanide",
    "Bromide",
    "Formaldehyde",
    "Hydrochloric acid",
    "Arsenite",
    "Iodine",
    "Iodide",
    "Tungsten",
    "Hydroxide",
    "Manganese",
    "Selenium",
    "Chromium",
    "Lead",
    "Methane",
    "Nitrous oxide",
    "Cyanate",
    "Carbon tetrachloride",
    "Perflutren",
    "Perchloroethylene",
    "Chloroform",
    "Chloromethane",
    "Trichloroethylene"
  )

metabolic_network_edge <-
  metabolic_network_edge %>%
  dplyr::filter(!from_compound_name %in% remove_compund &
                  !to_compound_name %in% remove_compund)

metabolic_network_edge <-
  metabolic_network_edge %>%
  dplyr::filter(from_compound_LAB_ID != to_compound_LAB_ID)

####remove some high connection node
high_degree_name <-
  c(
    metabolic_network_edge$from_compound_LAB_ID,
    metabolic_network_edge$to_compound_LAB_ID
  ) %>%
  data.frame(name = .) %>%
  dplyr::count(name) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(n >= 100) %>%
  pull(name)

high_degree_name

hmdb_ms1@spectra.info[match(high_degree_name, hmdb_ms1@spectra.info$HMDB.ID),
                      c("HMDB.ID", "KEGG.ID", "Compound.name")]

# get_remove_index <-
#   function(rpair_database,
#            metabolite_name = "Adenosine triphosphate",
#            paired_metabolites = c("ADP",
#                                   "Adenosine monophosphate",
#                                   "Acetyl adenylate",
#                                   "Propinol adenylate")) {
#     remove_idx <-
#       seq_len(nrow(rpair_database)) %>%
#       purrr::map(function(i) {
#         if (!rpair_database$from_compound_name[i] %in% metabolite_name &
#             !rpair_database$to_compound_name[i] %in% metabolite_name) {
#           return(NULL)
#         } else{
#           # cat(i, " ")
#           x <- c(rpair_database$from_compound_name[i],
#                  rpair_database$to_compound_name[i])
#           x <- x[x != metabolite_name]
#           if (!x %in% paired_metabolites) {
#             return(i)
#           } else{
#             return(NULL)
#           }
#         }
#       }) %>%
#       unlist()
#     remove_idx
#   }



####Node data
metabolic_network_node <-
  rbind(
    metabolic_network_edge %>%
      dplyr::select(contains("from_compound")) %>%
      dplyr::rename(
        LAB_ID = from_compound_LAB_ID,
        HMDB_ID = from_compound_HMDB_ID,
        KEGG_ID = from_compound_KEGG_ID,
        Compound_name = from_compound_name,
        mz = from_compound_mz
      ),
    metabolic_network_edge %>%
      dplyr::select(contains("to_compound")) %>%
      dplyr::rename(
        LAB_ID = to_compound_LAB_ID,
        HMDB_ID = to_compound_HMDB_ID,
        KEGG_ID = to_compound_KEGG_ID,
        Compound_name = to_compound_name,
        mz = to_compound_mz
      )
  ) %>%
  dplyr::distinct(LAB_ID, .keep_all = TRUE)

dim(metabolic_network_edge)
dim(metabolic_network_node)

library(tidygraph)

metabolic_network_edge <-
  metabolic_network_edge %>%
  dplyr::mutate(from = from_compound_LAB_ID,
                to = to_compound_LAB_ID) %>%
  dplyr::select(from, to, dplyr::everything())

metabolic_network_node <-
  metabolic_network_node %>%
  dplyr::mutate(name = LAB_ID) %>%
  dplyr::select(name, dplyr::everything())

metabolic_network <-
  tidygraph::tbl_graph(nodes = metabolic_network_node,
                       edges = metabolic_network_edge,
                       directed = FALSE) %>%
  activate(what = "nodes") %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

save(metabolic_network_edge, file = "metabolic_network_edge")
save(metabolic_network_node, file = "metabolic_network_node")
save(metabolic_network, file = "metabolic_network")

#####visulization
library(ggraph)

node_data <-
  metabolic_network %>%
  activate(what = "nodes") %>%
  as.data.frame() %>%
  tibble::as_tibble()

edge_data <-
  metabolic_network %>%
  activate(what = "edges") %>%
  as.data.frame() %>%
  tibble::as_tibble()

node_data[which(node_data$degree == 1),]

# metabolic_network_edge[which(metabolic_network_edge$from == "HMDB0000538"),]
# which(metabolic_network_edge$from == "HMDB0001341")
library(graphlayouts)

metabolic_network2 <-
  metabolic_network
# activate(what = "nodes") %>%
# filter(degree > 1 & degree < 20)

plot <-
  metabolic_network2 %>%
  ggraph::ggraph(layout = "stress") +
  geom_edge_link0(color = "black") +
  geom_node_point(size = 1) +
  ggraph::theme_graph()

#####new metabolic network with enzyme
metabolic_network_edge %>%
  dim()

metabolic_network_node %>%
  dim()

metabolic_network_edge_with_protein <-
  seq_len(nrow(metabolic_network_edge)) %>%
  purrr::map(function(i) {
    if ((i / 10) %in% seq(1, 100000, 1)) {
      cat(i, " ")
    }
    if (is.na(metabolic_network_edge$reaction_Enzyme_UNIPROT_ID[i])) {
      x <-
        metabolic_network_edge[i, ] %>%
        dplyr::select(
          from = from_compound_LAB_ID,
          to = to_compound_LAB_ID,
          rpair_id,
          from_ID = from_compound_LAB_ID,
          to_ID = to_compound_LAB_ID,
          reaction_database,
          reaction_KEGG_ID,
          reaction_RHEA_ID
        ) %>%
        dplyr::mutate(from_class = "protein",
                      to_class = "metabolite")
      return(x)
    }

    enzyme <-
      stringr::str_split(metabolic_network_edge$reaction_Enzyme_UNIPROT_ID[i],
                         "\\{\\}")[[1]]

    lapply(enzyme, function(x) {
      rbind(
        data.frame(
          from = x,
          to = metabolic_network_edge$from[i],
          rpair_id = paste(x, metabolic_network_edge$from[i], sep = "_"),
          from_ID = x,
          to_ID = metabolic_network_edge$from_compound_LAB_ID[i],
          reaction_database = metabolic_network_edge$reaction_database[i],
          reaction_KEGG_ID = metabolic_network_edge$reaction_KEGG_ID[i],
          reaction_RHEA_ID = metabolic_network_edge$reaction_RHEA_ID[i],
          from_class = "protein",
          to_class = "metabolite"
        ),
        data.frame(
          from = x,
          to = metabolic_network_edge$to[i],
          rpair_id = paste(x, metabolic_network_edge$to[i], sep = "_"),
          from_ID = x,
          to_ID = metabolic_network_edge$to_compound_LAB_ID[i],
          reaction_database = metabolic_network_edge$reaction_database[i],
          reaction_KEGG_ID = metabolic_network_edge$reaction_KEGG_ID[i],
          reaction_RHEA_ID = metabolic_network_edge$reaction_RHEA_ID[i],
          from_class = "protein",
          to_class = "metabolite"
        )
      )
    }) %>%
      dplyr::bind_rows()
  }) %>%
  dplyr::bind_rows()

save(metabolic_network_edge_with_protein, file = "metabolic_network_edge_with_protein")


metabolic_network_node_with_protein <-
  data.frame(metabolic_network_node,
             UNIPROT_ID = NA,
             class = "metabolite")

protein_id <-
  metabolic_network_edge$reaction_Enzyme_UNIPROT_ID[!is.na(metabolic_network_edge$reaction_Enzyme_UNIPROT_ID)] %>%
  stringr::str_split("\\{\\}") %>%
  unlist() %>%
  unique()

metabolic_network_node_with_protein2 <-
  data.frame(
    name = protein_id,
    LAB_ID = protein_id,
    HMDB_ID = NA,
    KEGG_ID = NA,
    Compound_name = NA,
    mz = NA,
    UNIPROT_ID = protein_id,
    class = "protein"
  )

metabolic_network_node_with_protein <-
  rbind(metabolic_network_node_with_protein,
        metabolic_network_node_with_protein2) %>%
  tibble::as_tibble()

save(metabolic_network_node_with_protein, file = "metabolic_network_node_with_protein")

metabolic_network_with_protein <-
  tidygraph::tbl_graph(nodes = metabolic_network_node_with_protein,
                       edges = metabolic_network_edge_with_protein,
                       directed = FALSE)

save(metabolic_network_with_protein, file = "metabolic_network_with_protein")
