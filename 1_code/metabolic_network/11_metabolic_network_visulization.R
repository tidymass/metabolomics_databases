no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

library(tidyverse)

source("1_code/tools.R")

setwd("2_data/metabolic_network/")

load("metabolic_network_with_protein.rda")

load("metabolic_network.rda")

#####visualization
library(ggraph)
library(tidygraph)
library(tidyverse)

metabolic_network_with_protein <-
  metabolic_network_with_protein %>%
  activate(what = "nodes") %>%
  dplyr::filter(!is.na(name)) %>% 
  dplyr::filter(name != "") %>% 
  dplyr::mutate(degree = tidygraph::centrality_degree())

node_data <-
  metabolic_network_with_protein %>%
  activate(what = "nodes") %>%
  as.data.frame() %>%
  tibble::as_tibble()

edge_data <-
  metabolic_network_with_protein %>%
  activate(what = "edges") %>%
  as.data.frame() %>%
  tibble::as_tibble()

edge_data$from <-
  node_data$name[edge_data$from]

edge_data$to <-
  node_data$name[edge_data$to]

node_data[which(node_data$degree > 2000),]

######property of network
node_data_protein <-
  node_data %>%
  dplyr::filter(class == "protein")

node_data_metabolite <-
  node_data %>%
  dplyr::filter(class == "metabolite")

library(ggforce)

plot <-
  node_data_protein %>%
  ggplot(aes(degree)) +
  geom_histogram(binwidth = 10,
                 color = "black",
                 fill = omics_color["protein"]) +
  labs(x = "Degree of proteins (Enzymes)",
       y = "Counts") +
  ggforce::facet_zoom(xlim = c(-10, 110)) +
  my_theme()
plot
# ggsave(plot,
#        filename = "protein_degree_distribution.pdf",
#        width = 9,
#        height = 9)


node_data_protein %>% 
  dplyr::filter(degree > 1000)


range(node_data_metabolite$degree)

plot <-
  node_data_metabolite %>%
  ggplot(aes(degree)) +
  geom_histogram(binwidth = 10,
                 color = "black",
                 fill = omics_color["metabolite"]) +
  labs(x = "Degree of metbolites (Metabolic features)",
       y = "Counts") +
  ggforce::facet_zoom(xlim = c(-10, 110)) +
  my_theme()
plot
# ggsave(plot,
#        filename = "metabolite_degree_distribution.pdf",
#        width = 9,
#        height = 9)


node_data_metabolite %>% 
  dplyr::filter(degree > 6000)

library(graphlayouts)

node_data <-
  node_data %>%
  dplyr::filter(degree < 100)

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% node_data$name &
                  to %in% node_data$name)

node_data <-
  node_data %>%
  dplyr::filter(name %in% c(edge_data$from, edge_data$to))

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% node_data$name,
                to %in% node_data$name)


###output for cytoscape
# write.csv(node_data, file = "node_data.csv", row.names = FALSE)
# write.csv(edge_data, file = "edge_data.csv", row.names = FALSE)
# 
# write.table(
#   node_data,
#   file = "node_data.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )
# write.table(
#   edge_data,
#   file = "edge_data.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )


######
library(ggraph)
plot <-
metabolic_network %>% 
  activate("nodes") %>% 
  dplyr::mutate(degree2 = centrality_degree()) %>% 
  dplyr::filter(degree2 > 10 & degree2 < 80) %>% 
  ggraph(layout = "fr") +
  geom_edge_link(edge_colour = "grey") +
  geom_node_point(aes(size = degree),
                  alpha = 0.8,
                  color = "#ff595e") +
  theme_void()
  
plot

# ggsave(plot,
#        filename = "metabolic_network.pdf",
#        width = 9,
#        height = 7)

metabolic_network
# 
# library(igraph)
# 
# distance_table <-
#   distance_table(metabolic_network, directed = FALSE)
# 
# all_shortest_paths(
#   graph = g2,
#   from = "1",
#   to = "10",
#   mode = "all"
# )
# 
# distance_data <-
#   seq_along(V(metabolic_network))[-length(V(metabolic_network))] %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     temp_distance <-
#       shortest_paths(
#         graph = metabolic_network,
#         from = V(metabolic_network)[i],
#         to = V(metabolic_network)[(i + 1):length(V(metabolic_network))],
#         mode = "all"
#       )$vpath %>%
#       lapply(function(x) {
#         length(as.character(x)) - 1
#       }) %>%
#       unlist()
#     data.frame(
#       from = vertex_attr(metabolic_network)$name[i],
#       to = vertex_attr(metabolic_network)$name[(i + 1):length(V(metabolic_network))],
#       distance = temp_distance
#     )
#   }) %>%
#   dplyr::bind_rows()
