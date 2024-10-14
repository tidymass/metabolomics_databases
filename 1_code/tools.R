omics_color <-
  c(
    "protein" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1],
    "metabolite" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[2]
  )



my_theme <-
  function(){
    theme_bw() +
      theme(panel.grid.minor = element_blank()
      )
  }
