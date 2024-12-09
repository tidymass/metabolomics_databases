convert_bigg_universal2metid <-
  function(data,
           path = ".",
           threads = 5) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    ####Formula and mz
    data[which(data == "", arr.ind = TRUE)] <- NA

    data <-
      data %>%
      dplyr::filter(!is.na(formula))

    ####add formula to H or remove H based on charge
    message("Extract formula...")
    pb <-
      progress::progress_bar$new(total = nrow(data))

    formula <-
      seq_len(nrow(data)) %>%
      purrr::map(function(idx) {
        cat(idx, " ")
        pb$tick()
        charge = data$charges[idx]
        formula <- data$formula[idx]

        if (is.na(charge)) {
          return(data$formula[idx])
        }

        if (charge == "0") {
          return(data$formula[idx])
        }

        if (length(grep("\\{\\}", charge)) == 0) {
          charge <- as.numeric(charge)
          adduct = paste0(ifelse(charge > 0, "M-", "M+"), abs(charge), "H")
          formula <-
            masstools::sum_formula(formula = formula, adduct = adduct)
          return(formula)
        }

        if (length(grep("\\{\\}", charge)) > 0) {
          charge <- as.numeric(stringr::str_split(charge, "\\{\\}")[[1]])
          formula <- stringr::str_split(formula, "\\{\\}")[[1]]
          if (length(charge) != length(formula)) {
            return(NA)
          }

          if (any(charge) == 0) {
            return(formula[which(charge == 0)])
          }

          adduct = paste0(ifelse(charge[1] > 0, "M-", "M+"), abs(charge[1]), "H")
          formula <-
            masstools::sum_formula(formula = formula[1], adduct = adduct)
          return(formula)
        }
      }) %>%
      unlist()

    data$formula <- formula

    data$HMDB <-
      data$HMDB %>%
      purrr::map(function(x) {
        if (is.na(x)) {
          return(x)
        }

        x <- stringr::str_split(x, "\\{\\}")[[1]]
        x <- stringr::str_replace(x, "HMDB", "HMDB00")
        x <- paste(x, collapse = "{}")
        x
      }) %>%
      unlist()

    data <-
      data %>%
      dplyr::filter(!is.na(formula))

    ###add mz
    message("Calculating mz...")
    pb <-
      progress::progress_bar$new(total = nrow(data))
    mz <-
      seq_len(nrow(data)) %>%
      purrr::map(function(i) {
        pb$tick()
        temp <-
          tryCatch(
            Rdisop::getMass(Rdisop::getMolecule(data$formula[i])),
            error = function(e)
              NA
          )
      }) %>%
      unlist()

    data$mz <- mz

    data <-
      data %>%
      dplyr::filter(!is.na(mz))

    data <-
      data %>%
      dplyr::rename(
        BIGG.ID = bigg_id,
        Compound.name = name,
        Formula = formula,
        BIOCYC.ID = BioCyc,
        CHEBI.ID = CHEBI,
        HMDB.ID = HMDB,
        INCHIKEY.ID = InChI_Key,
        KEGG.ID = KEGG,
        METANETX.ID = MetaNetX,
        REACTOME.ID = Reactome,
        SEED.ID = SEED,
        KEGG_DRUG.ID = KEGG_Drug,
        KEGG_GLYCAN.ID = KEGG_Glycan,
        LIPIDMAPS.ID = LipidMaps
      ) %>%
      dplyr::mutate(
        Lab.ID = BIGG.ID,
        CAS.ID = NA,
        RT = NA,
        mz.pos = NA,
        mz.neg = NA,
        Submitter = "BIGG"
      ) %>%
      dplyr::select(
        Lab.ID,
        Compound.name,
        mz,
        RT,
        CAS.ID,
        HMDB.ID,
        KEGG.ID,
        Formula,
        mz.pos,
        mz.neg,
        Submitter,
        everything()
      )

    data$Synonyms <-
      data$Compound.name %>%
      stringr::str_replace_all("; ", "\\{\\}")

    data$Compound.name <-
      data$Synonyms %>%
      stringr::str_split("\\{\\}") %>%
      lapply(function(x) {
        x[1]
      }) %>%
      unlist()

    temp_file <- tempfile()
    dir.create(temp_file, showWarnings = FALSE)
    readr::write_csv(x = data, file = file.path(temp_file, "data.csv"))


    bigg_ms1 =
      metid::construct_database(
        path = temp_file,
        version = as.character(Sys.Date()),
        metabolite.info.name = "data.csv",
        source = "BIGG",
        link = "http://bigg.ucsd.edu/",
        creater = "Xiaotao Shen",
        email = "shenxt@stanford.edu",
        rt = FALSE,
        threads = threads
      )

    unlink(file.path(temp_file, "data.csv"))
    unlink(temp_file)

    save(bigg_ms1, file = file.path(path, "bigg_ms1"))
    invisible(bigg_ms1)
  }