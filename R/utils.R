safe_setup_plan <- function(workers) {
  os <- .Platform$OS.type
  tryCatch({
    if (os == "windows") {
      future::plan(future::multisession, workers = workers)
    } else {
      future::plan(future::multicore, workers = workers)
    }
  }, error = function(e) {
    message("Failed to setup parallel backend. Falling back to sequential. Reason: ", e$message)
    future::plan(future::sequential)
  })
}

check_ewas_cpg_match <- function(ewas_data, eqtm, threshold = 0.8) {
  ewas_ids <- ewas_data[[1]]
  eqtm_cpgs <- getData(eqtm)$cpg
  match_rate <- mean(ewas_ids %in% eqtm_cpgs, na.rm = TRUE)
  match_rate >= threshold
}

generate_k_grid_fdr_guided <- function(ewas_data,
                                       rank_column,
                                       grid_size = 5,
                                       fdr_cutoff = 0.05,
                                       expand_factor = exp(1),
                                       verbose = FALSE) {
  pvals <- ewas_data[[rank_column]]
  fdr_vals <- p.adjust(pvals, method = "BH")

  sig_indices <- which(fdr_vals <= fdr_cutoff)
  n_sig <- length(sig_indices)

  if (n_sig < 10) {
    stop(paste0(
      "FDR-guided k_grid generation aborted: only ", n_sig,
      " CpGs passed FDR < ", fdr_cutoff, ".\n",
      "You may explicitly specify `k_grid` manually if you wish to continue, ",
      "but enrichment results may be unreliable due to extremely weak signal."
    ))
  }

  max_k <- min(length(pvals), ceiling(n_sig * expand_factor))
  min_k <- floor(n_sig * 0.25)

  # Use log scale to spread grid between min_k and max_k
  k_grid <- round(exp(seq(log(min_k), log(max_k), length.out = grid_size)))

  if (verbose) {
    message("FDR-guided k_grid: ", paste(k_grid, collapse = ", "),
            " (n_sig = ", n_sig, ")")
  }

  return(k_grid)
}

jaccard_dist <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

select_gene_lists_entropy_auto <- function(gene_lists, grid_size = 5, overlap_threshold = 0.7, verbose = FALSE) {
  all_genes <- unlist(gene_lists)
  gene_freq <- table(all_genes)

  entropy_score <- sapply(gene_lists, function(gset) {
    sum(log(1 / gene_freq[gset]), na.rm = TRUE)
  })

  instability <- sapply(seq_along(gene_lists), function(i) {
    mean(sapply(setdiff(seq_along(gene_lists), i), function(j) 1 - jaccard_dist(gene_lists[[i]], gene_lists[[j]])))
  })

  gene_probe_count <- table(unlist(gene_lists))
  penalty_score <- sapply(gene_lists, function(gset) {
    sum(log1p(gene_probe_count[gset]), na.rm = TRUE)
  })

  total_score <- scale(entropy_score) - scale(instability) - scale(penalty_score)

  remaining <- seq_along(gene_lists)
  selected <- c()

  while (length(remaining) > 0) {
    best_idx <- remaining[which.max(total_score[remaining])]
    selected <- c(selected, best_idx)

    overlap <- sapply(remaining, function(i) {
      jaccard_dist(gene_lists[[best_idx]], gene_lists[[i]])
    })

    remaining <- setdiff(remaining, remaining[overlap >= overlap_threshold])

    if (verbose) {
      message(sprintf("Selected gene list %d, removed %d overlapping lists.", best_idx, sum(overlap >= overlap_threshold)))
    }
  }

  return(gene_lists[selected])
}

# Prepares user data for enricher in a compliant and performant way
prepare_enrichment_data <- function(databases, organism = "human", verbose = FALSE) {

  user_data_list <- list()

  if ("GO" %in% databases) {
    if (verbose) message("Preparing GO annotation data...")
    # Get GO to Entrez Gene mappings
    go_data <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                                     columns = c("GOALL", "ONTOLOGYALL"),
                                     keytype = "ENTREZID")

    term2gene <- go_data[, c("GOALL", "ENTREZID")]
    term2name <- go_data[, c("GOALL", "ONTOLOGYALL")] # We just need a placeholder name, ID is fine

    # Get GO term descriptions
    go_terms <- AnnotationDbi::select(GO.db::GO.db, keys=unique(term2gene$GOALL), columns=c("TERM"), keytype="GOID")
    term2name <- go_terms[, c("GOID", "TERM")]
    colnames(term2name) <- c("TERM", "NAME")

    user_data_list$GO <- list(TERM2GENE = term2gene, TERM2NAME = term2name)
  }

  if ("KEGG" %in% databases) {
    if (verbose) message("Preparing KEGG annotation data...")
    # Use the exported function to download/load KEGG data
    kegg_data <- clusterProfiler::download_KEGG(species = "hsa")
    user_data_list$KEGG <- list(TERM2GENE = kegg_data$KEGGPATHID2EXTID, TERM2NAME = kegg_data$KEGGPATHID2NAME)
  }

  if ("Reactome" %in% databases) {
    if (verbose) message("Preparing Reactome annotation data...")
    # Use the reactome.db package for compliant data access
    term2gene <- AnnotationDbi::select(reactome.db,
                                       keys = keys(reactome.db, "ENTREZID"),
                                       columns = "PATHID",
                                       keytype = "ENTREZID")
    colnames(term2gene) <- c("GENE", "TERM")
    term2gene <- term2gene[, c("TERM", "GENE")] # Ensure correct column order

    term2name <- AnnotationDbi::select(reactome.db,
                                       keys = unique(term2gene$TERM),
                                       columns = "PATHNAME",
                                       keytype = "PATHID")
    colnames(term2name) <- c("TERM", "NAME")

    user_data_list$Reactome <- list(TERM2GENE = term2gene, TERM2NAME = term2name)
  }

  return(user_data_list)
}

run_enrichment <- function(gene_list,
                           databases = c("Reactome", "GO", "KEGG"),
                           universe = NULL,
                           readable = FALSE,
                           preloaded_data = list()) {
  enrich_results <- list()

  if (length(gene_list) == 0) return(enrich_results)

  for (db in databases) {
    if (is.null(preloaded_data[[db]])) next

    # Use the generic, performant enricher function with pre-loaded data
    res <- clusterProfiler::enricher(
      gene = gene_list,
      universe = universe,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500,
      TERM2GENE = preloaded_data[[db]]$TERM2GENE,
      TERM2NAME = preloaded_data[[db]]$TERM2NAME
    )

    if (!is.null(res)) {
      if (readable) {
        # Ensure the OrgDb is available for ID mapping
        res <- tryCatch(setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) res)
      }
      enrich_results[[db]] <- as.data.frame(res)
    } else {
      enrich_results[[db]] <- NULL
    }
  }

  return(enrich_results)
}

generate_gene_lists_grid <- function(eQTM, stat_grid, distance_grid, verbose = FALSE) {
  if (!inherits(eQTM, "eQTM")) stop("Input must be an eQTM object.")
  data <- getData(eQTM)

  gene_lists <- list()
  params <- list()

  for (r in stat_grid) {
    for (d in distance_grid) {
      subset <- data[abs(data$statistics) > r & data$distance < d, ]
      if (nrow(subset) == 0) next

      map <- split(subset$entrez, subset$cpg)
      map <- lapply(map, function(ids) unique(na.omit(as.character(ids))))
      map <- map[lengths(map) > 0]

      gene_list <- unique(unlist(map))
      if (length(gene_list) == 0) next

      gene_lists[[length(gene_lists) + 1]] <- gene_list
      params[[length(params) + 1]] <- c(stat = r, d = d)

      if (verbose) {
        message(sprintf("Generated: %d genes, %d CpGs (|stat| > %.2f, dist < %d)",
                        length(gene_list), length(map), r, d))
      }
    }
  }

  return(list(gene_lists = gene_lists, params = params))
}

harmonic_mean_p <- function(p_values) {
  valid_p <- p_values[!is.na(p_values) & p_values >= 0 & p_values <= 1]
  if (length(valid_p) == 0) return(1)
  if (length(valid_p) == 1) return(valid_p[1])
  return(length(valid_p) / sum(1 / valid_p))
}

combine_enrichment_results <- function(enrich_results, databases, verbose = FALSE) {
  if (verbose) message("Combining enrichment results...")

  all_pathways <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Extracting pathways for %s...", db))
    all_pathways[[db]] <- unique(unlist(lapply(enrich_results, function(x) {
      if (!is.null(x[[db]]) && is.data.frame(x[[db]]) && "ID" %in% colnames(x[[db]])) {
        x[[db]]$ID
      } else {
        NULL
      }
    })))
  }

  p_value_matrices <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Building p-value matrix for %s...", db))
    pathways <- all_pathways[[db]]
    if (length(pathways) == 0 || all(is.na(pathways))) {
      warning("No pathways found for database: ", db)
      next
    }

    mat <- matrix(NA, nrow = length(pathways), ncol = length(enrich_results))
    rownames(mat) <- pathways

    for (j in seq_along(enrich_results)) {
      current_df <- enrich_results[[j]][[db]]
      if (!is.null(current_df) && is.data.frame(current_df) && nrow(current_df) > 0) {
        for (pathway in pathways) {
          if (pathway %in% current_df$ID) {
            mat[pathway, j] <- current_df$pvalue[current_df$ID == pathway]
          } else {
            mat[pathway, j] <- 1
          }
        }
      } else {
        mat[, j] <- 1
      }
    }
    p_value_matrices[[db]] <- mat
  }

  combined_p <- list()
  for (db in databases) {
    mat <- p_value_matrices[[db]]
    if (is.null(mat) || nrow(mat) == 0) {
      warning("Skipping database ", db, ": no p-values to combine.")
      combined_p[[db]] <- numeric(0)
      next
    }

    if (verbose) message(sprintf("  Combining p-values for %s using HMP...", db))
    combined_vec <- c()
    for (pathway in rownames(mat)) {
      combined_vec[pathway] <- harmonic_mean_p(mat[pathway, ])
    }

    if (length(combined_vec) == 0 || all(is.na(combined_vec))) {
      warning("No valid combined p-values for database: ", db)
      combined_p[[db]] <- numeric(0)
    } else {
      combined_vec <- combined_vec[order(combined_vec)]
      combined_p[[db]] <- p.adjust(combined_vec, method = "BH")
    }
  }

  pathway_info <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Collecting pathway info for %s...", db))
    pathway_info[[db]] <- list()
    for (i in seq_along(enrich_results)) {
      current_df <- enrich_results[[i]][[db]]
      if (!is.null(current_df) && is.data.frame(current_df) && nrow(current_df) > 0) {
        df <- current_df[, c("ID", "Description", "geneID")]
        for (j in 1:nrow(df)) {
          pathway_id <- df$ID[j]
          if (is.null(pathway_info[[db]][[pathway_id]])) {
            pathway_info[[db]][[pathway_id]] <- list(Description = df$Description[j], geneIDs = character())
          }
          genes <- unlist(strsplit(df$geneID[j], "/"))
          pathway_info[[db]][[pathway_id]]$geneIDs <- union(pathway_info[[db]][[pathway_id]]$geneIDs, genes)
        }
      }
    }
  }

  result_tables <- list()
  for (db in databases) {
    if (verbose) message(sprintf("  Generating result table for %s...", db))
    if (length(combined_p[[db]]) == 0) {
      result_tables[[db]] <- data.frame()
      next
    }

    pathway_df <- data.frame(
      ID = names(pathway_info[[db]]),
      Description = sapply(pathway_info[[db]], function(x) x$Description),
      geneID = sapply(pathway_info[[db]], function(x) paste(x$geneIDs, collapse = "/")),
      stringsAsFactors = FALSE
    )

    result_tables[[db]] <- data.frame(
      ID = names(combined_p[[db]]),
      p.adjust = combined_p[[db]],
      row.names = NULL
    )

    result_tables[[db]] <- merge(result_tables[[db]], pathway_df, by = "ID", all.x = TRUE)
    result_tables[[db]] <- result_tables[[db]][order(result_tables[[db]]$p.adjust), ]
  }

  if (verbose) message("Result combination completed.")
  return(result_tables)
}

prune_pathways_by_vote <- function(enrich_results,
                                   fixed_prune = 3,
                                   min_genes = 2,
                                   verbose = FALSE) {
  n_runs <- length(enrich_results)
  min_vote_support <- if (is.null(fixed_prune)) {
    max(1, floor(n_runs^(1/3)))
  } else {
    fixed_prune
  }

  if (verbose) {
    message(sprintf("Pathway pruning: min votes = %d, min genes = %d", min_vote_support, min_genes))
  }

  all_hits <- unlist(lapply(enrich_results, function(x) {
    unlist(lapply(x, function(df) {
      if (is.null(df) || !"ID" %in% colnames(df) || !"Count" %in% colnames(df)) return(character(0))
      df$ID[df$Count >= min_genes]
    }))
  }))

  total_before <- length(unique(all_hits))
  freq_table <- table(all_hits)
  keep_ids <- names(freq_table)[freq_table >= min_vote_support]
  total_after <- length(keep_ids)

  if (verbose) {
    message(sprintf("Pruned pathways: %d retained out of %d total enriched pathways", total_after, total_before))
  }

  pruned_results <- lapply(enrich_results, function(x) {
    lapply(x, function(df) {
      if (is.null(df) || !"ID" %in% colnames(df) || !"Count" %in% colnames(df)) return(df)
      df[df$ID %in% keep_ids & df$Count >= min_genes, , drop = FALSE]
    })
  })

  return(pruned_results)
}
