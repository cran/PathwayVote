run_enrichment <- function(gene_list, databases, readable = FALSE, verbose = FALSE) {
  enrich_results <- list()
  if ("Reactome" %in% databases) {
    if (verbose) message("  Analyzing Reactome...")
    enrich_results$Reactome <- as.data.frame(
      enrichPathway(gene = gene_list, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = readable)
    )
  }
  if ("GO" %in% databases) {
    if (verbose) message("  Analyzing GO...")
    enrich_results$GO <- as.data.frame(
      enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = readable)
    )
  }
  if ("KEGG" %in% databases) {
    if (verbose) message("  Analyzing KEGG...")
    kegg_enrich <- enrichKEGG(gene = gene_list, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
    if (readable) {
      kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }
    enrich_results$KEGG <- as.data.frame(kegg_enrich)
  }
  if (verbose) message("  Enrichment completed.")
  return(enrich_results)
}

gene_filter <- function(eQTM, stat_threshold, distance) {
  if (!inherits(eQTM, "eQTM")) {
    stop("eQTM must be an eQTM object")
  }
  eQTM_data <- getData(eQTM)

  eQTM_filtered <- eQTM_data[abs(eQTM_data$statistics) >= stat_threshold & eQTM_data$distance <= distance, ]

  return(list(entrez = unique(na.omit(eQTM_filtered$entrez)), p_values = eQTM_filtered$p_value))
}

filter_gene_lists <- function(eQTM, stat_grid, distance_grid, overlap_threshold = 0.5, verbose = FALSE) {
  if (!inherits(eQTM, "eQTM")) stop("Input must be an eQTM object.")
  data <- getData(eQTM)

  gene_lists <- list()
  params <- list()

  for (r in stat_grid) {
    for (d in distance_grid) {
      subset <- data[abs(data$statistics) > r & data$distance < d, ]
      if (nrow(subset) == 0) next

      map <- split(subset$entrez, subset$cpg)
      map <- lapply(map, function(ids) unique(na.omit(as.numeric(ids))))
      map <- map[lengths(map) > 0]

      if (length(map) == 0) next

      gene_list <- unique(unlist(map))
      gene_list <- gene_list[!is.na(gene_list)]
      if (length(gene_list) == 0) next

      overlap_pass <- TRUE
      for (existing in gene_lists) {
        jaccard <- length(intersect(existing, gene_list)) / length(union(existing, gene_list))
        if (jaccard >= overlap_threshold) {
          overlap_pass <- FALSE
          break
        }
      }
      if (!overlap_pass) next

      if (verbose) {
        message(sprintf("  Accepted: %d genes, %d CpGs (stat > %.2f, dist < %d)",
                        length(gene_list), length(map), r, d))
      }

      gene_lists[[length(gene_lists) + 1]] <- gene_list
      params[[length(params) + 1]] <- c(stat = r, d = d)
    }
  }

  return(list(
    gene_lists = gene_lists,
    params = params
  ))
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
    all_pathways[[db]] <- unique(unlist(lapply(enrich_results, function(x) x[[db]]$ID)))
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
      for (pathway in pathways) {
        mat[pathway, j] <- ifelse(pathway %in% current_df$ID,
                                  current_df$pvalue[current_df$ID == pathway],
                                  1)
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
      if (nrow(enrich_results[[i]][[db]]) > 0) {
        df <- enrich_results[[i]][[db]][, c("ID", "Description", "geneID")]
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
                                   min_genes = 2,
                                   prune_strategy = c("cuberoot", "fixed"),
                                   fixed_value = 3,
                                   verbose = TRUE) {
  prune_strategy <- match.arg(prune_strategy)

  n_runs <- length(enrich_results)
  min_vote_support <- switch(
    prune_strategy,
    cuberoot = max(1, floor(n_runs^(1/3))),
    fixed = fixed_value
  )

  if (verbose) {
    message(sprintf("Prune strategy = '%s', total runs = %d", prune_strategy, n_runs))
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
    message(sprintf("Pruned pathways: %d retained out of %d total (min votes = %d, min genes = %d)",
                    total_after, total_before, min_vote_support, min_genes))
  }

  pruned_results <- lapply(enrich_results, function(x) {
    lapply(x, function(df) {
      if (is.null(df) || !"ID" %in% colnames(df) || !"Count" %in% colnames(df)) return(df)
      df[df$ID %in% keep_ids & df$Count >= min_genes, , drop = FALSE]
    })
  })

  return(pruned_results)
}
