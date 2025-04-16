run_benchmark_simulations <- function(
    n_repeats = 100,
    seed_base = 20240407,
    hierarchy_table,
    pathway2gene,
    eQTM_db_PathwayVote,
    voting_params = list(
      k_values = seq(50, 450, 100),
      stat_grid = seq(0.05, 0.6, 0.15),
      distance_grid = seq(10000, 100000, 30000),
      overlap_threshold = 0.8,
      min_vote_support = 3,
      min_genes_per_hit = 2
    ),
    fdr_cutoff = 0.05,
    verbose = TRUE
) {
  library(dplyr)
  library(tibble)
  library(purrr)

  all_metrics <- list()
  per_pathway_log <- list()

  # ----------------------
  extract_additional_metrics <- function(result_df, sim_data, pathway2gene, fdr_cutoff = 0.05) {
    detected_pathways <- result_df$ID[result_df$p.adjust <= fdr_cutoff]
    pathway_sizes <- sapply(detected_pathways, function(pid) length(pathway2gene[[pid]]))
    avg_pathway_size <- mean(pathway_sizes)

    covered_genes <- unique(unlist(strsplit(result_df$geneID[result_df$p.adjust <= fdr_cutoff], "/")))
    coverage_rate <- length(intersect(covered_genes, sim_data$signal_entrez_genes)) / length(sim_data$signal_entrez_genes)

    tibble(
      avg_pathway_size = avg_pathway_size,
      signal_gene_coverage = coverage_rate
    )
  }
  # ----------------------

  for (i in seq_len(n_repeats)) {
    seed <- seed_base + i
    if (verbose) message("🧪 Running simulation ", i, " (seed = ", seed, ")")

    # get top tier pathways with sufficient genes but unrelated pathways
    pathway_info <- as.list(select_diverse_top_pathways(hierarchy_table, pathway2gene, seed = seed))
    sim_data <- simulate_benchmark_dataset(
      pathway_info = pathway_info,
      pathway2gene = pathway2gene,
      eQTM_df = getData(eQTM_db_PathwayVote),
      signal_genes_per_pathway = 25,
      n_noise = 10000,
      n_per_gene = 3,
      signal_score_mean = 3,
      signal_score_sd = 0.5,
      correlation_threshold = 0,
      p_value_threshold = 1e-7,
      seed = seed
    )

    # --- Conventional ---
    conventional_result <- cpg2gene_enrichment(sim_data$ewas_data, top_k = max(voting_params$k_values))
    eval_c <- benchmark_enrichment_results(
      result_df = conventional_result,
      truth_pathway_ids = sim_data$truth_pathway_ids,
      hierarchy_table = hierarchy_table,
      fdr_cutoff = fdr_cutoff,
      include_descendants = TRUE,
      filter_descendants_by_gene = TRUE,
      signal_entrez_genes = sim_data$signal_entrez_genes,
      pathway2gene = pathway2gene
    )

    # --- PathwayVote ---
    vote_result <- pathway_vote(
      ewas_data = sim_data$ewas_data,
      eQTM = eQTM_db_PathwayVote,
      k_values = voting_params$k_values,
      stat_grid = voting_params$stat_grid,
      distance_grid = voting_params$distance_grid,
      overlap_threshold = voting_params$overlap_threshold,
      rank_column = "score",
      use_abs = TRUE,
      rank_decreasing = TRUE,
      min_vote_support = voting_params$min_vote_support,
      min_genes_per_hit = voting_params$min_genes_per_hit,
      databases = "Reactome",
      readable = FALSE,
      verbose = FALSE
    )

    eval_v <- benchmark_enrichment_results(
      result_df = vote_result$result_tables$Reactome,
      truth_pathway_ids = sim_data$truth_pathway_ids,
      hierarchy_table = hierarchy_table,
      fdr_cutoff = fdr_cutoff,
      include_descendants = TRUE,
      filter_descendants_by_gene = TRUE,
      signal_entrez_genes = sim_data$signal_entrez_genes,
      pathway2gene = pathway2gene
    )

    # ---- per-pathway log ----
    per_pathway_c <- eval_c$per_pathway %>%
      mutate(method = "Conventional", iteration = i)

    per_pathway_v <- eval_v$per_pathway %>%
      mutate(method = "PathwayVote", iteration = i)

    per_pathway_log[[i]] <- bind_rows(per_pathway_c, per_pathway_v)

    # ---- summary metrics ----
    base_metrics <- bind_rows(
      eval_c$overall %>% mutate(method = "Conventional"),
      eval_v$overall %>% mutate(method = "PathwayVote")
    )

    # ---- additional metrics ----
    extra_c <- extract_additional_metrics(conventional_result, sim_data, pathway2gene)
    extra_v <- extract_additional_metrics(vote_result$result_tables$Reactome, sim_data, pathway2gene)

    all_metrics[[i]] <- bind_rows(
      bind_cols(base_metrics[1, ], extra_c, iteration = i),
      bind_cols(base_metrics[2, ], extra_v, iteration = i)
    )
  }

  return(list(
    all_metrics = bind_rows(all_metrics),
    per_pathway_log = bind_rows(per_pathway_log)
  ))
}


cpg2gene_enrichment <- function(ewas_data,
                                top_k = 1000,
                                rank_column = "score",
                                rank_decreasing = TRUE,
                                use_abs = TRUE) {
  suppressMessages({
    library(minfi)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(ReactomePA)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  })

  ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

  if (!rank_column %in% colnames(ewas_data)) {
    stop("Column '", rank_column, "' not found in ewas_data")
  }

  ranking_values <- if (use_abs) abs(ewas_data[[rank_column]]) else ewas_data[[rank_column]]
  top_cpgs <- head(ewas_data$cpg[order(ranking_values, decreasing = rank_decreasing)], top_k)

  top_genes <- ann[top_cpgs, "UCSC_RefGene_Name"]
  gene_symbols <- unique(unlist(strsplit(na.omit(top_genes), ";")))
  gene_symbols <- gene_symbols[gene_symbols != ""]

  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = gene_symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")

  enrich_result <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    readable = FALSE
  )

  return(as.data.frame(enrich_result))
}

simulate_benchmark_dataset <- function(pathway_info,
                                       pathway2gene,
                                       eQTM_df,
                                       signal_genes_per_pathway = 10,
                                       n_noise = 20000,
                                       n_per_gene = 3,
                                       signal_score_mean = 3,
                                       signal_score_sd = 0.5,
                                       noise_score_mean = 0,
                                       noise_score_sd = 1,
                                       correlation_threshold = 0.15,
                                       p_value_threshold = 1e-4,
                                       seed = 1029) {
  set.seed(seed)

  truth_pathway_ids <- list()
  signal_entrez_genes <- c()
  signal_gene_map <- list()

  for (pw in names(pathway_info)) {
    pw_id <- pathway_info[[pw]]
    gene_pool <- pathway2gene[[pw_id]]
    message(sprintf("Pathway %s (%s) contains %s genes.", pw, pw_id, length(gene_pool)))

    if (is.null(gene_pool) || length(gene_pool) < signal_genes_per_pathway) {
      warning(sprintf("Pathway %s (%s) skipped due to insufficient gene count.", pw, pw_id))
      next
    }

    selected_genes <- sample(gene_pool, signal_genes_per_pathway)
    truth_pathway_ids[[pw]] <- pw_id
    signal_entrez_genes <- c(signal_entrez_genes, selected_genes)
    signal_gene_map[[pw]] <- selected_genes
  }

  signal_entrez_genes <- unique(signal_entrez_genes)

  # ---- screen signal CpG ----
  signal_subset <- eQTM_df[
    eQTM_df$entrez %in% signal_entrez_genes &
      eQTM_df$p_value <= p_value_threshold &
      abs(eQTM_df$statistics) >= correlation_threshold, ]

  # ---- each gene keep max n CpG ----
  signal_subset <- do.call(rbind, lapply(split(signal_subset, signal_subset$entrez), function(df) {
    if (nrow(df) <= n_per_gene) return(df)
    return(df[sample(nrow(df), n_per_gene), ])
  }))

  signal_cpgs <- unique(signal_subset$cpg)
  noise_candidates <- setdiff(unique(eQTM_df$cpg), signal_cpgs)
  noise_cpgs <- sample(noise_candidates, min(n_noise, length(noise_candidates)))

  # score
  signal_scores <- abs(rnorm(length(signal_cpgs), mean = signal_score_mean, sd = signal_score_sd))
  noise_scores <- abs(rnorm(length(noise_cpgs), mean = noise_score_mean, sd = noise_score_sd))

  ewas_data <- data.frame(
    cpg = c(signal_cpgs, noise_cpgs),
    score = c(signal_scores, noise_scores),
    is_signal = c(rep(TRUE, length(signal_cpgs)), rep(FALSE, length(noise_cpgs))),
    stringsAsFactors = FALSE
  )

  ewas_data <- ewas_data[sample(nrow(ewas_data)), ]
  rownames(ewas_data) <- NULL

  return(list(
    ewas_data = ewas_data,
    truth_pathway_ids = truth_pathway_ids,
    signal_entrez_genes = signal_entrez_genes,
    signal_gene_map = signal_gene_map
  ))
}

extract_pathways_by_level <- function(hierarchy_table, level, root) {
  current <- root
  all_levels <- list()
  all_levels[[1]] <- current

  for (i in 2:level) {
    next_level <- unique(hierarchy_table$child[hierarchy_table$parent %in% all_levels[[i - 1]]])
    if (length(next_level) == 0) break
    all_levels[[i]] <- next_level
  }

  target_ids <- unique(unlist(all_levels[level]))

  name_df <- as.data.frame(AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME))
  names(target_ids) <- name_df$PATHNAME[match(target_ids, name_df$PATHID)]

  return(target_ids)
}

sample_random_pathways_by_level <- function(hierarchy_table,
                                            n_root = 7,
                                            target_level = 3,
                                            min_depth = 5,
                                            seed = 1029) {
  set.seed(seed)

  # ---- decide how many layers of a pathway ----
  has_n_levels <- function(root_id, depth, hierarchy_table) {
    current <- root_id
    for (i in 2:depth) {
      current <- unique(hierarchy_table$child[hierarchy_table$parent %in% current])
      if (length(current) == 0) return(FALSE)
    }
    return(TRUE)
  }

  # ---- get all top level pathways and filter ----
  all_root_ids <- unique(hierarchy_table$parent[!hierarchy_table$parent %in% hierarchy_table$child])
  root_ids <- Filter(function(pid) has_n_levels(pid, min_depth, hierarchy_table), all_root_ids)

  sampled_roots <- sample(root_ids, n_root)

  name_df <- as.data.frame(AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME))

  pathway_info <- list()

  for (root_id in sampled_roots) {
    current_ids <- extract_pathways_by_level(hierarchy_table, level = target_level, root = root_id)
    selected <- sample(current_ids, 1)
    matched_idx <- match(selected, name_df$DB_ID)

    if (!is.na(matched_idx)) {
      pname <- name_df$path_name[matched_idx]
      if (!is.na(pname) && pname != "") {
        pathway_info[[pname]] <- selected
      }
    }
  }

  return(pathway_info)
}

benchmark_enrichment_results <- function(result_df,
                                         truth_pathway_ids,
                                         hierarchy_table,
                                         fdr_cutoff = 0.1,
                                         top_k = 10,
                                         include_descendants = TRUE,
                                         filter_descendants_by_gene = TRUE,
                                         signal_entrez_genes = NULL,
                                         pathway2gene = NULL,
                                         debug = FALSE) {
  if (!all(c("ID", "p.adjust") %in% colnames(result_df))) {
    stop("result_df must contain 'ID' and 'p.adjust'")
  }
  if (!is.list(truth_pathway_ids)) stop("truth_pathway_ids must be a named list.")
  if (filter_descendants_by_gene && is.null(signal_entrez_genes)) {
    stop("You must provide signal_entrez_genes if filter_descendants_by_gene = TRUE.")
  }
  if (filter_descendants_by_gene && is.null(pathway2gene)) {
    stop("You must provide pathway2gene map if filter_descendants_by_gene = TRUE.")
  }

  # -- Helper: BFS all descendants --
  get_all_descendants <- function(pid) {
    visited <- character()
    queue <- as.character(pid)
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      if (current %in% visited) next
      visited <- c(visited, current)
      children <- hierarchy_table$child[hierarchy_table$parent == current]
      queue <- c(queue, setdiff(children, visited))
    }
    setdiff(visited, pid)
  }

  # -- Signal pathway full ID set --
  signal_pathway_groups <- lapply(truth_pathway_ids, function(pid) {
    pid <- as.character(pid)
    related <- pid
    if (include_descendants) {
      desc <- get_all_descendants(pid)
      if (filter_descendants_by_gene) {
        desc <- desc[desc %in% names(pathway2gene)]
        desc <- desc[vapply(desc, function(p) {
          genes <- pathway2gene[[p]]
          length(intersect(genes, signal_entrez_genes)) > 0
        }, logical(1))]
      }
      related <- c(pid, desc)
    }
    related
  })

  all_signal_ids <- unique(unlist(signal_pathway_groups, use.names = FALSE))

  if (debug) {
    message("All effective signal-related pathway IDs:")
    print(all_signal_ids)
  }

  # -- Sort & annotate --
  result_df <- result_df[!is.na(result_df$p.adjust), , drop = FALSE]
  result_df <- result_df[order(result_df$p.adjust), ]
  result_df$rank <- seq_len(nrow(result_df))
  result_df$is_significant <- result_df$p.adjust <= fdr_cutoff

  # -- Per-pathway TP/FN assignment --
  per_pathway <- purrr::map_df(names(truth_pathway_ids), function(name) {
    pid <- truth_pathway_ids[[name]]
    group <- signal_pathway_groups[[name]]
    detected <- any(result_df$ID %in% group & result_df$is_significant)
    tibble::tibble(
      pathway = name,
      truth_id = paste(pid, collapse = ","),
      detected = detected
    )
  })

  TP <- sum(per_pathway$detected)
  FN <- nrow(per_pathway) - TP

  # -- FP / TN among non-signal pathways --
  non_signal_ids <- setdiff(unique(result_df$ID), all_signal_ids)
  non_signal_df <- result_df[result_df$ID %in% non_signal_ids, , drop = FALSE]
  FP <- sum(non_signal_df$is_significant)
  TN <- sum(!non_signal_df$is_significant)

  summary <- tibble::tibble(
    TP = TP, FN = FN, FP = FP, TN = TN
  )

  return(list(
    per_pathway = per_pathway,
    overall = summary
  ))
}

select_diverse_top_pathways <- function(hierarchy_table, pathway2gene, n = 7, min_genes = 200, max_jaccard = 0.2, seed = 1029) {
  set.seed(seed)

  top_ids <- unique(hierarchy_table$parent[!hierarchy_table$parent %in% hierarchy_table$child])
  top_ids <- top_ids[!grepl("R-HSA-9909396", top_ids)]

  name_df <- as.data.frame(AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME))
  top_names <- name_df$path_name[match(top_ids, name_df$DB_ID)]
  names(top_ids) <- top_names

  get_all_descendants <- function(pid) {
    children <- hierarchy_table$child[hierarchy_table$parent %in% pid]
    if (length(children) == 0) return(pid)
    return(unique(c(pid, unlist(lapply(children, get_all_descendants)))))
  }

  top_gene_sets <- lapply(top_ids, function(pid) {
    all_paths <- get_all_descendants(pid)
    gene_sets <- pathway2gene[intersect(all_paths, names(pathway2gene))]
    unique(unlist(gene_sets))
  })

  valid <- lengths(top_gene_sets) >= min_genes
  top_gene_sets <- top_gene_sets[valid]
  top_ids <- top_ids[valid]

  if (length(top_gene_sets) < n) {
    stop("Not enough top-level pathways with sufficient genes.")
  }

  all_names <- names(top_gene_sets)
  n_pathways <- length(top_gene_sets)
  jaccard_matrix <- matrix(0, nrow = n_pathways, ncol = n_pathways)
  rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- all_names

  for (i in seq_len(n_pathways)) {
    for (j in i:n_pathways) {
      g1 <- top_gene_sets[[i]]
      g2 <- top_gene_sets[[j]]
      intersect_size <- length(intersect(g1, g2))
      union_size <- length(union(g1, g2))
      score <- if (union_size == 0) 0 else intersect_size / union_size
      jaccard_matrix[i, j] <- jaccard_matrix[j, i] <- score
    }
  }

  avg_jaccard <- rowMeans(jaccard_matrix)
  start_pool <- names(sort(avg_jaccard))[1:5]
  selected <- sample(start_pool, 1)
  remaining <- setdiff(all_names, selected)

  while (length(selected) < n && length(remaining) > 0) {
    candidates <- sapply(remaining, function(cand) {
      mean(sapply(selected, function(sel) jaccard_matrix[cand, sel]))
    })

    valid_candidates <- names(candidates)[candidates <= max_jaccard]
    if (length(valid_candidates) == 0) break

    next_pool <- valid_candidates[order(candidates[valid_candidates])[1:min(3, length(valid_candidates))]]
    next_cand <- sample(next_pool, 1)

    selected <- c(selected, next_cand)
    remaining <- setdiff(remaining, next_cand)
  }

  if (length(selected) < n) {
    warning(sprintf("Only %d pathways selected under max_jaccard = %.2f", length(selected), max_jaccard))
  }

  return(top_ids[selected])
}

