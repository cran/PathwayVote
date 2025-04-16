#' @importFrom stats na.omit p.adjust
#' @importFrom utils head
#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly availableCores
#' @importFrom ReactomePA enrichPathway
#' @importFrom clusterProfiler enrichGO enrichKEGG setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @import dplyr

#' @title Pathway Vote Algorithm for eQTM Data (Auto Parallel)
#'
#' @description Performs pathway enrichment analysis using a voting-based approach on eQTM data.
#'
#' @param ewas_data A data.frame with columns: cpg and a ranking column (e.g., p_value, score).
#' @param eQTM An eQTM object containing eQTM data.
#' @param k_values A numeric vector of top k CpGs to select (e.g., c(10, 50, 100)).
#' @param stat_grid A numeric vector of statistics thresholds.
#' @param distance_grid A numeric vector of distance thresholds.
#' @param overlap_threshold A numeric value for gene list overlap threshold.
#' @param databases A character vector of pathway databases (e.g., "Reactome").
#' @param rank_column A character string indicating which column in `ewas_data` to use for ranking.
#' @param rank_decreasing Logical. If TRUE (default), sorts CpGs from high to low based on `rank_column`.
#' @param use_abs Logical. Whether to apply `abs()` to the ranking column before sorting CpGs.
#' @param prune_strategy Character, either "cuberoot" or "fixed". If "cuberoot", the minimum vote support is computed as (N)^(1/3) where N is the number of enrichment combinations. If "fixed", uses the value provided by `fixed_value`.
#' @param fixed_value Integer, used only if `prune_strategy = "fixed"`.
#' @param min_genes_per_hit Minimum number of genes (`Count`) a pathway must include to be considered.
#' @param workers Optional integer. Number of parallel workers. If NULL, use 2 logical cores by default.
#' @param readable Logical. whether to convert Entrez IDs to gene symbols in enrichment results.
#' @param verbose Logical. whether to print progress messages.
#'
#' @return A named list of data.frames, each containing enrichment results (pathway ID, p.adjust, Description, geneID) for one database (e.g., Reactome, KEGG).
#'
#' @examples
#' data <- data.frame(
#'   cpg = c("cg000001", "cg000002", "cg000003"),
#'   statistics = c(2.5, -1.8, 3.2),
#'   p_value = c(0.01, 0.03, 0.005),
#'   distance = c(50000, 80000, 30000),
#'   entrez = c("673", "1956", "5290")
#' )
#' eqtm_obj <- create_eQTM(data)
#' \donttest{
#' results <- pathway_vote(
#'   ewas_data = data,
#'   eQTM = eqtm_obj,
#'   k_values = c(2),
#'   stat_grid = c(1.5),
#'   distance_grid = c(1e5),
#'   overlap_threshold = 0.3,
#'   databases = c("KEGG"),
#'   rank_column = "p_value",
#'   rank_decreasing = FALSE,
#'   use_abs = FALSE,
#'   worker = 1, # If not specified, will use 2 cores by default
#'   verbose = FALSE
#' )
#' }
#'
#' @export
#'
pathway_vote <- function(ewas_data, eQTM, k_values, stat_grid, distance_grid,
                         overlap_threshold, databases = c("Reactome"),
                         rank_column = "p_value",
                         rank_decreasing = FALSE,
                         use_abs = FALSE,
                         prune_strategy = "cuberoot",
                         fixed_value = 3,
                         min_genes_per_hit = 3,
                         readable = FALSE,
                         workers = NULL,
                         verbose = FALSE) {

  # ---- Load and setup parallel environment ----
  required_pkgs <- c("PathwayVote", "furrr", "future", "ReactomePA", "clusterProfiler", "org.Hs.eg.db")
  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Please install it first.")
    }
  })
  suppressMessages({
    lapply(required_pkgs, library, character.only = TRUE)
  })

  # Set workers
  available_cores <- parallelly::availableCores(logical = TRUE)

  if (is.null(workers)) {
    user_specified_workers = FALSE
    # By default, use 2 cores per CRAN's requirement
    workers <- min(workers, 2)
  } else {
    if (!is.numeric(workers) || length(workers) != 1 || workers < 1) {
      stop("`workers` must be a positive integer")
    }
    workers <- min(workers, available_cores)
    user_specified_workers = TRUE
  }

  # Set future plan
  current_plan <- future::plan()
  if (!inherits(current_plan, "multisession") || future::nbrOfWorkers() != workers) {
    future::plan(future::multisession, workers = workers)
  }

  if (verbose) {
    if (workers == 1) {
      message("Using 1 worker (single-thread mode).")
    } else if (!user_specified_workers) {
      message("Using ", workers, " parallel workers (auto-detected).")
    } else {
      message("Using ", workers, " parallel workers (user-specified).")
    }
  }

  # ---- Input checks ----
  if (!inherits(eQTM, "eQTM")) stop("eQTM must be an eQTM object")
  if (!"cpg" %in% colnames(ewas_data)) stop("ewas_data must contain a 'cpg' column.")
  if (!rank_column %in% colnames(ewas_data)) stop("Column '", rank_column, "' not found in ewas_data")
  if (all(is.na(getData(eQTM)$entrez))) stop("Entrez IDs are required for pathway analysis")

  # ---- Sort EWAS data ----
  if (!rank_column %in% colnames(ewas_data)) stop("Column '", rank_column, "' not found in ewas_data")

  ranking_values <- if (use_abs) abs(ewas_data[[rank_column]]) else ewas_data[[rank_column]]
  ewas_data <- ewas_data[order(ranking_values, decreasing = rank_decreasing), ]

  # ---- Filter gene lists for each k ----
  if (verbose) message("Generating gene lists for all k values...")
  all_gene_sets <- list()
  valid_combination_count <- 0  # 统计有效组合数量

  for (k in k_values) {
    if (verbose) message(sprintf("Processing top %d CpGs...", k))
    selected_cpgs <- head(ewas_data$cpg, k)
    eQTM_subset <- new("eQTM", data = getData(eQTM)[getData(eQTM)$cpg %in% selected_cpgs, ],
                       metadata = getMetadata(eQTM))

    filtered_results <- filter_gene_lists(eQTM_subset, stat_grid, distance_grid, overlap_threshold, verbose = FALSE)

    for (i in seq_along(filtered_results$gene_lists)) {
      gene_list_i <- filtered_results$gene_lists[[i]]
      if (length(gene_list_i) == 0) next  # skip invalid combo

      valid_combination_count <- valid_combination_count + 1

      all_gene_sets[[length(all_gene_sets) + 1]] <- list(
        gene_list = gene_list_i,
        param = filtered_results$params[[i]],
        k = k
      )
    }
  }

  if (verbose) message(sprintf("Gene filtering completed. %d valid combinations retained.", valid_combination_count))

  # ---- Run enrichment analysis in parallel ----
  if (verbose) {
    if (workers == 1) {
      message("Running enrichment analysis (single thread)...")
    } else {
      message("Running enrichment analysis (parallel)...")
    }
  }

  enrich_results <- furrr::future_map(
    all_gene_sets,
    function(x) {
      tryCatch({
        run_enrichment(
          gene_list = x$gene_list,
          databases = databases,
          readable = readable,
          verbose = FALSE
        )
      }, error = function(e) {
        warning("Enrichment failed for one gene list: ", conditionMessage(e))
        return(NULL)
      })
    },
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("PathwayVote", "ReactomePA", "clusterProfiler", "org.Hs.eg.db")
    )
  )

  # ---- Prune pathways by vote ----
  enrich_results <- prune_pathways_by_vote(
    enrich_results,
    prune_strategy = prune_strategy,
    fixed_value = fixed_value,
    min_genes = min_genes_per_hit,
    verbose = verbose
  )

  # ---- Combine and return results ----
  result_tables <- combine_enrichment_results(enrich_results, databases, verbose = verbose)
  return(result_tables)
}

