#' @importFrom stats na.omit p.adjust quantile
#' @importFrom utils head
#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly availableCores
#' @importFrom AnnotationDbi select keys
#' @importFrom reactome.db reactome.db
#' @importFrom clusterProfiler enricher download_KEGG setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db

#' @title Pathway Voting-Based Enrichment Analysis
#'
#' @description
#' Performs pathway enrichment analysis using a voting-based framework that integrates
#' CpG–gene regulatory information from expression quantitative trait methylation (eQTM) data.
#' For a grid of top-ranked CpGs and filtering thresholds, gene sets are generated and refined using
#' an entropy-based pruning strategy that balances information richness, stability, and probe bias correction.
#' In particular, gene lists dominated by genes with disproportionately high numbers of CpG mappings
#' are penalized to mitigate active probe bias—a common artifact in methylation data analysis.
#' Enrichment results across parameter combinations are then aggregated using a voting scheme,
#' prioritizing pathways that are consistently recovered under diverse settings and robust to parameter perturbations.
#'
#' @param ewas_data A data.frame containing CpG-level association results. The first column must contain CpG probe IDs,
#' which will be matched against the eQTM object. The second column should contain a numeric ranking metric, such as a p-value, t-statistic, or
#' feature importance score.
#' @param eQTM An \code{eQTM} object containing CpG–gene linkage information, created by the \code{create_eQTM()} function. This object provides
#' the CpG-to-gene mapping used for pathway inference.
#' @param databases A character vector of pathway databases. Supporting: "Reactome", "KEGG", and "GO".
#' @param k_grid A numeric vector of top-k CpGs used for gene set construction. If NULL, the grid is automatically inferred using a log-scaled range
#' guided by the number of CpGs passing FDR < 0.05. Note: This requires that \code{ewas_data} contains raw p-values (second column) from an association analysis;
#' for other metrics (e.g., t-statistic or importance scores), \code{k_grid} must be provided manually.
#' @param stat_grid A numeric vector of eQTM statistic thresholds. If NULL, generated based on quantiles of the observed distribution.
#' @param distance_grid A numeric vector of CpG-gene distance thresholds (in base pairs). If NULL, generated similarly.
#' @param fixed_prune Integer or NULL. Minimum number of votes to retain a pathway. If NULL, will use cuberoot(N) where N is the number of enrichment runs.
#' @param grid_size Integer. Number of values in each grid when auto-generating. Default is 5.
#' @param min_genes_per_hit Minimum number of genes (`Count`) a pathway must include to be considered. Default is 2.
#' @param overlap_threshold Numeric between 0 and 1. Controls the maximum allowed Jaccard similarity between gene lists during redundancy filtering.
#' Defualt is 0.7, which provides robust and stable results across a variety of simulation scenarios.
#' @param workers Optional integer. Number of parallel workers. If NULL, use 2 logical cores.
#' @param readable Logical. whether to convert Entrez IDs to gene symbols in enrichment results.
#' @param verbose Logical. whether to print progress messages.
#'
#' @return A named list of data.frames, each corresponding to a selected pathway database
#' (e.g., `Reactome`, `KEGG`, `GO`). Each data.frame contains enriched pathways with
#' columns: `ID`, `p.adjust`, `Description`, and `geneID`.
#'
#' @examples
#' set.seed(123)
#'
#' # Simulated EWAS result: a mix of signal and noise
#' n_cpg <- 500
#' ewas <- data.frame(
#'   cpg = paste0("cg", sprintf("%08d", 1:n_cpg)),
#'   p_value = c(runif(n_cpg*0.1, 1e-9, 1e-5), runif(n_cpg*0.2, 1e-3, 0.05), runif(n_cpg*0.7, 0.05, 1))
#' )
#'
#' # Corresponding eQTM mapping (some of these CpGs have gene links)
#' signal_genes <- c("5290", "673", "1956", "7157", "7422")
#' background_genes <- as.character(1000:9999)
#' entrez_signal <- sample(signal_genes, n_cpg * 0.1, replace = TRUE)
#' entrez_background <- sample(setdiff(background_genes, signal_genes), n_cpg * 0.9, replace = TRUE)
#'
#' eqtm_data <- data.frame(
#'   cpg = ewas$cpg,
#'   statistics = rnorm(n_cpg, mean = 2, sd = 1),
#'   p_value = runif(n_cpg, min = 0.001, max = 0.05),
#'   distance = sample(1000:100000, n_cpg, replace = TRUE),
#'   entrez = c(entrez_signal, entrez_background),
#'   stringsAsFactors = FALSE
#' )
#' eqtm_obj <- create_eQTM(eqtm_data)
#'
#' # Run pathway voting with minimal settings
#' \dontrun{
#' results <- pathway_vote(
#'   ewas_data = ewas,
#'   eQTM = eqtm_obj,
#'   databases = c("GO", "KEGG", "Reactome"),
#'   readable = TRUE,
#'   verbose = TRUE
#' )
#' head(results$GO)
#' head(results$KEGG)
#' head(results$Reactome)
#' }
#'
#' @export
#'
pathway_vote <- function(ewas_data, eQTM,
                         databases = c("Reactome"),
                         k_grid = NULL,
                         stat_grid = NULL,
                         distance_grid = NULL,
                         fixed_prune = NULL,
                         grid_size = 5,
                         min_genes_per_hit = 2,
                         overlap_threshold = 0.7,
                         readable = FALSE,
                         workers = NULL,
                         verbose = FALSE) {

  if (verbose) {
    message("==== PathwayVote Start ====")
    message("Input CpGs: ", nrow(ewas_data))
    message("Databases: ", paste(databases, collapse = ", "))
  }

  required_pkgs <- c("PathwayVote", "purrr", "furrr", "future", "ReactomePA", "clusterProfiler", "org.Hs.eg.db")
  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Please install it first.")
    }
  })
  suppressMessages({
    lapply(required_pkgs, library, character.only = TRUE)
  })

  available_cores <- parallelly::availableCores(logical = TRUE)
  user_specified_workers <- !is.null(workers)

  if (is.null(workers)) {
    workers <- min(2, available_cores)
  } else {
    if (!is.numeric(workers) || length(workers) != 1 || workers < 1) {
      stop("`workers` must be a positive integer")
    }
    workers <- min(workers, available_cores)
  }

  safe_setup_plan(workers)

  if (verbose) {
    message("Using ", workers, ifelse(workers == 1, " worker", " workers"))
  }

  if (!inherits(eQTM, "eQTM")) stop("eQTM must be an eQTM object")
  if (!check_ewas_cpg_match(ewas_data, eQTM)) stop("First column of `ewas_data` does not match CpG IDs in eQTM object. Please verify your input.")
  if (all(is.na(getData(eQTM)$entrez))) stop("Entrez IDs are required for pathway analysis")
  if (ncol(ewas_data) < 2) stop("ewas_data must have at least two columns: CpG ID (first column) and a ranking column")

  # Automatically use the second column as ranking
  rank_column <- colnames(ewas_data)[2]
  rank_values <- ewas_data[[rank_column]]

  if (!is.numeric(rank_values)) {
    stop("The second column of ewas_data must be numeric for ranking purposes.")
  }

  is_pval_like <- all(rank_values >= 0 & rank_values <= 1, na.rm = TRUE)
  use_abs <- !is_pval_like
  rank_decreasing <- !is_pval_like

  if (verbose) {
    message(sprintf("Auto-selected ranking: %s (decreasing = %s, absolute = %s)",
                    rank_column, rank_decreasing, use_abs))
  }

  if (is.null(k_grid)) {
    if (!is_pval_like) {
      stop("Automatic k_grid generation is only supported when ewas_data contains p-values (second column). ",
           "Please provide k_grid manually for other ranking metrics.")
    }
    k_grid <- generate_k_grid_fdr_guided(ewas_data, rank_column, grid_size, verbose = verbose)
  }

  if (is.null(stat_grid)) {
    stat_vals <- abs(getData(eQTM)$statistics)
    range_stat <- quantile(stat_vals, probs = c(0.05, 0.95), na.rm = TRUE)
    stat_grid <- round(seq(range_stat[1], range_stat[2], length.out = grid_size), 2)
    if (verbose) message("Auto-selected stat_grid: ", paste(stat_grid, collapse = ", "))
  }

  if (is.null(distance_grid)) {
    dist_vals <- getData(eQTM)$distance
    range_dist <- quantile(dist_vals, probs = c(0.05, 0.95), na.rm = TRUE)
    distance_grid <- round(seq(range_dist[1], range_dist[2], length.out = grid_size), -3)
    if (verbose) message("Auto-selected distance_grid: ", paste(distance_grid, collapse = ", "))
  }

  ranking_values <- if (use_abs) abs(rank_values) else rank_values
  ewas_data <- ewas_data[order(ranking_values, decreasing = rank_decreasing), ]

  if (verbose) message("Generating gene lists for all k values...")
  all_gene_sets <- list()
  valid_combination_count <- 0

  for (k in k_grid) {
    if (verbose) message(sprintf("Processing top %d CpGs...", k))
    selected_cpgs <- head(ewas_data$cpg, k)
    eQTM_subset <- new("eQTM", data = getData(eQTM)[getData(eQTM)$cpg %in% selected_cpgs, ],
                       metadata = getMetadata(eQTM))

    # Step 1: Generate all candidate gene lists
    raw_results <- generate_gene_lists_grid(eQTM_subset, stat_grid, distance_grid, verbose = verbose)

    # Step 2: Apply entropy + stability-based + probe bias correction pruning
    entropy_filtered_lists <- select_gene_lists_entropy_auto(
      gene_lists = raw_results$gene_lists,
      grid_size = grid_size,
      overlap_threshold = overlap_threshold,
      verbose = verbose
    )

    # Step 3: Match back to their corresponding stat/distance params
    kept_indices <- which(vapply(raw_results$gene_lists, function(x)
      any(sapply(entropy_filtered_lists, function(y) setequal(x, y))), logical(1)))

    for (i in kept_indices) {
      gene_list_i <- raw_results$gene_lists[[i]]
      all_gene_sets[[length(all_gene_sets) + 1]] <- list(
        gene_list = gene_list_i,
        param = raw_results$params[[i]],
        k = k
      )
      valid_combination_count <- valid_combination_count + 1
    }
  }

  if (verbose) message(sprintf("Gene filtering completed. %d valid combinations retained.", valid_combination_count))

  gene_lists <- lapply(all_gene_sets, function(x) x$gene_list)
  universe_genes <- unique(na.omit(getData(eQTM)$entrez))

  if (verbose) message("Preparing all annotation data...")
  preloaded_data <- prepare_enrichment_data(databases, verbose = verbose)

  if (verbose) message("Running enrichment analysis...")
  enrich_results <- furrr::future_map(
    gene_lists,
    function(glist) {
      tryCatch({
        # 2. Use the fast function with preloaded data
        run_enrichment(
          gene_list = glist,
          databases = databases,
          universe = universe_genes,
          readable = readable,
          preloaded_data = preloaded_data
        )
      }, error = function(e) {
        warning("Enrichment failed for gene list ", paste(head(glist, 5), collapse = ","), ": ", conditionMessage(e))
        NULL
      })
    },
    .options = furrr::furrr_options(seed = TRUE),
    .progress = verbose
  )

  # Will use cube-root(N) if fixed_prune is NULL
  enrich_results <- prune_pathways_by_vote(
    enrich_results,
    fixed_prune = fixed_prune,
    min_genes = min_genes_per_hit,
    verbose = verbose
  )

  result_tables <- combine_enrichment_results(enrich_results, databases, verbose = verbose)
  return(result_tables)
}
