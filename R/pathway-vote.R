#' @importFrom stats glm predict quasipoisson setNames na.omit p.adjust quantile
#' @importFrom harmonicmeanp p.hmp
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
#' CpG-gene regulatory information from expression quantitative trait methylation (eQTM) data.
#' For a grid of top-ranked CpGs and filtering thresholds, gene sets are generated and refined using
#' an entropy-based pruning strategy that balances information richness, stability, and probe bias correction.
#' In particular, gene lists dominated by genes with disproportionately high numbers of CpG mappings
#' are penalized to mitigate active probe bias, a common artifact in methylation data analysis.
#' Enrichment results across parameter combinations are then aggregated using a voting scheme,
#' prioritizing pathways that are consistently recovered under diverse settings and robust to parameter perturbations.
#'
#' @param cpg_input A data.frame containing CpG-level results or identifiers. The first column must contain CpG IDs,
#' which can be Illumina probe IDs (e.g., "cg00000029") for array-based data, or genomic coordinates
#' (e.g., "chr1:10468" or "chr1:10468:+") for sequencing-based data. These IDs will be matched against
#' the eQTM object. Optionally, a second column may provide a ranking metric. If supplied, this must be:
#' (i) the complete set of raw p-values from association tests (required for automatic \code{k_grid}
#' generation), or (ii) an alternative metric such as t-statistics or feature importance scores, in which
#' case \code{k_grid} must be specified manually. If no ranking information is provided, all input CpGs
#' are used directly and \code{k_grid} is ignored.
#' @param eQTM An \code{eQTM} object containing CpG-gene linkage information, created by the \code{create_eQTM()} function. This object provides
#' the CpG-to-gene mapping used for pathway inference. Please make sure the CpG IDs used here match those in \code{cpg_input}.
#' @param databases A character vector of pathway databases. Supporting: "Reactome", "KEGG", and "GO".
#' @param k_grid A numeric vector specifying the top-k CpGs used for gene set construction. If \code{NULL}, the grid is inferred automatically,
#' but this requires that \code{cpg_input} contains: (i) the complete set of CpGs tested (first column), and (ii) raw p-values from the association
#' test (second column). If these conditions are not satisfied, or if alternative ranking metrics are provided (e.g., t-statistics, feature importance scores),
#' then \code{k_grid} must be specified manually.
#' @param stat_grid A numeric vector of eQTM statistic thresholds. If NULL, generated based on quantiles of the observed distribution.
#' @param distance_grid A numeric vector of CpG-gene distance thresholds (in base pairs). If NULL, generated based on quantiles of the observed distribution.
#' @param grid_size Integer. Number of values in each grid when auto-generating. Default is 5.
#' @param overlap_threshold Numeric between 0 and 1. Controls the maximum allowed Jaccard similarity between gene lists during redundancy filtering.
#' Default is 0.7, which provides robust and stable results across a variety of simulation scenarios.
#' @param fixed_prune Integer or NULL. Minimum number of votes to retain a pathway. If NULL, will use cuberoot(N) where N is the number of total enrichment runs.
#' @param min_genes_per_hit Minimum number of genes a pathway must include to be considered. Default is 2.
#' @param workers Optional integer. Number of parallel workers. If NULL, use 2 logical cores.
#' @param readable Logical. Whether to convert Entrez IDs to gene symbols in enrichment results.
#' @param verbose Logical. Whether to print progress messages.
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
#'   cpg_input = ewas,
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
pathway_vote <- function(cpg_input, eQTM,
                         databases = c("Reactome"),
                         k_grid = NULL,
                         stat_grid = NULL,
                         distance_grid = NULL,
                         grid_size = 5,
                         overlap_threshold = 0.7,
                         fixed_prune = NULL,
                         min_genes_per_hit = 2,
                         readable = FALSE,
                         workers = NULL,
                         verbose = FALSE) {

  if (is.vector(cpg_input) && is.character(cpg_input)) {
    cpg_input <- data.frame(cpg = cpg_input, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(cpg_input)) {
    stop("`cpg_input` must be a data.frame or a character vector of CpG IDs.")
  }
  if (!"cpg" %in% colnames(cpg_input)) {
    colnames(cpg_input)[1] <- "cpg"
  }
  cpg_input$cpg <- as.character(cpg_input$cpg)
  cpg_input <- cpg_input[!is.na(cpg_input$cpg) & nzchar(cpg_input$cpg), , drop = FALSE]
  cpg_input <- cpg_input[!duplicated(cpg_input$cpg), , drop = FALSE]

  if (nrow(cpg_input) == 0) {
    stop("No valid CpG IDs remain after filtering NA/empty/duplicates in `cpg_input`.")
  }

  if (verbose) {
    message("==== PathwayVote Start ====")
    message("Input rows: ", nrow(cpg_input))
    message("Databases: ", paste(databases, collapse = ", "))
  }

  # required_pkgs <- c("PathwayVote", "purrr", "furrr", "future", "ReactomePA", "clusterProfiler", "org.Hs.eg.db")
  # lapply(required_pkgs, function(pkg) {
  #   if (!requireNamespace(pkg, quietly = TRUE)) {
  #     stop("Package '", pkg, "' is required. Please install it first.")
  #   }
  # })
  # suppressMessages({
  #   lapply(required_pkgs, library, character.only = TRUE)
  # })

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
  if (verbose) message("Using ", workers, ifelse(workers == 1, " worker", " workers"))

  if (!inherits(eQTM, "eQTM")) stop("eQTM must be an eQTM object")
  if (!check_cpg_match(cpg_input, eQTM, verbose)) {
    stop("First column of `cpg_input` does not match CpG IDs in eQTM object. Please verify your input.")
  }
  if (all(is.na(getData(eQTM)$entrez))) stop("Entrez IDs are required for pathway analysis")

  # ------------------------------
  # Detect mode: CpG-only vs Ranking
  # ------------------------------
  mode_cpg_only <- TRUE
  rank_column <- NULL
  rank_values <- NULL

  if (ncol(cpg_input) >= 2) {
    rank_column <- colnames(cpg_input)[2]
    rank_values <- cpg_input[[rank_column]]
    if (is.numeric(rank_values)) {
      mode_cpg_only <- FALSE
    }
  }

  if (mode_cpg_only) {
    if (verbose) message("Mode detected: CpG-only (no ranking). `k_grid` will be ignored.")
  } else {
    # Traditional ranking mode checks
    is_pval_like <- all(rank_values >= 0 & rank_values <= 1, na.rm = TRUE)
    use_abs <- !is_pval_like
    rank_decreasing <- !is_pval_like

    if (verbose) {
      if (is_pval_like) {
        message(sprintf("Mode detected: Ranking mode using column `%s` as raw p-values.", rank_column))
        message("  Ranking rule: smaller p-values rank higher (ascending order).")
      } else {
        message(sprintf("Mode detected: Ranking mode using column `%s` as numeric scores.", rank_column))
        message(sprintf("  Ranking rule: larger absolute values of `%s` rank higher (descending order).", rank_column))
      }
    }

    if (is.null(k_grid)) {
      if (!is_pval_like) {
        stop("Automatic k_grid generation is only supported when the second column is p-values. ",
             "Please provide k_grid manually for other ranking metrics, or provide CpG-only input.")
      }
      k_grid <- generate_k_grid_fdr_guided(cpg_input, rank_column, grid_size, verbose = verbose)
    }

    ranking_values <- if (use_abs) abs(rank_values) else rank_values
    cpg_input <- cpg_input[order(ranking_values, decreasing = rank_decreasing), , drop = FALSE]
  }

  # ------------------------
  # Grids for stat/distance
  # ------------------------
  if (is.null(stat_grid)) {
    stat_vals <- abs(getData(eQTM)$statistics)
    range_stat <- stats::quantile(stat_vals, probs = c(0.05, 0.95), na.rm = TRUE)
    stat_grid <- round(seq(range_stat[1], range_stat[2], length.out = grid_size), 2)
    if (verbose) message("Auto-selected stat_grid: ", paste(stat_grid, collapse = ", "))
  }
  if (is.null(distance_grid)) {
    dist_vals <- getData(eQTM)$distance
    range_dist <- stats::quantile(dist_vals, probs = c(0.05, 0.95), na.rm = TRUE)
    distance_grid <- round(seq(range_dist[1], range_dist[2], length.out = grid_size), -3)
    if (verbose) message("Auto-selected distance_grid: ", paste(distance_grid, collapse = ", "))
  }

  cpg_count <- get_cpg_count_per_gene(eQTM)

  # ------------------------------
  # Generate candidate gene lists
  # ------------------------------
  if (verbose) message("==== Generating gene lists ====")

  total_initial <- 0L
  total_retained <- 0L
  all_gene_sets <- list()

  if (mode_cpg_only) {
    # Single pass using all provided CpGs
    selected_cpgs <- unique(cpg_input$cpg)
    if (verbose) message("CpG-only: using all provided CpGs (n = ", length(selected_cpgs), ") in a single pass.")
    eQTM_subset <- new("eQTM", data = getData(eQTM)[getData(eQTM)$cpg %in% selected_cpgs, ],
                       metadata = getMetadata(eQTM))

    if (nrow(getData(eQTM_subset)) == 0)
      stop("None of the input CpG IDs matched the eQTM object. Check ID formats (probe IDs vs genomic coordinates).")

    raw_results <- generate_gene_lists_grid(eQTM_subset, stat_grid, distance_grid)
    initial_k <- length(raw_results$gene_lists)
    total_initial <- total_initial + initial_k

    entropy_filtered_lists <- select_gene_lists_entropy_auto(
      gene_lists = raw_results$gene_lists,
      cpg_count = cpg_count,
      overlap_threshold = overlap_threshold
    )

    kept_indices <- which(vapply(raw_results$gene_lists, function(x)
      any(sapply(entropy_filtered_lists, function(y) setequal(x, y))), logical(1)))

    retained_k <- length(kept_indices)
    total_retained <- total_retained + retained_k
    if (verbose) {
      pct_k <- if (initial_k > 0) 100 * retained_k / initial_k else 0
      message(sprintf("CpG-only: retained %d of %d candidates (%.1f%%).", retained_k, initial_k, pct_k))
    }

    for (i in kept_indices) {
      all_gene_sets[[length(all_gene_sets) + 1]] <- list(
        gene_list = raw_results$gene_lists[[i]],
        param = raw_results$params[[i]],
        k = NA_integer_
      )
    }

  } else {
    # loop over k_grid
    for (k in k_grid) {
      if (verbose) message(sprintf("Processing top %d CpGs...", k))
      selected_cpgs <- utils::head(cpg_input$cpg, k)
      eQTM_subset <- new("eQTM", data = getData(eQTM)[getData(eQTM)$cpg %in% selected_cpgs, ],
                         metadata = getMetadata(eQTM))

      if (nrow(getData(eQTM_subset)) == 0) {
        if (verbose) message("No CpG matched eQTM under current top-k selection; skipping this k.")
        next
      }

      raw_results <- generate_gene_lists_grid(eQTM_subset, stat_grid, distance_grid)
      initial_k <- length(raw_results$gene_lists)
      total_initial <- total_initial + initial_k

      entropy_filtered_lists <- select_gene_lists_entropy_auto(
        gene_lists = raw_results$gene_lists,
        cpg_count = cpg_count,
        overlap_threshold = overlap_threshold
      )

      kept_indices <- which(vapply(raw_results$gene_lists, function(x)
        any(sapply(entropy_filtered_lists, function(y) setequal(x, y))), logical(1)))

      retained_k <- length(kept_indices)
      total_retained <- total_retained + retained_k
      if (verbose) {
        pct_k <- if (initial_k > 0) 100 * retained_k / initial_k else 0
        message(sprintf("k = %d: retained %d of %d candidates (%.1f%%).",
                        k, retained_k, initial_k, pct_k))
      }

      for (i in kept_indices) {
        all_gene_sets[[length(all_gene_sets) + 1]] <- list(
          gene_list = raw_results$gene_lists[[i]],
          param = raw_results$params[[i]],
          k = k
        )
      }
    }
  }

  if (verbose) {
    pct_all <- if (total_initial > 0) 100 * total_retained / total_initial else 0
    message(sprintf(
      "Gene list candidates filtering completed. %d valid candidates out of %d initial candidates retained (%.1f%%).",
      total_retained, total_initial, pct_all
    ))
  }

  if (length(all_gene_sets) == 0) {
    stop("No candidate gene lists were retained after entropy pruning. ",
         "Try relaxing thresholds (e.g., lower stat_grid / increase distance_grid) ",
         "or check CpG-eQTM mapping coverage.")
  }

  gene_lists <- lapply(all_gene_sets, function(x) x$gene_list)
  universe_genes <- unique(stats::na.omit(getData(eQTM)$entrez))

  if (verbose) message("==== Preparing all annotation data ====")
  preloaded_data <- prepare_enrichment_data(databases, verbose = verbose)

  if (verbose) message("==== Running enrichment analysis ====")
  enrich_results <- furrr::future_map(
    gene_lists,
    function(glist) {
      tryCatch({
        run_enrichment(
          gene_list = glist,
          databases = databases,
          universe = universe_genes,
          readable = readable,
          preloaded_data = preloaded_data
        )
      }, error = function(e) {
        warning("Enrichment failed for gene list ", paste(utils::head(glist, 5), collapse = ","), ": ", conditionMessage(e))
        NULL
      })
    },
    .options = furrr::furrr_options(seed = TRUE),
    .progress = verbose
  )

  enrich_results <- prune_pathways_by_vote(
    enrich_results,
    fixed_prune = fixed_prune,
    min_genes = min_genes_per_hit,
    verbose = verbose
  )

  result_tables <- combine_enrichment_results(enrich_results, databases, verbose = verbose)
  return(result_tables)
}
