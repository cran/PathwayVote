#' @importFrom methods setClass setGeneric setMethod new validObject
#' @importFrom AnnotationDbi mapIds
#' @import dplyr

#' @title eQTM Class
#'
#' @description A class to store eQTM data for pathway analysis. eQTM stands for Expression Quantitative Trait Methylation.
#' @slot data A data.frame containing eQTM data with columns: cpg, statistics, p_value, distance, and at least one of entrez or ensembl.
#' @slot metadata A list of metadata (e.g., data source, time point). Reserved for future use.
#' @export
#' @name eQTM-class
setClass(
  "eQTM",
  slots = list(
    data = "data.frame",
    metadata = "list"
  ),
  prototype = list(
    data = data.frame(
      cpg = character(),
      gene = character(),
      statistics = numeric(),
      p_value = numeric(),
      distance = numeric(),
      entrez = character(),
      ensembl = character(),
      stringsAsFactors = FALSE
    ),
    metadata = list()
  ),
  validity = function(object) {
    required_cols <- c("cpg", "statistics", "p_value", "distance")
    if (!all(required_cols %in% colnames(object@data))) {
      return("Missing required columns: cpg, statistics, p_value, distance")
    }
    if (!"entrez" %in% colnames(object@data) && !"ensembl" %in% colnames(object@data)) {
      return("At least one of entrez or ensembl must be provided")
    }
    if (!is.character(object@data$cpg)) {
      return("cpg must be a character vector")
    }
    if (!is.numeric(object@data$statistics) || !is.numeric(object@data$p_value) || !is.numeric(object@data$distance)) {
      return("statistics, p_value, and distance must be numeric")
    }
    if (any(object@data$p_value < 0 | object@data$p_value > 1, na.rm = TRUE)) {
      return("p_value must be between 0 and 1")
    }
    if (any(object@data$distance < 0, na.rm = TRUE)) {
      return("distance must be non-negative")
    }
    TRUE
  }
)

#' @title Create an eQTM object
#'
#' @param data A data.frame containing eQTM data with columns: cpg, statistics, p_value, distance, and at least one of entrez or ensembl.
#' \describe{
#'   \item{cpg}{Character. CpG probe ID (e.g., "cg00000029"), representing a methylation site.}
#'   \item{statistics}{Numeric. Test statistic from eQTM association analysis (e.g., correlation coefficient, r-square, regression coefficient, or t-statistic). Can be positive or negative.}
#'   \item{p_value}{Numeric. P-value associated with the test statistic, must be between 0 and 1.}
#'   \item{distance}{Numeric. Genomic distance (in base pairs) between the CpG and the associated gene's transcription start site (TSS). Must be non-negative.}
#'   \item{entrez}{Character. Entrez gene ID of the associated gene. At least one of entrez or ensembl must be provided.}
#'   \item{ensembl}{Character. Ensembl gene ID of the associated gene. At least one of entrez or ensembl must be provided.}
#' }
#' @param metadata A list of metadata (optional).
#'
#' @return An eQTM object.
#'
#' @examples
#' data <- data.frame(
#'   cpg = c("cg000001", "cg000002"),
#'   statistics = c(2.5, -1.8),
#'   p_value = c(0.01, 0.03),
#'   distance = c(50000, 80000),
#'   entrez = c("673", "1956")
#' )
#' eqtm_obj <- create_eQTM(data)
#'
#' @export
#'
create_eQTM <- function(data, metadata = list()) {
  if (!is.data.frame(data)) {
    stop("Input data must be a data.frame")
  }

  required_cols <- c("cpg", "statistics", "p_value", "distance")
  if (!all(required_cols %in% colnames(data))) {
    stop("Missing required columns: cpg, statistics, p_value, distance")
  }

  if (!"gene" %in% colnames(data)) data$gene <- NA_character_
  if (!"entrez" %in% colnames(data)) data$entrez <- NA_character_
  if (!"ensembl" %in% colnames(data)) data$ensembl <- NA_character_

  if (all(is.na(data$entrez)) && all(is.na(data$ensembl))) {
    stop("At least one of entrez or ensembl must be provided")
  }

  data$cpg <- as.character(data$cpg)
  data$gene <- as.character(data$gene)
  data$entrez <- as.character(data$entrez)
  data$ensembl <- as.character(data$ensembl)

  if (all(is.na(data$entrez)) && !all(is.na(data$ensembl))) {
    message("Converting Ensembl IDs to Entrez IDs and gene symbols...")

    suppressWarnings({
      ensembl_to_entrez <- mapIds(org.Hs.eg.db, keys = data$ensembl, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
      ensembl_to_symbol <- mapIds(org.Hs.eg.db, keys = data$ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    })

    data$entrez <- ensembl_to_entrez[match(data$ensembl, names(ensembl_to_entrez))]
    data$gene <- ensembl_to_symbol[match(data$ensembl, names(ensembl_to_symbol))]

    n_missing <- sum(is.na(data$entrez))
    if (n_missing > 0) {
      warning(sprintf("Entrez ID conversion failed for %d Ensembl IDs (%.2f%%).", n_missing, 100 * n_missing / nrow(data)))
    }

  } else if (all(is.na(data$ensembl)) && !all(is.na(data$entrez))) {
    message("Converting Entrez IDs to Ensembl IDs and gene symbols...")

    suppressWarnings({
      entrez_to_ensembl <- mapIds(org.Hs.eg.db, keys = data$entrez, column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first")
      entrez_to_symbol <- mapIds(org.Hs.eg.db, keys = data$entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    })

    data$ensembl <- entrez_to_ensembl[match(data$entrez, names(entrez_to_ensembl))]
    data$gene <- entrez_to_symbol[match(data$entrez, names(entrez_to_symbol))]

    n_missing <- sum(is.na(data$ensembl))
    if (n_missing > 0) {
      warning(sprintf("Ensembl ID conversion failed for %d Entrez IDs (%.2f%%).", n_missing, 100 * n_missing / nrow(data)))
    }
  } else {
    # Both provided: do nothing
    message("Using provided Ensembl and Entrez IDs without conversion.")
  }

  eQTM_obj <- new("eQTM", data = data, metadata = metadata)
  validObject(eQTM_obj)
  return(eQTM_obj)
}

#' @title Get eQTM Data
#' @description Retrieve the eQTM data.frame from an eQTM object.
#' @param object An eQTM object.
#' @return A data.frame stored in the object.
#' @name getData
#' @export
NULL

#' @title Get eQTM Metadata
#' @description Retrieve the metadata list from an eQTM object.
#' @param object An eQTM object.
#' @return A list containing metadata.
#' @name getMetadata
#' @export
NULL

#' @rdname getData
#' @export
setGeneric("getData", function(object) standardGeneric("getData"))

#' @rdname getData
#' @export
setMethod("getData", "eQTM", function(object) object@data)

#' @rdname getMetadata
#' @export
setGeneric("getMetadata", function(object) standardGeneric("getMetadata"))

#' @rdname getMetadata
#' @export
setMethod("getMetadata", "eQTM", function(object) object@metadata)

