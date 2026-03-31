#' RaCInG: R Analysis of Cell Interactions and Graphs
#'
#' `RaCInG` provides functions to prepare ligand-receptor input matrices,
#' generate patient-specific communication graphs, derive kernel or Monte Carlo
#' features, and perform downstream statistical analysis.
#'
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
"_PACKAGE"

utils::globalVariables(
  c(
    "Adj_P_value", "Celltype", "Datasets", "Effect", "Expr_Ligand",
    "Expr_Receptor", "Feature", "Ligand", "Receiver", "Receptor",
    "Sender", "Significant", "category_intercell_source",
    "category_intercell_target"
  )
)

#' Check that optional packages are installed
#'
#' Internal helper used by high-level workflows that depend on optional
#' preprocessing packages.
#'
#' @param packages Character vector of package names.
#'
#' @return Invisibly returns `TRUE` when all packages are available.
#' @keywords internal
.check_installed_packages <- function(packages) {
  missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(
      "The following optional packages are required for this step but are not installed: ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them first or supply precomputed inputs.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Resolve cell-type expression profiles
#'
#' Internal helper that either returns a supplied expression-profile matrix or
#' computes it from an `estimate_expression_profiles()` function available in the
#' current R session.
#'
#' @param counts Gene-by-sample count matrix.
#' @param deconv Deconvolution result matrix.
#' @param cell_expr_profile Optional precomputed cell-type expression profile.
#'
#' @return A data frame of cell-type expression profiles.
#' @keywords internal
.resolve_expression_profiles <- function(counts, deconv, cell_expr_profile = NULL) {
  if (!is.null(cell_expr_profile)) {
    return(cell_expr_profile)
  }

  estimator <- get0("estimate_expression_profiles", mode = "function")
  if (is.null(estimator)) {
    stop(
      "`cell_expr_profile` must be supplied, or a compatible ",
      "`estimate_expression_profiles()` function must be available in the session.",
      call. = FALSE
    )
  }

  expr_counts <- estimator(counts, deconv)
  rbind(sapply(expr_counts, function(x) colMeans(x))) %>% as.data.frame()
}
