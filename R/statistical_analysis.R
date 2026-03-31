#' Run Wilcoxon tests across network features
#'
#' @param data_matrix Numeric matrix or data frame with patients in rows and features in columns.
#' @param groups Vector of group labels with length matching `nrow(data_matrix)`.
#' @param p_adjust_method Multiple-testing correction method passed to [stats::p.adjust()].
#'
#' @return A data frame with test statistics and adjusted p-values.
#' @export
wilcox_group_test <- function(data_matrix, groups, p_adjust_method = "fdr") {
  # data_matrix: numeric matrix or data.frame, rows = patients, cols = features
  # groups: factor or character vector of group labels (length = nrow(data_matrix))
  # p_adjust_method: method for multiple testing correction ("fdr", "bonferroni", etc.)

  if(length(groups) != nrow(data_matrix)) {
    stop("Length of groups must match the number of rows in data_matrix")
  }

  groups <- as.factor(groups)

  # Initialize vectors to store results
  features <- colnames(data_matrix)
  p_values <- numeric(length(features))
  stats <- numeric(length(features))

  for(i in seq_along(features)) {
    feature_values <- data_matrix[, i]

    # Wilcoxon test between two groups
    test_res <- stats::wilcox.test(feature_values ~ groups)

    p_values[i] <- test_res$p.value
    stats[i] <- test_res$statistic
  }

  # Adjust p-values for multiple testing
  p_adj <- stats::p.adjust(p_values, method = p_adjust_method)

  # Return results as a data.frame
  result_df <- data.frame(
    Feature = features,
    Wilcox_statistic = stats,
    P_value = p_values,
    Adjusted_P_value = p_adj,
    stringsAsFactors = FALSE
  )

  # Sort by adjusted p-value
  result_df <- result_df[order(result_df$Adjusted_P_value), ]
  return(result_df)
}

#' Create a volcano plot from Wilcoxon results
#'
#' @param wilcox_results Output of [wilcox_group_test()].
#' @param top_labels Number of top significant features to label.
#' @param p_threshold Adjusted p-value threshold used to mark significance.
#' @param title Plot title.
#'
#' @return A `ggplot2` object.
#' @export
volcano_plot <- function(wilcox_results, top_labels = 10,
                                     p_threshold = 0.05, title = "Volcano Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for `volcano_plot()`. Please install it first.", call. = FALSE)
  }

  # Make sure the required columns exist
  if(!all(c("Feature", "Adjusted_P_value", "Wilcox_statistic") %in% colnames(wilcox_results))) {
    stop("wilcox_results must contain columns: Feature, Wilcox_statistic, Adjusted_P_value")
  }

  # Build plot dataframe
  plot_df <- data.frame(
    Feature = wilcox_results$Feature,
    Effect = wilcox_results$Wilcox_statistic,   # or replace with your real effect size if available
    Adj_P_value = wilcox_results$Adjusted_P_value
  )

  # Significant features
  plot_df$Significant <- ifelse(plot_df$Adj_P_value < p_threshold, "Yes", "No")

  # Identify top features by significance
  top_features <- utils::head(plot_df[plot_df$Significant == "Yes", ][order(plot_df$Adj_P_value), "Feature"], top_labels)

  # Plot
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Effect, y = -log10(Adj_P_value))) +
    ggplot2::geom_point(ggplot2::aes(color = Significant), shape = 4, alpha = 0.7) +
    ggplot2::geom_text(
      data = subset(plot_df, Feature %in% top_features),
      ggplot2::aes(label = Feature),
      size = 3,
      vjust = -0.5,
      check_overlap = TRUE
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), color = "red", linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("No" = "gainsboro", "Yes" = "blue")) +
    ggplot2::xlim(min(plot_df$Effect) - 1, max(plot_df$Effect) + 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Wilcox Statistic / Effect Size",
      y = "-log10(Adjusted P-value)",
      title = title,
      color = "Significant"
    ) +
    ggplot2::theme(legend.position = "top")
}
