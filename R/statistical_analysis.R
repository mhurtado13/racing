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
    test_res <- wilcox.test(feature_values ~ groups)

    p_values[i] <- test_res$p.value
    stats[i] <- test_res$statistic
  }

  # Adjust p-values for multiple testing
  p_adj <- p.adjust(p_values, method = p_adjust_method)

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

volcano_plot <- function(wilcox_results, top_labels = 10,
                                     p_threshold = 0.05, title = "Volcano Plot") {
  library(ggplot2)

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
  top_features <- head(plot_df[plot_df$Significant == "Yes", ][order(plot_df$Adj_P_value), "Feature"], top_labels)

  # Plot
  ggplot(plot_df, aes(x = Effect, y = -log10(Adj_P_value))) +
    geom_point(aes(color = Significant), shape = 4, alpha = 0.7) +
    geom_text(data = subset(plot_df, Feature %in% top_features),
              aes(label = Feature), size = 3, vjust = -0.5, check_overlap = TRUE) +
    geom_hline(yintercept = -log10(p_threshold), color = "red", linetype = "dashed") +
    scale_color_manual(values = c("No" = "gainsboro", "Yes" = "blue")) +
    xlim(min(plot_df$Effect) - 1, max(plot_df$Effect) + 1) +
    theme_minimal() +
    labs(x = "Wilcox Statistic / Effect Size",
         y = "-log10(Adjusted P-value)",
         title = title,
         color = "Significant") +
    theme(legend.position = "top")
}
