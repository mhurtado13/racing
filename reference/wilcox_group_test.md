# Run Wilcoxon tests across network features

Run Wilcoxon tests across network features

## Usage

``` r
wilcox_group_test(data_matrix, groups, p_adjust_method = "fdr")
```

## Arguments

- data_matrix:

  Numeric matrix or data frame with patients in rows and features in
  columns.

- groups:

  Vector of group labels with length matching `nrow(data_matrix)`.

- p_adjust_method:

  Multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A data frame with test statistics and adjusted p-values.
