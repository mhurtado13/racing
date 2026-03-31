# Resolve cell-type expression profiles

Internal helper that either returns a supplied expression-profile matrix
or computes it from an `estimate_expression_profiles()` function
available in the current R session.

## Usage

``` r
.resolve_expression_profiles(counts, deconv, cell_expr_profile = NULL)
```

## Arguments

- counts:

  Gene-by-sample count matrix.

- deconv:

  Deconvolution result matrix.

- cell_expr_profile:

  Optional precomputed cell-type expression profile.

## Value

A data frame of cell-type expression profiles.
