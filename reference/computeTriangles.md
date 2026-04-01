# Compute triangle features from kernel matrices

Compute triangle features from kernel matrices

## Usage

``` r
computeTriangles(
  kernel,
  cell_names,
  patient_names,
  unifKernel = NULL,
  norm = FALSE,
  bundle = TRUE
)
```

## Arguments

- kernel:

  Kernel array from
  [`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md).

- cell_names:

  Character vector of cell-type names.

- patient_names:

  Character vector of patient names.

- unifKernel:

  Optional normalized baseline kernel.

- norm:

  Logical; if `TRUE`, divide by the baseline triangle scores.

- bundle:

  Logical; if `TRUE`, aggregate directionally equivalent triangles.

## Value

A patient-by-feature data frame of triangle scores.
