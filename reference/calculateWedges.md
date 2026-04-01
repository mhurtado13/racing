# Calculate wedge features from a kernel

Calculate wedge features from a kernel

## Usage

``` r
calculateWedges(kernel, unifKernel = NULL, cells, bundle = TRUE)
```

## Arguments

- kernel:

  Kernel array returned by
  [`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md).

- unifKernel:

  Optional normalized baseline kernel.

- cells:

  Character vector of cell-type names.

- bundle:

  Logical; if `TRUE`, combine directionally equivalent wedges.

## Value

A patient-by-feature data frame of wedge scores.
