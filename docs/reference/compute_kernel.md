# Compute the RaCInG kernel for one or more patients

Compute the RaCInG kernel for one or more patients

## Usage

``` r
compute_kernel(liglist, reclist, Cmatrix, LRmatrix, normalize = FALSE)
```

## Arguments

- liglist:

  Cell-by-ligand compatibility matrix.

- reclist:

  Cell-by-receptor compatibility matrix.

- Cmatrix:

  Patient-by-cell-type abundance matrix.

- LRmatrix:

  Ligand-by-receptor-by-patient interaction tensor.

- normalize:

  Logical; if `TRUE`, also compute a uniformized baseline kernel.

## Value

Either a 3D kernel array or a list with `kernel` and `kernel_norm`.
