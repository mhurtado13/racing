# Load RaCInG input matrices from disk

Load RaCInG input matrices from disk

## Usage

``` r
generateInput(file_name, output_folder, read_signs = FALSE)
```

## Arguments

- file_name:

  File stem used for the `Lmatrix_`, `Rmatrix_`, `Cmatrix_`, and
  `LRmatrix_` CSV files.

- output_folder:

  Directory containing the exported input files.

- read_signs:

  Logical; if `TRUE`, attempts to read `Sign_matrix_<file_name>.csv`
  from `output_folder`.

## Value

A named list containing the input matrices and their labels.
