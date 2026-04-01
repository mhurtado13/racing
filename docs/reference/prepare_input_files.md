# Build RaCInG input files from raw count data

This function combines the preprocessing and input-loading steps into a
single call. It generates the `L`, `R`, `C`, and `LR` CSV files from raw
counts, then reads them back to produce the normalised matrices and 3-D
tensor required by the kernel and Monte Carlo workflows.

## Usage

``` r
prepare_input_files(
  counts,
  output_folder = "Results/",
  deconv = NULL,
  cc_network = NULL,
  fun_LR = min,
  cell_expr_profile = NULL,
  source = "source_genesymbol",
  target = "target_genesymbol",
  deconv_method = "Quantiseq",
  cbsx.name = NULL,
  cbsx.token = NULL,
  file_name = NULL,
  signed = FALSE
)
```

## Arguments

- counts:

  Gene-by-sample count matrix.

- output_folder:

  Directory where the generated `L`, `R`, `C`, and `LR` files are
  written.

- deconv:

  Optional deconvolution matrix. If omitted, the function will try to
  compute it.

- cc_network:

  Optional ligand-receptor prior network.

- fun_LR:

  Function used to combine ligand and receptor expression values.

- cell_expr_profile:

  Optional cell-type expression profile matrix.

- source, target:

  Column names to use as ligand and receptor identifiers when
  `cc_network` is supplied.

- deconv_method:

  Deconvolution method passed to
  [`multideconv::compute.deconvolution()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.html).

- cbsx.name, cbsx.token:

  Optional credentials forwarded to the deconvolution workflow.

- file_name:

  File stem used when exporting the generated CSV files.

- signed:

  Logical; if `TRUE`, also try to load a sign matrix from
  `output_folder`.

## Value

A named list with the processed input matrices and their labels:
`Lmatrix`, `Rmatrix`, `Cmatrix` (normalised), `LRmatrix` (3-D tensor),
`celltypes`, `ligands`, `receptors`, `Sign_matrix`, and `CC_table`.
