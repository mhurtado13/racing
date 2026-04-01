# Run Monte Carlo simulations for one or more patients

Run Monte Carlo simulations for one or more patients

## Usage

``` r
runSim(
  Lmatrix,
  Rmatrix,
  Cmatrix,
  LRmatrix,
  cells,
  communication_type,
  pats = "all",
  N = 10000,
  itNo = 100,
  av = 20,
  output_folder = NULL,
  file.name = NULL,
  norm = FALSE,
  patient_idx = NULL
)
```

## Arguments

- Lmatrix:

  Cell-by-ligand compatibility matrix.

- Rmatrix:

  Cell-by-receptor compatibility matrix.

- Cmatrix:

  Patient-by-cell-type abundance matrix.

- LRmatrix:

  Ligand-receptor-by-patient tensor.

- cells:

  Character vector of cell-type names.

- communication_type:

  Feature family to simulate (`"D"`, `"W"`, `"TT"`, `"CT"`, or
  `"GSCC"`).

- pats:

  Number of patients to process, or `"all"`.

- N:

  Number of cells per graph.

- itNo:

  Number of Monte Carlo iterations.

- av:

  Target average degree.

- output_folder:

  Directory used to write the `.out` files.

- file.name:

  Output filename stem.

- norm:

  Logical; if `TRUE`, use a uniformized LR baseline.

- patient_idx:

  Optional single patient index to simulate.

## Value

Invisibly writes the simulation outputs to disk.
