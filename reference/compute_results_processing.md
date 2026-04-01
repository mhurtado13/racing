# Convert raw simulation outputs into feature matrices

Convert raw simulation outputs into feature matrices

## Usage

``` r
compute_results_processing(
  celltypes,
  patient_names,
  triangle_type,
  remove_direction = TRUE,
  normalized = TRUE,
  sim_raw_file = NULL,
  sim_norm_file = NULL,
  output_folder = NULL,
  file.name = NULL
)
```

## Arguments

- celltypes:

  Character vector of cell-type names.

- patient_names:

  Character vector of patient identifiers.

- triangle_type:

  Communication type (`"D"`, `"W"`, `"TT"`, `"CT"`, or `"GSCC"`).

- remove_direction:

  Logical; if `TRUE`, merge directionally equivalent features.

- normalized:

  Logical; if `TRUE`, divide raw values by a uniform baseline.

- sim_raw_file:

  Path to the raw simulation output file.

- sim_norm_file:

  Optional path to the normalized simulation output file.

- output_folder:

  Directory where the processed CSV files are written.

- file.name:

  Output filename stem.

## Value

A list with processed mean and standard-deviation data frames.
