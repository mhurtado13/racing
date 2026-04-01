# Getting started with RaCInG

## Overview

`RaCInG` reconstructs patient-specific cell-cell communication networks
from bulk RNA-seq data. The package supports two complementary
workflows:

- a **kernel-based** approach for fast deterministic feature extraction,
  and
- a **Monte Carlo** approach for simulation-based network summaries.

This vignette shows the recommended workflow and the most important
entry points for new users.

## Installation

### Install from GitHub

``` r
# install.packages("remotes")
remotes::install_github("mhurtado13/racing")
library(RaCInG)
```

### Install from a local checkout

``` r
# install.packages("devtools")
devtools::install(".")
library(RaCInG)
```

If you want to build the RaCInG input matrices directly from raw counts,
install the optional preprocessing dependencies used by
[`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md):

``` r
install.packages(c("ggplot2", "OmnipathR"))
# Additional optional packages: ADImpute, multideconv, liana
```

## Workflow at a glance

| Goal                                 | Function                                                                                                    | Output                                                  |
|--------------------------------------|-------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| Build input matrices from raw counts | [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)             | Named list with `L`, `R`, `C`, `LR` matrices and labels |
| Compute deterministic features       | [`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md)         | Kernel arrays + feature matrix                          |
| Compute simulation summaries         | [`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md) | Processed Monte Carlo results                           |
| Compare clinical groups              | [`wilcox_group_test()`](https://mhurtado13.github.io/racing/reference/wilcox_group_test.md)                 | Statistics table and volcano plot                       |

## Recommended workflow

### 1. Build input matrices from raw counts

Use
[`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)
to generate, save, and load the input matrices in a single call.

``` r
input <- prepare_input_files(
  counts = counts_matrix,
  output_folder = "Results/",
  file_name = "SKCM"
)

str(input)
```

### 2. Run the kernel method

The kernel method is the fastest way to derive direct, wedge, triangle,
or GSCC features across patients. You can pass `counts` to let the
function compute inputs automatically, or supply previously computed
matrices via `input_data`.

``` r
# Option A: from raw counts (runs prepare_input_files internally)
kernel_res <- compute_racing_kernel(
  counts = counts_matrix,
  output_folder = tempdir(),
  file_name = "SKCM",
  communication_type = "W",
  norm = TRUE
)

# Option B: from pre-computed input matrices (skips preprocessing)
kernel_res <- compute_racing_kernel(
  input_data = input,
  communication_type = "W",
  norm = TRUE
)

head(kernel_res$features[, 1:5])
```

### 3. Run the Monte Carlo method

Use the Monte Carlo workflow when you want simulation-based summaries or
uncertainty estimates from repeated graph realizations. The same
`input_data` shortcut is available here.

``` r
set.seed(1)
mc_res <- compute_racing_montecarlo(
  counts = counts_matrix,
  output_folder = tempdir(),
  deconv_method = "Quantiseq",
  file_name = "SKCM",
  nPatients = 3,
  communication_type = "W",
  Ncells = 100,
  Ngraphs = 10,
  Ndegree = 3,
  norm = TRUE
)

# Or from pre-computed inputs:
mc_res <- compute_racing_montecarlo(
  input_data = input,
  output_folder = tempdir(),
  file_name = "SKCM",
  communication_type = "W",
  Ncells = 100,
  Ngraphs = 10,
  Ndegree = 3,
  norm = TRUE
)

head(mc_res$output$mean[, 1:5])
```

### 4. Perform statistical testing

Once features are available in a patient-by-feature matrix, use the
built-in Wilcoxon workflow to compare groups.

``` r
grouping <- c("Responder", "Responder", "Non-responder", "Non-responder")
wilcox_results <- wilcox_group_test(kernel_res$features, grouping)
head(wilcox_results)
volcano_plot(wilcox_results, top_labels = 15)
```

## Notes

- [`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md)
  and
  [`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md)
  are the main entry points.
- Both accept an `input_data` argument with pre-computed matrices (as
  returned by
  [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)),
  which skips all preprocessing.
- [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)
  requires extra packages for deconvolution and prior-network assembly.
- The original Python implementation is available at
  <https://github.com/SysBioOncology/RaCInG>.

## Running with the bundled example data

The package ships with `skcm_example`, a pre-processed named list
derived from 10 TCGA SKCM melanoma patients. It can be passed directly
to the kernel or Monte Carlo workflows via the `input_data` parameter.

``` r
library(RaCInG)
data(skcm_example)
str(skcm_example, max.level = 1)
#> List of 8
#>  $ Lmatrix    : num [1:9, 1:276] 1 0 1 1 1 1 1 1 0 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ Rmatrix    : num [1:9, 1:298] 1 0 1 1 0 1 1 1 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ Cmatrix    : num [1:10, 1:9] 0.00374 0.01407 0.02119 0 0.00313 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ LRmatrix   : num [1:276, 1:298, 1:10] 0.000228 0 0 0 0 ...
#>  $ celltypes  : chr [1:9] "B" "CAF" "CD8+ T" "DC" ...
#>  $ ligands    : chr [1:276] "LGALS9" "ADAM10" "TNFSF12" "ICOSLG" ...
#>  $ receptors  : chr [1:298] "PTPRC" "MET" "CD44" "LRP1" ...
#>  $ Sign_matrix: num [1:276, 1:298] 0 0 0 0 0 0 0 0 0 0 ...
```

### Kernel method on the example data

``` r
kernel_res <- compute_racing_kernel(
  input_data   = skcm_example,
  output_folder = tempdir(),
  communication_type = "W",
  norm = TRUE
)

head(kernel_res$features[, 1:5])
```

### Monte Carlo method on the example data

``` r
set.seed(42)
mc_res <- compute_racing_montecarlo(
  input_data   = skcm_example,
  output_folder = tempdir(),
  file_name     = "skcm_example",
  nPatients     = 2,
  communication_type = "W",
  Ncells  = 100,
  Ngraphs = 5,
  Ndegree = 3,
  norm    = TRUE
)
```
