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

| Goal                              | Function                                                                                                    | Output                            |
|-----------------------------------|-------------------------------------------------------------------------------------------------------------|-----------------------------------|
| Build input files from raw counts | [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)             | `L`, `R`, `C`, and `LR` matrices  |
| Load prepared inputs              | [`generateInput()`](https://mhurtado13.github.io/racing/reference/generateInput.md)                         | Named list of matrices and labels |
| Compute deterministic features    | [`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md)         | Kernel arrays + feature matrix    |
| Compute simulation summaries      | [`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md) | Processed Monte Carlo results     |
| Compare clinical groups           | [`wilcox_group_test()`](https://mhurtado13.github.io/racing/reference/wilcox_group_test.md)                 | Statistics table and volcano plot |

## Recommended workflow

### 1. Start from prepared input matrices

If you already have precomputed `L`, `R`, `C`, and `LR` matrices, load
them with
[`generateInput()`](https://mhurtado13.github.io/racing/reference/generateInput.md).

``` r
input <- generateInput(
  file_name = "SKCM",
  output_folder = "Results/"
)

str(input)
```

### 2. Run the kernel method

The kernel method is the fastest way to derive direct, wedge, triangle,
or GSCC features across patients.

``` r
kernel_res <- compute_racing_kernel(
  counts = counts_matrix,
  output_folder = tempdir(),
  file_name = "SKCM",
  communication_type = "W",
  norm = TRUE
)

head(kernel_res$features[, 1:5])
```

### 3. Run the Monte Carlo method

Use the Monte Carlo workflow when you want simulation-based summaries or
uncertainty estimates from repeated graph realizations.

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
- [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)
  requires extra packages for deconvolution and prior-network assembly.
- The original Python implementation is available at
  <https://github.com/SysBioOncology/RaCInG>.
