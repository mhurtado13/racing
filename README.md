# RaCInG — R package

> **R**andom **C**ell-cell **I**nteraction **G**enerator: an R package for constructing
> patient-specific cell-cell interaction networks from bulk RNA-seq data.

---

## Overview

`RaCInG` (R Analysis of Cell Interactions and Graphs) is an R package that implements
the **RaCInG** random-graph model for reconstructing personalized cell-cell interaction
(CCI) networks from bulk RNA sequencing data. It combines patient-specific cellular
deconvolution estimates and ligand-receptor interaction prior knowledge to generate
probabilistic network ensembles and extract biologically meaningful network features.

This package is a native **R implementation** of the original Python-based pipeline
described in:

> van Santvoort M, Lapuente-Santana Ó, Zopoglou M, Zackl C, Finotello F,
> van der Hoorn P & Eduati F. **Mathematically mapping the network of cells in the tumor
> microenvironment.** *Cell Reports Methods*, 5(2):100985, 2025.  
> DOI: [10.1016/j.crmeth.2025.100985](https://doi.org/10.1016/j.crmeth.2025.100985) | PMID: 39954673

The original Python source code and supplementary materials are available at
[SysBioOncology/RaCInG](https://github.com/SysBioOncology/RaCInG).

---

## Model Description

The RaCInG model operates in four main steps:

1. **Input preparation** — Bulk RNA-seq data are transformed into four matrices:
   - **C-matrix**: cell-type abundance fractions (from deconvolution)
   - **LR-matrix**: ligand-receptor interaction pair activation per patient
   - **L-matrix**: ligand-to-cell-type compatibility
   - **R-matrix**: receptor-to-cell-type compatibility

2. **Network generation** — Patient-specific CCI networks are constructed using one of
   two methods:
   - **Kernel Method** — Efficiently computes the expected interaction strength between
     all pairs of cell types via a closed-form kernel matrix.
   - **Monte Carlo Method** — Generates an ensemble of random network realizations by
     sampling cells and ligand-receptor pairs according to the input distributions.

3. **Feature extraction** — Network fingerprints (direct connections, wedges, trust
   triangles, cycle triangles, and GSCC) are extracted from the ensemble and averaged
   to produce per-patient biomarker values.

4. **Statistical analysis** — Extracted features are tested for associations with
   patient-level metadata (e.g., immunotherapy response) using the Wilcoxon rank-sum
   test with multiple testing correction.

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("mhurtado13/racing")
```

### Dependencies

| Package | Role |
|---|---|
| `data.table` | Fast CSV reading |
| `dplyr` | Data wrangling |
| `tibble` | Tidy data frames |
| `tidyr` | Data reshaping |
| `Matrix` | Sparse adjacency matrices |
| `pracma` | Root-finding utilities (Kernel Method) |
| `ggplot2` | Volcano plots |

---

## Quick Start

```r
library(racing)

# 1. Load input matrices from CSV files
res <- generateInput(
  weight_type = "min",
  cancer_type = "SKCM",
  folder      = "path/to/Example_input"
)

Lmatrix   <- res$CellLigList       # cell-ligand compatibility matrix
Rmatrix   <- res$CellRecList       # cell-receptor compatibility matrix
Cdistr    <- res$Dtypes            # cell-type abundances (patients × cell types)
LRdistr   <- res$DconnectionTensor # LR interaction tensor (ligands × receptors × patients)
cellTypes <- res$celltypes
ligs      <- res$ligands
recs      <- res$receptors

# 2. Generate a single network realization for patient 1
graph <- model1(
  N          = 100,
  avdeg      = 3,
  cellLigList = Lmatrix,
  cellRecList = Rmatrix,
  Dcelltype  = Cdistr[1, ],
  Dligrec    = LRdistr[, , 1],
  genRandom  = FALSE
)
V     <- graph$V     # cell-type label per vertex
E     <- graph$E     # edge list
types <- graph$types # LR type per edge

# 3. Compute the kernel for all patients
kernel <- compute_kernel(Lmatrix, Rmatrix, Cdistr, LRdistr)

# 4. Extract direct communication features
direct_df <- calculateDirect(kernel, cells = cellTypes)

# 5. Statistical analysis
groups <- c(rep("Responder", 10), rep("Non-responder", 10))
results <- wilcox_group_test(direct_df, groups)
volcano_plot(results)
```

---

## Function Reference

### Input generation (`RaCInG_input_generation.R`)

| Function | Description |
|---|---|
| `generateInput()` | Main wrapper: reads all four input matrices from CSV files |
| `createCellLigList()` | Reads the cell–ligand compatibility matrix from a CSV file |
| `createCellRecList()` | Reads the cell–receptor compatibility matrix from a CSV file |
| `createCellTypeDistr()` | Reads cell-type abundance fractions (deconvolution output) |
| `createInteractionDistr()` | Reads LR interaction probabilities and returns a 3-D tensor (ligands × receptors × patients) |
| `Read_Lig_Rec_Interaction()` | Reads ligand-receptor interaction sign/weight table |
| `sortPermute()` | Helper: alphabetically sorts a character vector and returns the permutation index |

---

### Network generation (`network_generation.R`)

| Function | Description |
|---|---|
| `model1()` | Generates a single graph realization of RaCInG given patient-specific distributions |
| `generateUniformLRGraph()` | Generates a graph with a uniform LR distribution (used for normalization) |
| `genRandomCellTypeDistr()` | Generates a random cell-type probability distribution (testing) |
| `genRandomLigRecDistr()` | Generates a random ligand-receptor probability matrix (testing) |
| `genRandomCellLigands()` | Generates a random binary cell-ligand compatibility matrix (testing) |
| `genRandomCellReceptors()` | Generates a random binary cell-receptor compatibility matrix (testing) |
| `genRandomCellTypeList()` | Samples a list of vertex cell-type labels from the cell-type distribution |

---

### Kernel Method (`Kernel_Method.R`)

| Function | Description |
|---|---|
| `compute_kernel()` | Computes the network kernel for each patient: a cell-type × cell-type array of expected interaction strengths |
| `calculateDirect()` | Extracts pairwise direct communication scores from the kernel (one row per patient) |
| `calculateWedges()` | Extracts wedge (two-hop path) scores for all cell-type triplets from the kernel |

---

### Monte Carlo Method (`Monte_Carlo_Method.R`)

| Function | Description |
|---|---|
| `countDirect()` | Counts direct edges by cell-type pair across Monte Carlo iterations |
| `countWedges()` | Counts wedges by cell-type triplet across Monte Carlo iterations |
| `countTrustTriangles()` | Counts trust triangles (v→w, v→u, w→u) by cell-type triplet across MC iterations |
| `countCycleTriangles()` | Counts cycle triangles (v→w→u→v) by cell-type triplet across Monte Carlo iterations |

---

### Feature extraction (`feature_extraction.R`)

| Function | Description |
|---|---|
| `Trust_Triangles()` | Finds all trust triangles in a directed graph and returns count + triangle list |
| `Cycle_Triangles()` | Finds all cycle triangles (directed 3-cycles) and returns count + triangle list |
| `Wedges()` | Finds all wedges (directed 2-paths) and returns count + wedge list |
| `Find_Number_Triangles()` | Counts all triangles (multi-edges counted) via matrix trace |
| `Find_Number_Triangles_Unique()` | Counts unique triangles (multi-edges ignored) |
| `Find_Number_Trust_Triangles_Unique()` | Counts unique trust triangles |
| `Find_Number_2Loops()` | Counts 2-loops (mutual edges, multi-edges counted) |
| `Find_Number_2Loops_Unique()` | Counts unique 2-loops |
| `Find_Number_Wedges()` | Counts all wedges (multi-edges counted) |
| `Find_Number_Wedges_Unique()` | Counts unique wedges |

---

### Utilities (`Utilities.R`)

| Function | Description |
|---|---|
| `EdgetoAdj()` | Converts an edge list to a sparse adjacency matrix (parallel edges summed) |
| `EdgetoAdj_No_loop()` | Converts an edge list to a sparse adjacency matrix, removing self-loops |
| `Count_Types()` | Aggregates a list of graph objects (wedges, triangles) into a count tensor indexed by cell type |

---

### Output processing (`txt_to_csv.R`)

| Function | Description |
|---|---|
| `Read_Sim_Output()` | Parses `.out` simulation output files produced by the pipeline (supports D, W, TT, CT, GSCC feature types) |

---

### Statistical analysis (`statistical_analysis.R`)

| Function | Description |
|---|---|
| `wilcox_group_test()` | Performs Wilcoxon rank-sum tests for each network feature between two patient groups, with FDR correction |
| `volcano_plot()` | Generates a volcano plot (effect size vs. −log₁₀ adjusted p-value) from Wilcoxon test results |

---

## Vignette

A full worked example using skin cutaneous melanoma (SKCM) data is provided in
[vignettes/RaCiNG.Rmd](vignettes/RaCiNG.Rmd). It walks through input loading,
network generation, kernel computation, feature extraction, and statistical analysis.

---

## Original Python Implementation

The original Python code (Monte Carlo, Kernel Method, Circos visualization, HPC CLI)
is available at:

**[SysBioOncology/RaCInG](https://github.com/SysBioOncology/RaCInG)**

---

## Citation

If you use this package, please cite the original publication:

```
van Santvoort M, Lapuente-Santana Ó, Zopoglou M, Zackl C, Finotello F,
van der Hoorn P & Eduati F (2025). Mathematically mapping the network of cells
in the tumor microenvironment. Cell Reports Methods, 5(2), 100985.
https://doi.org/10.1016/j.crmeth.2025.100985
```

BibTeX:

```bibtex
@article{vansantvoort2025racing,
  author  = {van Santvoort, Mike and Lapuente-Santana, \'{O}scar and Zopoglou, Maria
             and Zackl, Constantin and Finotello, Francesca and van der Hoorn, Pim
             and Eduati, Federica},
  title   = {Mathematically mapping the network of cells in the tumor microenvironment},
  journal = {Cell Reports Methods},
  year    = {2025},
  volume  = {5},
  number  = {2},
  pages   = {100985},
  doi     = {10.1016/j.crmeth.2025.100985},
  pmid    = {39954673}
}
```

---

## Authors

**Original model design and Python implementation:**

- Mike van Santvoort
- Óscar Lapuente-Santana
- Maria Zopoglou
- Constantin Zackl
- Francesca Finotello
- Pim van der Hoorn
- Federica Eduati

**R package development:**

- Marcelo Hurtado

---

## License

This package is distributed under the MIT License. See [LICENSE](LICENSE) for details.
