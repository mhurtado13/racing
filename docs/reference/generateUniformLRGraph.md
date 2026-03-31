# Generate a graph under a uniformized ligand-receptor baseline

Generate a graph under a uniformized ligand-receptor baseline

## Usage

``` r
generateUniformLRGraph(
  LRdistr,
  Lmatrix,
  Rmatrix,
  Cdistr,
  cellTypes,
  patient = 1,
  N = 20,
  avdeg = 2
)
```

## Arguments

- LRdistr:

  Ligand-receptor tensor.

- Lmatrix:

  Cell-by-ligand compatibility matrix.

- Rmatrix:

  Cell-by-receptor compatibility matrix.

- Cdistr:

  Patient-by-cell-type abundance matrix.

- cellTypes:

  Character vector of cell-type labels.

- patient:

  Patient index to simulate.

- N:

  Number of cells in the generated graph.

- avdeg:

  Target average degree.

## Value

A list containing the simulated graph and the uniform LR distribution
used.
