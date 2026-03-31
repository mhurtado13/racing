# Generate a single RaCInG graph realization

Generate a single RaCInG graph realization

## Usage

``` r
model1(
  N,
  avdeg,
  cellLigList = NULL,
  cellRecList = NULL,
  Dcelltype = NULL,
  Dligrec = NULL,
  Signmatrix = NULL,
  genRandom = TRUE
)
```

## Arguments

- N:

  Number of vertices (cells) in the graph.

- avdeg:

  Target average degree.

- cellLigList:

  Cell-by-ligand compatibility matrix.

- cellRecList:

  Cell-by-receptor compatibility matrix.

- Dcelltype:

  Cell-type abundance probabilities.

- Dligrec:

  Ligand-by-receptor probability matrix.

- Signmatrix:

  Optional ligand-receptor sign matrix.

- genRandom:

  Logical; if `TRUE`, generate random test inputs internally.

## Value

A list with vertex labels, an edge list, and ligand-receptor types.
