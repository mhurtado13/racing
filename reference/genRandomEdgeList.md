# Sample an edge list from ligand-receptor probabilities

Sample an edge list from ligand-receptor probabilities

## Usage

``` r
genRandomEdgeList(Dligrec, vertextypelist, structurelig, structurerec)
```

## Arguments

- Dligrec:

  Ligand-by-receptor probability matrix.

- vertextypelist:

  Integer vector assigning a cell type to each vertex.

- structurelig:

  Cell-by-ligand compatibility matrix.

- structurerec:

  Cell-by-receptor compatibility matrix.

## Value

A list with `cell_connection` and `ligrec_type` matrices.
