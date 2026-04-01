# Convert an edge list to a sparse adjacency matrix without self-loops

Convert an edge list to a sparse adjacency matrix without self-loops

## Usage

``` r
EdgetoAdj_No_loop(E, N)
```

## Arguments

- E:

  Two-column edge list with source and target vertex indices.

- N:

  Total number of vertices.

## Value

A sparse adjacency matrix with diagonal entries removed.
