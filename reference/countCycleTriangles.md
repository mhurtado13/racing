# Count cycle triangles across Monte Carlo graph simulations

Count cycle triangles across Monte Carlo graph simulations

## Usage

``` r
countCycleTriangles(Dcell, Dconn, lig, rec, cellnames, N, av, itNo)
```

## Arguments

- Dcell:

  Cell-type abundance vector for one patient.

- Dconn:

  Ligand-receptor probability matrix for one patient.

- lig:

  Cell-by-ligand compatibility matrix.

- rec:

  Cell-by-receptor compatibility matrix.

- cellnames:

  Character vector of cell-type names.

- N:

  Number of cells per simulated graph.

- av:

  Target average degree.

- itNo:

  Number of Monte Carlo iterations.

## Value

A list of average and standard-deviation cycle-triangle counts.
