# Read and normalize cell-type abundance estimates

Read and normalize cell-type abundance estimates

## Usage

``` r
createCellTypeDistr(cells, filename)
```

## Arguments

- cells:

  Character vector of expected cell types.

- filename:

  Path to a CSV file containing patient-by-cell-type abundances.

## Value

A list with the normalized `Dtypes` matrix plus labels.
