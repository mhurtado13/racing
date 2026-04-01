# Count graph motifs by cell-type combination

Count graph motifs by cell-type combination

## Usage

``` r
Count_Types(oblist, V, maxTypes = 0)
```

## Arguments

- oblist:

  Matrix of graph objects where each row contains vertex indices.

- V:

  Integer vector mapping vertices to cell types.

- maxTypes:

  Optional maximum number of types used to size the output array.

## Value

An array counting how often each type combination occurs.
