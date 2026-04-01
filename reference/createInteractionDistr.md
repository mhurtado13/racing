# Read ligand-receptor interaction probabilities

Read ligand-receptor interaction probabilities

## Usage

``` r
createInteractionDistr(filename, ligands, receptors)
```

## Arguments

- filename:

  Path to a CSV file containing patient-by-interaction weights.

- ligands:

  Character vector of ligand names.

- receptors:

  Character vector of receptor names.

## Value

A 3D array with dimensions ligand × receptor × patient.
