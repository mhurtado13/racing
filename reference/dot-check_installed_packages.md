# Check that optional packages are installed

Internal helper used by high-level workflows that depend on optional
preprocessing packages.

## Usage

``` r
.check_installed_packages(packages)
```

## Arguments

- packages:

  Character vector of package names.

## Value

Invisibly returns `TRUE` when all packages are available.
