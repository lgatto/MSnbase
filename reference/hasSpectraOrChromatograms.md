# Checks if raw data files have any spectra or chromatograms

Helper functions to check whether raw files contain spectra or
chromatograms.

## Usage

``` r
hasSpectra(files)

hasChromatograms(files)
```

## Arguments

- files:

  A [`character()`](https://rdrr.io/r/base/character.html) with raw data
  filenames.

## Value

A `logical(n)` where `n == length(x)` with `TRUE` if that files contains
at least one spectrum, `FALSE` otherwise.

## Author

Laurent Gatto

## Examples

``` r
f <- c(MsDataHub::MRM.standmix.5.mzML(),
       MsDataHub::MS3TMT10_01022016_32917.33481.mzML.gz())
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache
hasSpectra(f)
#> EH10329 EH10410 
#>   FALSE    TRUE 
hasChromatograms(f)
#> EH10329 EH10410 
#>    TRUE    TRUE 
```
