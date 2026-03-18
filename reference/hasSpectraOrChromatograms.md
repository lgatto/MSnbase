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
f <- msdata::proteomics(full.names = TRUE)[1:2]
hasSpectra(f)
#>                /__w/_temp/Library/msdata/proteomics/MRM-standmix-5.mzML.gz 
#>                                                                      FALSE 
#> /__w/_temp/Library/msdata/proteomics/MS3TMT10_01022016_32917-33481.mzML.gz 
#>                                                                       TRUE 
hasChromatograms(f)
#>                /__w/_temp/Library/msdata/proteomics/MRM-standmix-5.mzML.gz 
#>                                                                       TRUE 
#> /__w/_temp/Library/msdata/proteomics/MS3TMT10_01022016_32917-33481.mzML.gz 
#>                                                                       TRUE 
```
