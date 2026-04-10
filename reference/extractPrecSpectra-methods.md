# Extracts precursor-specific spectra from an 'MSnExp' object

Extracts the MSMS spectra that originate from the precursor(s) having
the same MZ value as defined in the`prec` argument.

A warning will be issued of one or several of the precursor MZ values in
`prec` are absent in the experiment precursor MZ values (i.e in
`precursorMz(object)`).

## Methods

- `signature(object = "MSnExp", prec = "numeric")`:

  Returns an
  `"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
  containing MSMS spectra whose precursor MZ values are in `prec`.

## Author

Laurent Gatto

## Examples

``` r
file <- dir(system.file(package="MSnbase",dir="extdata"),
            full.name=TRUE,pattern="mzXML$")
aa <- readMSData(file,verbose=FALSE)
my.prec <- precursorMz(aa)[1]
my.prec
#>    F1.S1 
#> 645.3741 
bb <- extractPrecSpectra(aa,my.prec)
precursorMz(bb)
#>    F1.S1    F1.S3 
#> 645.3741 645.3741 
processingData(bb)
#> - - - Processing information - - -
#> Data loaded: Fri Apr 10 14:44:22 2026 
#> 1 (2) precursors (spectra) extracted: Fri Apr 10 14:44:22 2026 
#>  MSnbase version: 2.37.2 
```
