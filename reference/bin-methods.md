# Bin 'MSnExp' or 'Spectrum' instances

This method aggregates individual spectra (`Spectrum` instances) or
whole experiments (`MSnExp` instances) into discrete bins. All intensity
values which belong to the same bin are summed together.

## Methods

- `signature(object = "MSnExp", binSize = "numeric", verbose = "logical")`:

  Bins all spectra in an `MSnExp` object. Use `binSize` to control the
  size of a bin (in Dalton, default is `1`). Displays a control bar if
  verbose set to `TRUE` (default). Returns a binned `MSnExp` instance.

- `signature(object = "Spectrum", binSize = "numeric", breaks = "numeric", msLevel. = "numeric")`:

  Bin the `Spectrum` object. Use `binSize` to control the size of a bin
  (in Dalton, default is `1`). Similar to
  [`hist`](https://rdrr.io/r/graphics/hist.html) you could use `breaks`
  to specify the breakpoints between m/z bins. `msLevel.` defines the
  level of the spectrum, and if `msLevel(object) != msLevel.`, cleaning
  is ignored. Only relevant when called from `OnDiskMSnExp` and is only
  relevant for developers.

  Returns a binned `Spectrum` instance.

## Author

Sebastian Gibb \<mail@sebastiangibb.de\>

## See also

[`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md),
[`pickPeaks`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md),
[`smooth`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md),
[`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
and
[`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
for other spectra processing methods.

## Examples

``` r
s <- new("Spectrum2", mz=1:10, intensity=1:10)
intensity(s)
#>  [1]  1  2  3  4  5  6  7  8  9 10
intensity(bin(s, binSize=2))
#> [1]  3  7 11 15 19

data(itraqdata)
sum(peaksCount(itraqdata))
#> [1] 106289
itraqdata2 <- bin(itraqdata, binSize=2)
sum(peaksCount(itraqdata2))
#> [1] 54186
processingData(itraqdata2)
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> Updated from version 0.3.0 to 0.3.1 [Fri Jul  8 20:23:25 2016] 
#> Spectra binned: Fri Apr 10 14:43:00 2026 
#>  MSnbase version: 1.1.22 
```
