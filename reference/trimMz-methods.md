# Trims 'MSnExp' or 'Spectrum' instances

This method selects a range of MZ values in a single spectrum
(`Spectrum` instances) or all the spectra of an experiment (`MSnExp`
instances). The regions to trim are defined by the range of `mz`
argument, such that MZ values \<= `min(mz)` and MZ values \>= `max(mz)`
are trimmed away.

## Methods

- `signature(object = "MSnExp", mz = "numeric", msLevel. = "numeric")`:

  Trims all spectra in `MSnExp` object according to `mz`. If `msLevel.`
  is defined, then only spectra of that level are trimmer.

- `signature(object = "Spectrum", mz = "numeric", msLevel. = "numeric")`:

  Trims the `Spectrum` object and retruns a new trimmed object.
  `msLevel.` defines the level of the spectrum, and if
  `msLevel(object) != msLevel.`, cleaning is ignored. Only relevant when
  called from `OnDiskMSnExp` and is only relevant for developers.

## Author

Laurent Gatto

## See also

[`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
and
[`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)
for other spectra processing methods.

## Examples

``` r
mz <- 1:100
sp1 <- new("Spectrum2",
           mz = mz,
           intensity = abs(rnorm(length(mz))))

sp2 <- trimMz(sp1, c(25, 75))
range(mz(sp1))
#> [1]   1 100
range(mz(sp2))
#> [1] 25 75

data(itraqdata)
itraqdata2 <- filterMz(itraqdata, c(113, 117))
range(mz(itraqdata))
#> [1]   99.99872 2069.27344
range(mz(itraqdata2))
#> [1] 113.0459 116.7362
processingData(itraqdata2)
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> Updated from version 0.3.0 to 0.3.1 [Fri Jul  8 20:23:25 2016] 
#> Filter: trim MZ [113..117] on MS level(s) 2. 
#>  MSnbase version: 1.1.22 
```
