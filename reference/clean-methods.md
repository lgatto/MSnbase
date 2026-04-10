# Clean 'MSnExp', 'Spectrum' or 'Chromatogram' instances

This method cleans out individual spectra (`Spectrum` instances),
chromatograms
([`Chromatogram`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.md)
instances) or whole experiments (`MSnExp` instances) of 0-intensity
peaks. Unless `all` is set to `FALSE`, original 0-intensity values are
retained only around peaks. If more than two 0's were separating two
peaks, only the first and last ones, those directly adjacent to the peak
ranges are kept. If two peaks are separated by only one 0-intensity
value, it is retained. An illustrative example is shown below.

## Methods

- `signature(object = "MSnExp", all = "logical", verbose = "logical")`:

  Cleans all spectra in `MSnExp` object. Displays a control bar if
  verbose set to `TRUE` (default). Returns a cleaned `MSnExp` instance.

- `signature(object = "Spectrum", all = "logical", msLevel. = "numeric")`:

  Cleans the `Spectrum` object. Returns a cleaned `Spectrum` instance.
  If `all = TRUE`, then all zeros are removed. `msLevel.` defines the
  level of the spectrum, and if `msLevel(object) != msLevel.`, cleaning
  is ignored. Only relevant when called from `OnDiskMSnExp` and is only
  relevant for developers.

- `signature(object = "Chromatogram", all = "logical", na.rm = "logical")`:

  Cleans the
  [`Chromatogram`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.md)
  instance and returns a cleaned
  [`Chromatogram`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.md)
  object. If `na.rm` is `TRUE` (default is `FALSE`) all `NA` intensities
  are removed before cleaning the chromatogram.

## Author

Laurent Gatto

## See also

[`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
and
[`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
for other spectra processing methods.

## Examples

``` r
int <- c(1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0)
sp1 <- new("Spectrum2",
           intensity=int,
           mz=1:length(int))
sp2 <- clean(sp1) ## default is all=FALSE
intensity(sp1)
#>  [1] 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0
intensity(sp2)
#>  [1] 1 0 0 1 1 1 0 0 1 1 0 0 1 0
intensity(clean(sp1, all = TRUE))
#> [1] 1 1 1 1 1 1 1

mz(sp1)
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> [26] 26 27 28 29 30 31 32
mz(sp2)
#>  [1]  1  2  8  9 10 11 12 16 17 18 19 28 29 30
mz(clean(sp1, all = TRUE))
#> [1]  1  9 10 11 17 18 29

data(itraqdata)
itraqdata2 <- clean(itraqdata)
sum(peaksCount(itraqdata))
#> [1] 106289
sum(peaksCount(itraqdata2))
#> [1] 59747
processingData(itraqdata2)
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> Updated from version 0.3.0 to 0.3.1 [Fri Jul  8 20:23:25 2016] 
#> Spectra cleaned: Fri Apr 10 15:48:55 2026 
#>  MSnbase version: 1.1.22 

## Create a simple Chromatogram object
chr <- Chromatogram(rtime = 1:12,
                    intensity = c(0, 0, 20, 0, 0, 0, 123, 124343, 3432, 0, 0, 0))

## Remove 0-intensity values keeping those adjacent to peaks
chr <- clean(chr)
intensity(chr)
#> [1]      0     20      0      0    123 124343   3432      0

## Remove all 0-intensity values
chr <- clean(chr, all = TRUE)
intensity(chr)
#> [1]     20    123 124343   3432

## Clean a Chromatogram with NAs.
chr <- Chromatogram(rtime = 1:12,
                    intensity = c(0, 0, 20, NA, NA, 0, 123, 124343, 3432, 0, 0, 0))
chr <- clean(chr, all = FALSE, na.rm = TRUE)
intensity(chr)
#> [1]      0     20      0    123 124343   3432      0
```
