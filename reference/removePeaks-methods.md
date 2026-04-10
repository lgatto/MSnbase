# Removes low intensity peaks

This method sets low intensity peaks from individual spectra (`Spectrum`
instances) or whole experiments (`MSnExp` instances) to 0. The intensity
threshold is set with the `t` parameter. Default is the `"min"`
character. The threshold is then set as the non-0 minimum intensity
found in the spectrum. Any other numeric values is valid. All peaks with
maximum intensity smaller or equal to `t` are set to 0.

If the spectrum is in profile mode, ranges of successive non-0 peaks \<=
`t` are set to 0. If the spectrum is centroided, then individual peaks
\<= `t` are set to 0. See the example below for an illustration.

Note that the number of peaks is not changed; the peaks below the
threshold are set to 0 and the object is not cleanded out (see
[`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)).
An illustrative example is shown below.

## Methods

- `signature(object = "MSnExp", t, verbose = "logical" )`:

  Removes low intensity peaks of all spectra in `MSnExp` object. `t`
  sets the minimum peak intensity. Default is "min", i.e the smallest
  intensity in each spectrum. Other `numeric` values are valid. Displays
  a control bar if verbose set to `TRUE` (default). Returns a new
  `MSnExp` instance.

- `signature(object = "Spectrum", t, msLevel. = "numeric")`:

  Removes low intensity peaks of `Spectrum` object. `t` sets the minimum
  peak intensity. Default is "min", i.e the smallest intensity in each
  spectrum. Other `numeric` values are valid. `msLevel.` defines the
  level of the spectrum, and if `msLevel(object) != msLevel.`, cleaning
  is ignored. Only relevant when called from `OnDiskMSnExp` and is only
  relevant for developers.

  Returns a new `Spectrum` instance.

## Author

Laurent Gatto

## See also

[`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)
and
[`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
for other spectra processing methods.

## Examples

``` r
int <- c(2, 0, 0, 0, 1, 5, 1, 0, 0, 1, 3, 1, 0, 0, 1, 4, 2, 1)
sp1 <- new("Spectrum2",
           intensity = int,
           mz = 1:length(int),
           centroided = FALSE)
sp2 <- removePeaks(sp1) ## no peaks are removed here
                        ## as min intensity is 1 and
                        ## no peak has a max int <= 1
sp3 <- removePeaks(sp1, 3)
intensity(sp1)
#>  [1] 2 0 0 0 1 5 1 0 0 1 3 1 0 0 1 4 2 1
intensity(sp2)
#>  [1] 2 0 0 0 1 5 1 0 0 1 3 1 0 0 1 4 2 1
intensity(sp3)
#>  [1] 0 0 0 0 1 5 1 0 0 0 0 0 0 0 1 4 2 1

peaksCount(sp1) == peaksCount(sp2)
#> [1] TRUE
peaksCount(sp3) <= peaksCount(sp1)
#> [1] TRUE

data(itraqdata)
itraqdata2 <- removePeaks(itraqdata, t = 2.5e5)
table(unlist(intensity(itraqdata)) == 0)
#> 
#> FALSE  TRUE 
#> 45628 60661 
table(unlist(intensity(itraqdata2)) == 0)
#> 
#> FALSE  TRUE 
#> 10280 96009 
processingData(itraqdata2)
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> Updated from version 0.3.0 to 0.3.1 [Fri Jul  8 20:23:25 2016] 
#> Curves <= 250000 set to '0': Fri Apr 10 14:45:13 2026 
#>  MSnbase version: 1.1.22 

## difference between centroided and profile peaks

int <- c(104, 57, 32, 33, 118, 76, 38, 39, 52, 140, 52, 88, 394, 71,
         408, 94, 2032)
sp <- new("Spectrum2",
          intensity = int,
          centroided = FALSE,
          mz = seq_len(length(int)))

## unchanged, as ranges of peaks <= 500 considered
intensity(removePeaks(sp, 500))
#>  [1]  104   57   32   33  118   76   38   39   52  140   52   88  394   71  408
#> [16]   94 2032
stopifnot(identical(intensity(sp), intensity(removePeaks(sp, 500))))

centroided(sp) <- TRUE
## different!
intensity(removePeaks(sp, 500))
#>  [1]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> [16]    0 2032
```
