# Smooths 'MSnExp' or 'Spectrum' instances

This method smooths individual spectra (`Spectrum` instances) or whole
experiments (`MSnExp` instances). Currently, the
Savitzky-Golay-Smoothing (`method = "SavitzkyGolay"`) and the
Moving-Average-Smoothing (`method = "MovingAverage"`) are available, as
implemented in the
[`MALDIquant::smoothIntensity`](https://rdrr.io/pkg/MALDIquant/man/smoothIntensity-methods.html)
function. Additional methods might be added at a later stage.

## Methods

- `signature(x = "MSnExp", method = "character", halfWindowSize = "integer", verbose = "logical", ...)`:

  Smooths all spectra in `MSnExp`. `method` could be `"SavitzkyGolay"`
  or `"MovingAverage"`. `"halfWindowSize"` controls the window size of
  the filter. The resulting window size is `2 * halfWindowSize + 1`. The
  best size differs depending on the selected `method`. For
  `method = "SavitzkyGolay"` it should be lower than *FWHM* of the peaks
  (full width at half maximum; please find details in Bromba and Ziegler
  1981). The arguments `...` are passed to the internal functions. For
  `method="MovingAverage"` there is an additional `weighted` argument
  (default: `FALSE`) to indicate if the average should be equal weight
  (default) or if it should have weights depending on the distance from
  the center as calculated as `1/2^abs(-halfWindowSize:halfWindowSize)`
  with the sum of all weigths normalized to 1. For
  `method="SavitzkyGolay"` an additonal argument is `polynomialOrder`
  (default: 3). It controls the polynomial order of the Savitzky-Golay
  Filter. This method displays a progress bar if `verbose = TRUE`.
  Returns an `MSnExp` instance with smoothed spectra.

- `signature(x = "Spectrum", method = "character", halfWindowSize = "integer", ...)`:

  Smooths the spectrum (`Spectrum` instance). This method is the same as
  above but returns a smoothed `Spectrum` instead of an `MSnExp` object.
  It has no `verbose` argument. Please read the details for the above
  `MSnExp` method.

## Author

Sebastian Gibb \<mail@sebastiangibb.de\>

## See also

[`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md),
[`pickPeaks`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md),
[`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
and
[`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
for other spectra processing methods.

## References

A. Savitzky and M. J. Golay. 1964. Smoothing and differentiation of data
by simplified least squares procedures. Analytical chemistry, 36(8),
1627-1639.

M. U. Bromba and H. Ziegler. 1981. Application hints for Savitzky-Golay
digital smoothing filters. Analytical Chemistry, 53(11), 1583-1586.

S. Gibb and K. Strimmer. 2012. MALDIquant: a versatile R package for the
analysis of mass spectrometry data. Bioinformatics 28: 2270-2271.
<http://strimmerlab.org/software/maldiquant/>

## Examples

``` r
sp1 <- new("Spectrum1",
           intensity = c(1:6, 5:1),
           mz = 1:11)
sp2 <- smooth(sp1, method = "MovingAverage", halfWindowSize = 2)
intensity(sp2)
#>  [1] 3.0 3.0 3.0 4.0 4.6 4.8 4.6 4.0 3.0 3.0 3.0

data(itraqdata)
itraqdata2 <- smooth(itraqdata, 
                     method = "MovingAverage", 
                     halfWindowSize = 2)
processingData(itraqdata2)
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> Updated from version 0.3.0 to 0.3.1 [Fri Jul  8 20:23:25 2016] 
#> Spectra smoothed (MovingAverage): Fri Apr 10 15:51:01 2026 
#>  MSnbase version: 1.1.22 
```
