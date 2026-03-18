# Noise Estimation for 'Spectrum' instances

This method performs a noise estimation on individual spectra
(`Spectrum` instances). There are currently two different noise
estimators, the Median Absolute Deviation (`method = "MAD"`) and
Friedman's Super Smoother (`method = "SuperSmoother"`), as implemented
in the
[`MALDIquant::detectPeaks`](https://rdrr.io/pkg/MALDIquant/man/detectPeaks-methods.html)
and
[`MALDIquant::estimateNoise`](https://rdrr.io/pkg/MALDIquant/man/estimateNoise-methods.html)
functions respectively.

## Methods

- `signature(object = "Spectrum", method = "character", ...)`:

  Estiamtes the noise in a non-centroided spectrum (`Spectrum`
  instance). `method` could be `"MAD"` or `"SuperSmoother"`. The
  arguments `...` are passed to the noise estimator functions
  implemented in
  [`MALDIquant::estimateNoise`](https://rdrr.io/pkg/MALDIquant/man/estimateNoise-methods.html).
  Currenlty only the `method = "SuperSmoother"` accepts additional
  arguments, e.g. `span`. Please see
  [`supsmu`](https://rdrr.io/r/stats/supsmu.html) for details. This
  method returns a two-column matrix with the m/z and intensity values
  in the first and the second column.

- `signature(object = "MSnExp", method = "character", ...)`:

  Estimates noise for all spectra in `object`.

## Author

Sebastian Gibb \<mail@sebastiangibb.de\>

## See also

[`pickPeaks`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md),
and the underlying method in `MALDIquant`: `estimateNoise`.

## References

S. Gibb and K. Strimmer. 2012. MALDIquant: a versatile R package for the
analysis of mass spectrometry data. Bioinformatics 28: 2270-2271.
<http://strimmerlab.org/software/maldiquant/>

## Examples

``` r
sp1 <- new("Spectrum1",
           intensity = c(1:6, 5:1),
           mz = 1:11,
           centroided = FALSE)
estimateNoise(sp1, method = "SuperSmoother")
#>       mz intensity
#>  [1,]  1      1.08
#>  [2,]  2      2.00
#>  [3,]  3      2.92
#>  [4,]  4      3.68
#>  [5,]  5      4.20
#>  [6,]  6      4.40
#>  [7,]  7      4.20
#>  [8,]  8      3.68
#>  [9,]  9      2.92
#> [10,] 10      2.00
#> [11,] 11      1.08
```
