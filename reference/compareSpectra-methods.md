# Compare Spectra of an 'MSnExp' or 'Spectrum' instances

This method compares spectra (`Spectrum` instances) pairwise or all
spectra of an experiment (`MSnExp` instances). Currently the comparison
is based on the number of common peaks `fun = "common"`, the Pearson
correlation `fun = "cor"`, the dot product `fun = "dotproduct"` or a
user-defined function.

For `fun = "common"` the `tolerance` (default `25e-6`) can be set and
the tolerance can be defined to be relative (default `relative = TRUE`)
or absolute (`relative = FALSE`). To compare spectra with `fun = "cor"`
and `fun = "dotproduct"`, the spectra need to be binned. The `binSize`
argument (in Dalton) controls the binning precision. Please see
[`bin`](https://lgatto.github.io/MSnbase/reference/bin-methods.md) for
details.

Instead of these three predefined functions for `fun` a user-defined
comparison function can be supplied. This function takes two
[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
objects as the first two arguments and `...` as third argument. The
function must return a single `numeric` value. See the example section.

## Methods

- `signature(x = "MSnExp", y = "missing", fun = "character", ...)`:

  Compares all spectra in an `MSnExp` object. The `...` arguments are
  passed to the internal functions. Returns a `matrix` of dimension
  `length(x)` by `length(x)`.

- `signature(x = "Spectrum", y = "Spectrum", fun = "character", ...)`:

  Compares two `Spectrum` objects. See the above explanation for `fun`
  and `...`. Returns a single `numeric` value.

## Author

Sebastian Gibb \<mail@sebastiangibb.de\>

## See also

[`bin`](https://lgatto.github.io/MSnbase/reference/bin-methods.md),
[`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md),
[`pickPeaks`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md),
[`smooth`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md),
[`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
and
[`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
for other spectra processing methods.

## References

Stein, S. E., & Scott, D. R. (1994). Optimization and testing of mass
spectral library search algorithms for compound identification. Journal
of the American Society for Mass Spectrometry, 5(9), 859-866. doi:
https://doi.org/10.1016/1044-0305(94)87009-8

Lam, H., Deutsch, E. W., Eddes, J. S., Eng, J. K., King, N., Stein, S.
E. and Aebersold, R. (2007) Development and validation of a spectral
library searching method for peptide identification from MS/MS.
Proteomics, 7: 655-667. doi: https://doi.org/10.1002/pmic.200600625

## Examples

``` r
s1 <- new("Spectrum2", mz=1:10, intensity=1:10)
s2 <- new("Spectrum2", mz=1:10, intensity=10:1)
compareSpectra(s1, s2)
#> [1] 10
compareSpectra(s1, s2, fun="cor", binSize=2)
#> [1] -1
compareSpectra(s1, s2, fun="dotproduct")
#> [1] 0.5714286

## define our own (useless) comparison function (it is just a basic example)
equalLength <- function(x, y, ...) {
  return(peaksCount(x)/(peaksCount(y)+.Machine$double.eps))
}
compareSpectra(s1, s2, fun=equalLength)
#> [1] 1
compareSpectra(s1, new("Spectrum2", mz=1:5, intensity=1:5), fun=equalLength)
#> [1] 2
compareSpectra(s1, new("Spectrum2"), fun=equalLength)
#> [1] 4.5036e+16

data(itraqdata)
compareSpectra(itraqdata[1:5], fun="cor")
#>            X1        X10       X11        X12        X13
#> X1         NA 0.07147120 0.3825327 0.51883926 0.45045283
#> X10 0.0714712         NA 0.0949519 0.09966693 0.07594562
#> X11 0.3825327 0.09495190        NA 0.47739388 0.46612990
#> X12 0.5188393 0.09966693 0.4773939         NA 0.56151180
#> X13 0.4504528 0.07594562 0.4661299 0.56151180         NA
```
