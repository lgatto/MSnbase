# Count the number of quantitfied features.

This function counts the number of quantified features, i.e non NA
quantitation values, for each group of features for all the samples in
an
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
object. The group of features are defined by a feature variable names,
i.e the name of a column of `fData(object)`.

## Usage

``` r
nQuants(x, groupBy)
```

## Arguments

- x:

  An instance of class
  `"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`.

- groupBy:

  An object of class `factor` defining how to summerise the features.
  (Note that this parameter was previously named `fcol` and referred to
  a feature variable label. This has been updated in version 1.19.12 for
  consistency with other functions.)

## Value

A `matrix` of dimensions `length(levels(groupBy))` by `ncol(x)`

A `matrix` of dimensions `length(levels(factor(fData(object)[, fcol])))`
by `ncol(object)` of integers.

## Details

This function is typically used after
[`topN`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md) and
before
[`combineFeatures`](https://lgatto.github.io/MSnbase/reference/combineFeatures.md),
when the summerising function is `sum`, or any function that does not
normalise to the number of features aggregated. In the former case, sums
of features might be the result of 0 (if no feature was quantified) to
`n` (if all `topN`'s `n` features were quantified) features, and one
might want to rescale the sums based on the number of non-NA features
effectively summed.

## Author

Laurent Gatto <lg390@cam.ac.uk>, Sebastian Gibb <mail@sebastiangibb.de>

## Examples

``` r
data(msnset)
n <- 2
msnset <- topN(msnset, groupBy = fData(msnset)$ProteinAccession, n)
m <- nQuants(msnset, groupBy = fData(msnset)$ProteinAccession)
msnset2 <- combineFeatures(msnset,
                           groupBy = fData(msnset)$ProteinAccession,
                           method = sum)
stopifnot(dim(n) == dim(msnset2))
head(exprs(msnset2))
#>         iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> BSA       2941.052   5174.721   8955.985  18023.287
#> ECA0172  17593.548  18545.620  19361.837  18328.237
#> ECA0435   9847.256  11115.636  11550.405  10158.590
#> ECA0452   1524.148   1399.897   1547.218   1563.230
#> ECA0469   2139.889   2071.378   2058.840   1999.391
#> ECA0621   1101.062   1124.167   1140.093   1191.805
head(exprs(msnset2) * (n/m))
#>         iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> BSA       2941.052   5174.721   8955.985  18023.287
#> ECA0172  35187.095  37091.239  38723.673  36656.473
#> ECA0435   9847.256  11115.636  11550.405  10158.590
#> ECA0452   3048.295   2799.794   3094.437   3126.460
#> ECA0469   2139.889   2071.378   2058.840   1999.391
#> ECA0621   2202.124   2248.334   2280.186   2383.611
```
