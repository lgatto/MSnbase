# Calculates coeffivient of variation for features

This function calculates the column-wise coefficient of variation (CV),
i.e. the ration between the standard deviation and the mean, for the
features in an
[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md).
The CVs are calculated for the groups of features defined by `groupBy`.
For groups defined by single features, `NA` is returned.

## Usage

``` r
featureCV(x, groupBy, na.rm = TRUE, norm = "none", suffix = NULL)
```

## Arguments

- x:

  An instance of class
  [`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md).

- groupBy:

  An object of class `factor` defining how to summarise the features.

- na.rm:

  A `logical(1)` defining whether missing values should be removed.

- norm:

  One of normalisation methods applied prior to CV calculation. See
  [`normalise()`](https://lgatto.github.io/MSnbase/reference/normalise-methods.md)
  for more details. Here, the default is `'none'`, i.e. no
  normalisation.

- suffix:

  A `character(1)` to be used to name the new CV columns. Default is
  `NULL` to ignore this. This argument should be set when CV values are
  already present in the
  [`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
  feature variables.

## Value

A `matrix` of dimensions `length(levels(groupBy))` by `ncol(x)` with the
respecive CVs. The column names are formed by pasting `CV.` and the
sample names of object `x`, possibly suffixed by `.suffix`.

## See also

[`combineFeatures()`](https://lgatto.github.io/MSnbase/reference/combineFeatures.md)

## Author

Laurent Gatto and Sebastian Gibb

## Examples

``` r
data(msnset)
msnset <- msnset[1:4]
gb <- factor(rep(1:2, each = 2))
featureCV(msnset, gb)
#>   CV.iTRAQ4.114 CV.iTRAQ4.115 CV.iTRAQ4.116 CV.iTRAQ4.117
#> 1     0.4116294   0.672121019     0.9798589     1.1049021
#> 2     0.1010699   0.005077863     0.1128963     0.2347332
featureCV(msnset, gb, suffix = "2")
#>   CV.iTRAQ4.114.2 CV.iTRAQ4.115.2 CV.iTRAQ4.116.2 CV.iTRAQ4.117.2
#> 1       0.4116294     0.672121019       0.9798589       1.1049021
#> 2       0.1010699     0.005077863       0.1128963       0.2347332
```
