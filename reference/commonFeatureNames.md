# Keep only common feature names

Subsets `MSnSet` instances to their common feature names.

## Usage

``` r
commonFeatureNames(x, y)
```

## Arguments

- x:

  An instance of class
  [`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
  or a `list` or `MSnSetList` with at least 2 `MSnSet` objects.

- y:

  An instance of class
  [`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md).
  Ignored if `x` is a `list`/`MSnSetList`.

## Value

An `linkS4class{MSnSetList}` composed of the input `MSnSet` containing
only common features in the same order. The names of the output are
either the names of the `x` and `y` input variables or the names of `x`
if a list is provided.

## Author

Laurent Gatto

## Examples

``` r
library("pRolocdata")
data(tan2009r1)
data(tan2009r2)
cmn <- commonFeatureNames(tan2009r1, tan2009r2)
#> 551 features in common
names(cmn)
#> [1] "tan2009r1" "tan2009r2"
## as a named list
names(commonFeatureNames(list(a = tan2009r1, b = tan2009r2)))
#> 551 features in common
#> [1] "a" "b"
## without message
suppressMessages(cmn <- commonFeatureNames(tan2009r1, tan2009r2))
## more than 2 instance
data(tan2009r3)
cmn <- commonFeatureNames(list(tan2009r1, tan2009r2, tan2009r3))
#> 404 features in common
length(cmn)
#> [1] 3
```
