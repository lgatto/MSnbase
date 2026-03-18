# Class `"FeatComp"`

Comparing feature names of two comparable `MSnSet` instances.

## Objects from the Class

Objects can be created with `compfnames`. The method compares the
feature names of two objects of class `"MSnSet"`. It prints a summary
matrix of common and unique feature names and invisibly returns a list
of `FeatComp` instances.

The function will compute the common and unique features for all feature
names of the two input objects (`featureNames(x)` and `feautreNames(y)`)
as well as distinct subsets as defined in the `fcol1` and `fcol2`
feautre variables.

## Slots

- `name`::

  Object of class `"character"` defining the name of the compared
  features. By convention, `"all"` is used when all feature names are
  used; otherwise, the respective levels of the feature variables
  `fcol1` and `fcol2`.

- `common`::

  Object of class `"character"` with the common feature names.

- `unique1`::

  Object of class `"character"` with the features unique to the first
  `MSnSet` (`x` in `compfname`).

- `unique2`::

  Object of class `"character"` with the features unique to the seconn
  `MSnSet` (`y` in `compfname`).

- `all`::

  Object of class `"logical"` defining if all features of only a subset
  were compared. One expects that `name == "all"` when `all` is `TRUE`.

## Methods

Accessors `names`, `common`, `unique1` and `unique2` can be used to
access the respective `FeatComp` slots.

- compfnames:

  `signature(x = "MSnSet", y = "MSnSet", fcol1 = "character", fcol2 = "character", simplify = "logical", verbose = "logical")`:
  creates the `FeatComp` comparison object for instances `x` and `y`.
  The feature variables to be considered to details feature comparison
  can be defined by `fcol1` (default is `"markers"` and `fcol2` for `x`
  and `y` respectively). Setting either to `NULL` will only consider all
  feature names; in such case, of `simplify` is `TRUE` (default), an
  `FeatComp` object is returned instead of a list of length 1. The
  `verbose` logical controls if a summary table needs to be printed
  (default is `TRUE`).

- compfnames:

  `signature(x = "list", y = "missing", ...)`: when `x` is a list of
  `MSnSet` instances, `compfnames` is applied to all element pairs of
  `x`. Additional parameters `fcol1`, `fcol2`, `simplify` and `verbose`
  are passed to the pairwise comparison method.

- show:

  `signature(object = "FeatComp")`: prints a summary of the `object`.

## Author

Laurent Gatto and Thomas Naake

## See also

[`averageMSnSet`](https://lgatto.github.io/MSnbase/reference/averageMSnSet.md)
to compuate an average `MSnSet`.

## Examples

``` r
library("pRolocdata")
#> 
#> This is pRolocdata version 1.49.0.
#> Use 'pRolocdata()' to list available data sets.
data(tan2009r1)
data(tan2009r2)
x <- compfnames(tan2009r1, tan2009r2)
x[[1]]
#> Object of class "FeatComp", 'all' features:
#>  Common feature: 551 
#>  Unique to 1: 337 
#>  Unique to 2: 320 
x[2:3]
#> [[1]]
#> Object of class "FeatComp", 'unknown' features:
#>  Common feature: 413 
#>  Unique to 1: 264 
#>  Unique to 2: 331 
#> 
#> [[2]]
#> Object of class "FeatComp", 'ER' features:
#>  Common feature: 21 
#>  Unique to 1: 7 
#>  Unique to 2: 0 
#> 
head(common(x[[1]]))
#> [1] "P20353" "P53501" "Q7KU78" "P04412" "Q7JZN0" "Q9VM65"

data(tan2009r3)
tanl <- list(tan2009r1, tan2009r2, tan2009r3)
xx <- compfnames(tanl, fcol1 = NULL)
length(xx)
#> [1] 15
tail(xx)
#> $`1v3`
#> Object of class "FeatComp", 'all' features:
#>  Common feature: 494 
#>  Unique to 1: 394 
#>  Unique to 2: 176 
#> 
#> $`1v2`
#> Object of class "FeatComp", 'all' features:
#>  Common feature: 551 
#>  Unique to 1: 337 
#>  Unique to 2: 320 
#> 
#> $`1v3`
#> Object of class "FeatComp", 'all' features:
#>  Common feature: 494 
#>  Unique to 1: 394 
#>  Unique to 2: 176 
#> 
#> $`3v2`
#> Object of class "FeatComp", 'all' features:
#>  Common feature: 477 
#>  Unique to 1: 193 
#>  Unique to 2: 394 
#> 
#> $`3v3`
#> Object of class "FeatComp", 'all' features:
#>  Common feature: 670 
#>  Unique to 1: 0 
#>  Unique to 2: 0 
#> 
#> $`2v3`
#> Object of class "FeatComp", 'all' features:
#>  Common feature: 477 
#>  Unique to 1: 394 
#>  Unique to 2: 193 
#> 

all.equal(xx[[15]],
          compfnames(tan2009r2, tan2009r3, fcol1 = NULL))
#> [1] TRUE
str(sapply(xx, common))
#> List of 15
#>  $ 1v2: chr [1:551] "P20353" "P53501" "Q7KU78" "P04412" ...
#>  $ 1v1: chr [1:888] "P20353" "P53501" "Q7KU78" "P04412" ...
#>  $ 1v3: chr [1:494] "Q7KU78" "Q7JZN0" "P15348" "Q00174" ...
#>  $ 1v2: chr [1:551] "P20353" "P53501" "Q7KU78" "P04412" ...
#>  $ 1v3: chr [1:494] "Q7KU78" "Q7JZN0" "P15348" "Q00174" ...
#>  $ 2v1: chr [1:551] "P20353" "P53501" "Q7KU78" "P04412" ...
#>  $ 2v3: chr [1:477] "Q7KU78" "Q7JZN0" "Q9VJ80" "P15348" ...
#>  $ 2v2: chr [1:871] "P20432" "P20353" "P53501" "Q7KU78" ...
#>  $ 2v3: chr [1:477] "Q7KU78" "Q7JZN0" "Q9VJ80" "P15348" ...
#>  $ 1v3: chr [1:494] "Q7KU78" "Q7JZN0" "P15348" "Q00174" ...
#>  $ 1v2: chr [1:551] "P20353" "P53501" "Q7KU78" "P04412" ...
#>  $ 1v3: chr [1:494] "Q7KU78" "Q7JZN0" "P15348" "Q00174" ...
#>  $ 3v2: chr [1:477] "Q7KU78" "Q7JZN0" "Q9VJ80" "P15348" ...
#>  $ 3v3: chr [1:670] "Q7KU78" "Q7JZN0" "M9PC85" "Q9VJ80" ...
#>  $ 2v3: chr [1:477] "Q7KU78" "Q7JZN0" "Q9VJ80" "P15348" ...
```
