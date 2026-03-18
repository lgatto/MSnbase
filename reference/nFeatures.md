# How many features in a group?

This function computes the number of features in the group defined by
the feature variable `fcol` and appends this information in the feature
data of `object`.

## Usage

``` r
nFeatures(object, fcol)
```

## Arguments

- object:

  An instance of class `MSnSet`.

- fcol:

  Feature variable defining the feature grouping structure.

## Value

An updated `MSnSet` with a new feature variable `fcol.nFeatures`.

## Author

Laurent Gatto

## Examples

``` r
library(pRolocdata)
data("hyperLOPIT2015ms3r1psm")
hyperLOPIT2015ms3r1psm <- nFeatures(hyperLOPIT2015ms3r1psm,
                                    "Protein.Group.Accessions")
i <- c("Protein.Group.Accessions", "Protein.Group.Accessions.nFeatures")
fData(hyperLOPIT2015ms3r1psm)[1:10, i]
#>     Protein.Group.Accessions Protein.Group.Accessions.nFeatures
#> 136                   Q99PL5                                108
#> 141                   P14901                                  7
#> 169                   O89090                                  3
#> 203                 O89032-3                                  1
#> 238                   Q8R2M2                                 13
#> 239                   Q9JIG7                                  3
#> 261                   Q3TWI9                                  5
#> 268                   Q8K0H5                                  2
#> 274                   Q99NH0                                  9
#> 276                   Q9CZD3                                 24
```
