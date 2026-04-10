# Features of Interest

The *Features of Interest* infrastructure allows to define a set of
features of particular interest to be used/matched against existing data
sets contained in
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`.
A specific set of features is stored as an `FeaturesOfInterest` object
and a collection of such non-redundant instances (for example for a
specific organism, project, ...) can be collected in a `FoICollection`.

## Objects from the Class

Objects can be created with the respective `FeaturesOfInterest` and
`FoICollection` constructors.

`FeaturesOfInterest` instances can be generated in two different ways:
the constructor takes either (1) a set of features names (a `character`
vector) and a description (`character` of length 1 - any subsequent
elements are silently ignored) or (2) feature names, a description and
an instance of class
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`.
In the latter case, we call such `FeaturesOfInterest` objects traceable,
because we can identify the origin of the feature names and thus their
validity. This is done by inspecting the `MSnSet` instance and recording
its dimensions, its name and a unique md5 hash tag (these are stores as
part of the optional `objpar` slot). In such cases, the feature names
passed to the `FeaturesOfInterest` constructor must also be present in
the `MSnSet`; if one or more are not, an error will be thrown. If your
features of interest to be recorded stem for an existing experiment and
have all been observed, it is advised to pass the 3 arguments to the
constructor to ensure that the feature names as valid. Otherwise, only
the third argument should be omitted.

`FoICollection` instances can be constructed by creating an empty
collection and serial additions of `FeaturesOfInterest` using
`addFeaturesOfInterest` or by passing a list of `FeaturesOfInterest`
instance.

## Slots

`FeaturesOfInterest` class:

- `description`::

  Object of class `"character"` describing the instance.

- `objpar`::

  Optional object of class `"list"` providing details about the `MSnSet`
  instance originally used to create the instance. See details section.

- `fnames`::

  Object of class `"character"` with the feature of interest names.

- `date`::

  Object of class `"character"` with the date the instance was first
  generated.

- `.__classVersion__`::

  Object of class `"Versions"` with the `FeaturesOfInterest` class
  version. Only relevant for development.

`FoICollection` class:

- `foic`::

  Object of class `"list"` with the `FeaturesOfInterest`.

- `.__classVersion__`::

  Object of class `"Versions"` with the `FoICollection` class version.
  Only relevant for development.

## Extends

Class `"Versioned"`, directly.

## Methods

`FeaturesOfInterest` class:

- description:

  `signature(object = "FeaturesOfInterest")`: returns the description of
  `object`.

- foi:

  `signature(object = "FeaturesOfInterest")`: returns the features of
  interests.

- length:

  `signature(x = "FeaturesOfInterest")`: returns the number of features
  of interest in `x`.

- show:

  `signature(object = "FeaturesOfInterest")`: displays `object`.

- fnamesIn:

  `signature(x = "FeaturesOfInterst", y = "MSnSet", count = "logical")`:
  if `count` is `FALSE` (default), return a logical indicating whether
  there is at least one feautre of interest present in `x`? Otherwise,
  returns the number of such features. Works also with matrices and
  data.frames.

- \[:

  Subsetting works like lists. Returns a new `FoICollection`.

- \[\[:

  Subsetting works like lists. Returns a new `FeatureOfInterest`.

`FoICollection` class:

- description:

  `signature(object = "FoICollection")`: returns the description of
  `object`.

- foi:

  `signature(object = "FoICollection")`: returns a list of
  `FeaturesOfInterest`.

- length:

  `signature(x = "FoICollection")`: returns the number of
  `FeaturesOfInterest` in the collection.

- lengths:

  `signature(x = "FoICollection")`: returns the number of features of
  interest in each `FeaturesOfInterest` in the collection `x`.

- addFeaturesOfInterest:

  `signature(x = "FeaturesOfInterest", y = "FoICollection")`: add the
  `FeaturesOfInterest` instance `x` to `FoICollection` `y`. If `x` is
  already present, a message is printed and `y` is returned unchanged.

- rmFeaturesOfInterest:

  `signature(object = "FoICollection", i = "numeric")`: removes the
  `i`th `FeatureOfInterest` in the collection `object`.

- show:

  `signature(object = "FoICollection")`: displays `object`.

## Author

Laurent Gatto

## Examples

``` r
library("pRolocdata")
data(tan2009r1)

x <- FeaturesOfInterest(description = "A traceable test set of features of interest",
                        fnames = featureNames(tan2009r1)[1:10],
                        object = tan2009r1)
x
#> Traceable object of class "FeaturesOfInterest"
#>  Created on Fri Apr 10 14:42:24 2026 
#>  Description:
#>   A traceable test set of features of interest
#>  10 features of interest:
#>    P20353, P53501  ...  Q9VCK0, Q9VIU7

description(x)
#> [1] "A traceable test set of features of interest"
foi(x)
#>  [1] "P20353" "P53501" "Q7KU78" "P04412" "Q7KJ73" "Q7JZN0" "Q7KLV9" "Q9VM65"
#>  [9] "Q9VCK0" "Q9VIU7"

y <- FeaturesOfInterest(description = "Non-traceable features of interest",
                        fnames = featureNames(tan2009r1)[111:113])
y
#> Object of class "FeaturesOfInterest"
#>  Created on Fri Apr 10 14:42:24 2026 
#>  Description:
#>   Non-traceable features of interest
#>  3 features of interest:
#>    Q9VT75, Q9VNA3, A8JNJ6

## an illegal FeaturesOfInterest
try(FeaturesOfInterest(description = "Won't work",
                       fnames = c("A", "Z", featureNames(tan2009r1)),
                       object = tan2009r1))
#> Error in FeaturesOfInterest(description = "Won't work", fnames = c("A",  : 
#>   2 feature(s) of interest absent from your object's feature names:
#>    A, Z.


FeaturesOfInterest(description = "This work, but not traceable",
                       fnames = c("A", "Z", featureNames(tan2009r1)))
#> Object of class "FeaturesOfInterest"
#>  Created on Fri Apr 10 14:42:24 2026 
#>  Description:
#>   This work, but not traceable
#>  890 features of interest:
#>    A, Z  ...  Q8SZM1, P07909



xx <- FoICollection()
xx
#> A collection of 0 features of interest.

xx <- addFeaturesOfInterest(x, xx)
xx <- addFeaturesOfInterest(y, xx)
names(xx) <- LETTERS[1:2]
xx
#> A collection of 2 features of interest.

## Sub-setting
xx[1]
#> A collection of 1 features of interest.
xx[[1]]
#> Traceable object of class "FeaturesOfInterest"
#>  Created on Fri Apr 10 14:42:24 2026 
#>  Description:
#>   A traceable test set of features of interest
#>  10 features of interest:
#>    P20353, P53501  ...  Q9VCK0, Q9VIU7
xx[["A"]]
#> Traceable object of class "FeaturesOfInterest"
#>  Created on Fri Apr 10 14:42:24 2026 
#>  Description:
#>   A traceable test set of features of interest
#>  10 features of interest:
#>    P20353, P53501  ...  Q9VCK0, Q9VIU7

description(xx)
#>                                              A 
#> "A traceable test set of features of interest" 
#>                                              B 
#>           "Non-traceable features of interest" 
foi(xx)
#> $A
#> Traceable object of class "FeaturesOfInterest"
#>  Created on Fri Apr 10 14:42:24 2026 
#>  Description:
#>   A traceable test set of features of interest
#>  10 features of interest:
#>    P20353, P53501  ...  Q9VCK0, Q9VIU7
#> 
#> $B
#> Object of class "FeaturesOfInterest"
#>  Created on Fri Apr 10 14:42:24 2026 
#>  Description:
#>   Non-traceable features of interest
#>  3 features of interest:
#>    Q9VT75, Q9VNA3, A8JNJ6
#> 

fnamesIn(x, tan2009r1)
#> [1] TRUE
fnamesIn(x, tan2009r1, count = TRUE)
#> [1] 10

rmFeaturesOfInterest(xx, 1)
#> A collection of 1 features of interest.
```
