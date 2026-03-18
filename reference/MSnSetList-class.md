# Storing multiple related MSnSets

A class for storing lists of
[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
instances.

## Details

There are two ways to store different sets of measurements pertaining an
experimental unit, such as replicated measures of different conditions
that were recorded over more than one MS acquisition. Without focusing
on any proteomics technology in particular, these multiple assays can be
recorded as

- A single combined `MSnSet` (see the section *Combining MSnSet
  instances* in the *MSnbase-demo* section). In such cases, the
  different experimental (phenotypical) conditions are recorded as an
  [`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)
  in the `phenoData` slots.

  Quantitative data for features that were missing in an assay are
  generally encode as missing with `NA` values. Alternatively, only
  features observed in all assays could be selected. See the
  [`commonFeatureNames`](https://lgatto.github.io/MSnbase/reference/commonFeatureNames.md)
  functions to select only common features among two or more `MSnSet`
  instance.

- Each set of measurements is stored in an `MSnSet` which are combined
  into one `MSnSetList`. Each `MSnSet` elements can have identical or
  different samples and features. Unless compiled directly manually by
  the user, one would expect at least one of these dimensions
  (features/rows or samples/columns) are conserved (i.e. all feature or
  samples names are identical). See `split`/`unsplit` below.

## Objects from the Class

Objects can be created and manipluated with:

- `MSnSetList(x, log, featureDAta)`:

  The class constructor that takes a list of valid `MSnSet` instances as
  input `x`, an optional logging `list`, and an optional feature
  metadata `data.frame`.

- `split(x, f)`:

  An `MSnSetList` can be created from an
  [`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
  instance. `x` is a single `MSnSet` and `f` is a `factor` or a
  `character` of length 1. In the latter case, `f` will be matched to
  the feature- and phenodata variable names (in that order). If a match
  is found, the respective variable is extracted, converted to a factor
  if it is not one already, and used to split `x` along the
  features/rows (`f` was a feature variable name) or samples/columns
  (`f` was a phenotypic variable name). If `f` is passed as a factor,
  its length will be matched to `nrow(x)` or `ncol(x)` (in that order)
  to determine if `x` will be split along the features (rows) or sample
  (columns). Hence, the length of `f` must match exactly to either
  dimension.

- `unsplit(value, f)`:

  The `unsplit` method reverses the effect of splitting the `value`
  `MSnSet` along the groups `f`.

- `as(x, "MSnSetList")`:

  Where `x` is an instance of class
  [MzTab](https://lgatto.github.io/MSnbase/reference/MzTab-class.md).
  See the class documentation for details.

## Slots

- `x`::

  Object of class `list` containing valid `MSnSet` instances. Can be
  extracted with the `msnsets()` accessor.

- `log`::

  Object of class `list` containing an object creation log, containing
  among other elements the call that generated the object. Can be
  accessed with `objlog()`.

- `featureData`::

  Object of class `DataFrame` that stores metadata for each object in
  the `x` slot. The number of rows of this `data.frame` must be equal to
  the number of items in the `x` slot and their respective (row)names
  must be identical.

- `.__classVersion__`::

  The version of the instance. For development purposes only.

## Methods

- `"[["`:

  Extracts a single `MSnSet` at position.

- `"["`:

  Extracts one of more `MSnSets` as `MSnSetList`.

- `length`:

  Returns the number of `MSnSets`.

- `names`:

  Returns the names of `MSnSets`, if available. The replacement method
  is also available.

- `show`:

  Display the object by printing a short summary.

- `lapply(x, FUN, ...)`:

  Apply function `FUN` to each element of the input `x`. If the
  application of `FUN` returns and `MSnSet`, then the return value is an
  `MSnSetList`, otherwise a `list`

- `sapply(x, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)`:

  A `lapply` wrapper that simplifies the ouptut to a vector, matric or
  array is possible. See
  [`?base::sapply`](https://rdrr.io/r/base/lapply.html) for details.

- `fData`:

  Returns the features metadata `featureData` slot.

- `fData<-`:

  Features metadata `featureData` replacement method.

## See also

The
[`commonFeatureNames`](https://lgatto.github.io/MSnbase/reference/commonFeatureNames.md)
function to select common features among `MSnSet` instances.

## Author

Laurent Gatto

## Examples

``` r
library("pRolocdata")
data(tan2009r1)
data(tan2009r2)

## The MSnSetList class
##  for an unnamed list, names are set to indices
msnl <- MSnSetList(list(tan2009r1, tan2009r2))
names(msnl)
#> [1] "1" "2"
##  a named example
msnl <- MSnSetList(list(A = tan2009r1, B = tan2009r2))
names(msnl)
#> [1] "A" "B"
msnsets(msnl)
#> $A
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 888 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: X114 X115 X116 X117
#>   varLabels: Fractions
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: P20353 P53501 ... P07909 (888 total)
#>   fvarLabels: FBgn Protein.ID ... markers.tl (16 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#>   pubMedIds: 19317464 
#> Annotation:  
#> - - - Processing information - - -
#> Added markers from  'mrk' marker vector. Thu Jul 16 22:53:44 2015 
#>  MSnbase version: 1.17.12 
#> 
#> $B
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 871 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: X114 X115 X116 X117
#>   varLabels: Fractions
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: P20432 P20353 ... Q9VIW3 (871 total)
#>   fvarLabels: FBgn Protein.ID ... markers (13 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#>   pubMedIds: 19317464 
#> Annotation:  
#> - - - Processing information - - -
#> Added markers from  'mrk' marker vector. Thu Jul 16 22:53:44 2015 
#>  MSnbase version: 1.17.12 
#> 
length(msnl)
#> [1] 2
objlog(msnl)
#> $call
#> MSnSetList(x = list(A = tan2009r1, B = tan2009r2))
#> 
msnl[[1]] ## an MSnSet
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 888 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: X114 X115 X116 X117
#>   varLabels: Fractions
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: P20353 P53501 ... P07909 (888 total)
#>   fvarLabels: FBgn Protein.ID ... markers.tl (16 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#>   pubMedIds: 19317464 
#> Annotation:  
#> - - - Processing information - - -
#> Added markers from  'mrk' marker vector. Thu Jul 16 22:53:44 2015 
#>  MSnbase version: 1.17.12 
msnl[1]   ## an MSnSetList of length 1
#> Instance of class 'MSnSetList' containig 1 objects.

## Iterating over the elements
lapply(msnl, dim) ## a list
#> $A
#> [1] 888   4
#> 
#> $B
#> [1] 871   4
#> 
lapply(msnl, normalise, method = "quantiles") ## an MSnSetList
#> Instance of class 'MSnSetList' containig 2 objects.

fData(msnl)
#> DataFrame with 2 rows and 0 columns
fData(msnl)$X <- sapply(msnl, nrow)
fData(msnl)
#> DataFrame with 2 rows and 1 column
#>           X
#>   <integer>
#> A       888
#> B       871

## Splitting and unsplitting
##  splitting along the columns/samples
data(dunkley2006)
head(pData(dunkley2006))
#>        membrane.prep fraction replicate
#> M1F1A              1        1         A
#> M1F4A              1        4         A
#> M1F7A              1        7         A
#> M1F11A             1       11         A
#> M1F2B              1        2         B
#> M1F5B              1        5         B
(splt <- split(dunkley2006, "replicate"))
#> Instance of class 'MSnSetList' containig 2 objects.
lapply(splt, dim) ## the number of rows and columns of the split elements
#> $A
#> [1] 689   8
#> 
#> $B
#> [1] 689   8
#> 
unsplt <- unsplit(splt, dunkley2006$replicate)
stopifnot(compareMSnSets(dunkley2006, unsplt))

##  splitting along the rows/features
head(fData(dunkley2006))
#>           assigned  evidence method   new pd.2013 pd.markers markers.orig
#> AT1G09210       ER predicted  PLSDA known      ER   ER lumen           ER
#> AT1G21750       ER predicted  PLSDA known      ER   ER lumen           ER
#> AT1G51760       ER   unknown  PLSDA   new      ER   ER lumen      unknown
#> AT1G56340       ER predicted  PLSDA known      ER   ER lumen           ER
#> AT2G32920       ER predicted  PLSDA known      ER   ER lumen           ER
#> AT2G47470       ER predicted  PLSDA known      ER   ER lumen           ER
#>            markers
#> AT1G09210 ER lumen
#> AT1G21750 ER lumen
#> AT1G51760 ER lumen
#> AT1G56340 ER lumen
#> AT2G32920 ER lumen
#> AT2G47470 ER lumen
(splt <- split(dunkley2006, "markers"))
#> Instance of class 'MSnSetList' containig 10 objects.
unsplt <- unsplit(splt, factor(fData(dunkley2006)$markers))
simplify2array(lapply(splt, dim))
#>      ER lumen ER membrane Golgi Mitochondrion PM Plastid Ribosome TGN unknown
#> [1,]       14          45    28            55 46      20       19  13     428
#> [2,]       16          16    16            16 16      16       16  16      16
#>      vacuole
#> [1,]      21
#> [2,]      16
stopifnot(compareMSnSets(dunkley2006, unsplt))
```
