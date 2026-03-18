# Combines features in an `MSnSet` object

This function combines the features in an
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
instance applying a summarisation function (see `fun` argument) to sets
of features as defined by a factor (see `fcol` argument). Note that the
feature names are automatically updated based on the `groupBy`
parameter.

The coefficient of variations are automatically computed and collated to
the featureData slot. See `cv` and `cv.norm` arguments for details.

If NA values are present, a message will be shown. Details on how
missing value impact on the data aggregation are provided below.

## Arguments

- object:

  An instance of class
  `"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
  whose features will be summerised.

- groupBy:

  A `factor`, `character`, `numeric` or a `list` of the above defining
  how to summerise the features. The list must be of length
  `nrow(object)`. Each element of the list is a vector describing the
  feature mapping. If the list can be named, its names must match
  `fetureNames(object)`. See `redundancy.handler` for details about the
  latter.

- fun:

  Deprecated; use `method` instead.

- method:

  The summerising function. Currently, mean, median, weighted mean, sum,
  median polish, robust summarisation (using
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html), implemented in
  [`MsCoreUtils::robustSummary()`](https://rdrr.io/pkg/MsCoreUtils/man/robustSummary.html)),
  iPQF (see [`iPQF`](https://lgatto.github.io/MSnbase/reference/iPQF.md)
  for details) and NTR (see
  [`NTR`](https://lgatto.github.io/MSnbase/reference/normToReference.md)
  for details) are implemented, but user-defined functions can also be
  supplied. Note that the robust menthods assumes that the data are
  already log-transformed.

- fcol:

  Feature meta-data label (fData column name) defining how to summerise
  the features. It must be present in `fvarLabels(object)` and, if
  present, will be used to defined `groupBy` as `fData(object)[, fcol]`.
  Note that `fcol` is ignored if `groupBy` is present.

- redundancy.handler:

  If `groupBy` is a `list`, one of `"unique"` (default) or `"multiple"`
  (ignored otherwise) defining how to handle peptides that can be
  associated to multiple higher-level features (proteins) upon
  combination. Using `"unique"` will only consider uniquely matching
  features (features matching multiple proteins will be discarded).
  `"multiple"` will allow matching to multiple proteins and each feature
  will be repeatedly tallied for each possible matching protein.

- cv:

  A `logical` defining if feature coefficients of variation should be
  computed and stored as feature meta-data. Default is `TRUE`.

- cv.norm:

  A `character` defining how to normalise the feature intensitites prior
  to CV calculation. Default is `sum`. Use `none` to keep intensities as
  is. See
  [`featureCV`](https://lgatto.github.io/MSnbase/reference/featureCV.md)
  for more details.

- verbose:

  A `logical` indicating whether verbose output is to be printed out.

- ...:

  Additional arguments for the `fun` function.

## Value

A new
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
instance is returned with `ncol` (i.e. number of samples) is unchanged,
but `nrow` (i.e. the number od features) is now equals to the number of
levels in `groupBy`. The feature metadata (`featureData` slot) is
updated accordingly and only the first occurrence of a feature in the
original feature meta-data is kept.

## Details

Missing values have different effect based on the aggregation method
employed, as detailed below. See also examples below.

1.  When using either `"sum"`, `"mean"`, `"weighted.mean"` or
    `"median"`, any missing value will be propagated at the higher
    level. If `na.rm = TRUE` is used, then the missing value will be
    ignored.

2.  Missing values will result in an error when using `"medpolish"`,
    unless `na.rm = TRUE` is used.

3.  When using robust summarisation (`"robust"`), individual missing
    values are excluded prior to fitting the linear model by robust
    regression. To remove all values in the feature containing the
    missing values, use `filterNA`.

4.  The `"iPQF"` method will fail with an error if missing value are
    present, which will have to be handled explicitly. See below.

More generally, missing values often need dedicated handling such as
filtering (see
[`filterNA`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md))
or imputation (see `impute`).

## Author

Laurent Gatto with contributions from Martina Fischer for iPQF and
Ludger Goeminne, Adriaan Sticker and Lieven Clement for robust.

## References

iPQF: a new peptide-to-protein summarization method using peptide
spectra characteristics to improve protein quantification. Fischer M,
Renard BY. Bioinformatics. 2016 Apr 1;32(7):1040-7.
doi:10.1093/bioinformatics/btv675. Epub 2015 Nov 20. PubMed
PMID:26589272.

## See also

[`featureCV`](https://lgatto.github.io/MSnbase/reference/featureCV.md)
to calculate coefficient of variation,
[`nFeatures`](https://lgatto.github.io/MSnbase/reference/nFeatures.md)
to document the number of features per group in the feature data, and
the [`aggvar`](https://lgatto.github.io/MSnbase/reference/aggvar.md) to
explore variability within protein groups.

[`iPQF`](https://lgatto.github.io/MSnbase/reference/iPQF.md) for iPQF
summarisation.

[`NTR`](https://lgatto.github.io/MSnbase/reference/normToReference.md)
for normalisation to reference summarisation.

## Examples

``` r
data(msnset)
msnset <- msnset[11:15, ]
exprs(msnset)
#>     iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> X19  32838.044  37066.058  41429.627  39700.475
#> X2    3715.089   4254.323   4748.462   5249.904
#> X20  34509.686  34928.747  41911.032  42843.839
#> X21  21262.148  23168.729  25407.068  25949.954
#> X22   8635.316  10036.529   9254.432   7769.749

## arbitrary grouping into two groups
grp <- as.factor(c(1, 1, 2, 2, 2))
msnset.comb <- combineFeatures(msnset, groupBy = grp, method = "sum")
dim(msnset.comb)
#> [1] 2 4
exprs(msnset.comb)
#>   iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> 1   36553.13   41320.38   46178.09   44950.38
#> 2   64407.15   68134.01   76572.53   76563.54
fvarLabels(msnset.comb)
#>  [1] "spectrum"            "ProteinAccession"    "ProteinDescription" 
#>  [4] "PeptideSequence"     "file"                "retention.time"     
#>  [7] "precursor.mz"        "precursor.intensity" "charge"             
#> [10] "peaks.count"         "tic"                 "ionCount"           
#> [13] "ms.level"            "acquisition.number"  "collision.energy"   
#> [16] "CV.iTRAQ4.114"       "CV.iTRAQ4.115"       "CV.iTRAQ4.116"      
#> [19] "CV.iTRAQ4.117"      

## grouping with a list
grpl <- list(c("A", "B"), "A", "A", "C", c("C", "B"))
## optional naming
names(grpl) <- featureNames(msnset)
exprs(combineFeatures(msnset, groupBy = grpl, method = "sum", redundancy.handler = "unique"))
#>   iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> A   38224.78   39183.07   46659.49   48093.74
#> C   21262.15   23168.73   25407.07   25949.95
exprs(combineFeatures(msnset, groupBy = grpl, method = "sum", redundancy.handler = "multiple"))
#>   iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> A   71062.82   76249.13   88089.12   87794.22
#> B   41473.36   47102.59   50684.06   47470.22
#> C   29897.46   33205.26   34661.50   33719.70

## missing data
exprs(msnset)[4, 4] <-
    exprs(msnset)[2, 2] <- NA
exprs(msnset)
#>     iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> X19  32838.044   37066.06  41429.627  39700.475
#> X2    3715.089         NA   4748.462   5249.904
#> X20  34509.686   34928.75  41911.032  42843.839
#> X21  21262.148   23168.73  25407.068         NA
#> X22   8635.316   10036.53   9254.432   7769.749
## NAs propagate in the 115 and 117 channels
exprs(combineFeatures(msnset, grp, "sum"))
#> Your data contains missing values. Please read the relevant section in
#> the combineFeatures manual page for details on the effects of missing
#> values on data aggregation.
#>   iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> 1   36553.13         NA   46178.09   44950.38
#> 2   64407.15   68134.01   76572.53         NA
## NAs are removed before summing
exprs(combineFeatures(msnset, grp, "sum", na.rm = TRUE))
#> Your data contains missing values. Please read the relevant section in
#> the combineFeatures manual page for details on the effects of missing
#> values on data aggregation.
#>   iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> 1   36553.13   37066.06   46178.09   44950.38
#> 2   64407.15   68134.01   76572.53   50613.59

## using iPQF
data(msnset2)
anyNA(msnset2)
#> [1] FALSE
res <- combineFeatures(msnset2,
                       groupBy = fData(msnset2)$accession,
                       redundancy.handler = "unique",
                       method = "iPQF",
                       low.support.filter = FALSE,
                       ratio.calc = "sum",
                       method.combine = FALSE)
#> The following 1 proteins are only supported by 1 or 2 peptide spectra,
#> hence, protein quantification is not reliable and can only be calculated
#> by the 'mean' in these cases, corresponding protein accessions are:
#>   O95678
head(exprs(res))
#>        X114.ions X115.ions X116.ions X117.ions
#> O95678 0.2404726 0.2682764 0.2584247 0.2328263
#> P01766 0.2610278 0.2467206 0.2544715 0.2377801
#> P01776 0.2678859 0.2591250 0.2423396 0.2306495
#> P02749 0.2640340 0.2523566 0.2510357 0.2325736
#> P02763 0.2503318 0.2524583 0.2501628 0.2470472
#> P07225 0.2533961 0.2506013 0.2504353 0.2455673

## using robust summarisation
data(msnset) ## reset data
msnset <- log(msnset, 2) ## log2 transform

## Feature X46, in the ENO protein has one missig value
which(is.na(msnset), arr.ind = TRUE)
#>     row col
#> X2    2   2
#> X21   4   4
exprs(msnset["X46", ])
#> Error in (function (cond) .Internal(C_tryCatchHelper(addr, 1L, cond)))(structure(list(message = "subscript out of bounds",     call = orig[[nm]][i, , ..., drop = drop], object = structure(c(15.0030805728121,     11.8591811816279, 15.0747137350561, 14.3759997002879, 13.0760333393114,     15.177811062346, NA, 15.0921272743665, 14.4998912898294,     13.2929728095096, 15.3383752007826, 12.2132444458243, 15.3550424297539,     14.6329423029583, 13.1759286811449, 15.2768686547183, 12.3580752059455,     15.3868001257316, NA, 12.9236522919322), dim = 5:4, dimnames = list(        c("X19", "X2", "X20", "X21", "X22"), c("iTRAQ4.114",         "iTRAQ4.115", "iTRAQ4.116", "iTRAQ4.117"))), subscript = 1L,     index = "X46"), class = c("subscriptOutOfBoundsError", "error", "condition"))): error in evaluating the argument 'object' in selecting a method for function 'exprs': subscript out of bounds
## Only the missing value in X46 and iTRAQ4.116 will be ignored
res <- combineFeatures(msnset,
                       fcol = "ProteinAccession",
                       method = "robust")
#> Your data contains missing values. Please read the relevant section in
#> the combineFeatures manual page for details on the effects of missing
#> values on data aggregation.
tail(exprs(res))
#>         iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> ECA1032   15.00308   15.17781   15.33838   15.27687
#> ECA1104   14.37600   14.49989   14.63294         NA
#> ECA1294   11.85918         NA   12.21324   12.35808
#> ECA3356   13.07603   13.29297   13.17593   12.92365
#> ECA4514   15.07471   15.09213   15.35504   15.38680

msnset2 <- filterNA(msnset) ## remove features with missing value(s)
res2 <- combineFeatures(msnset2,
                        fcol = "ProteinAccession",
                        method = "robust")
## Here, the values for ENO are different because the whole feature
## X46 that contained the missing value was removed prior to fitting.
tail(exprs(res2))
#>         iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> ECA1032   15.00308   15.17781   15.33838   15.27687
#> ECA3356   13.07603   13.29297   13.17593   12.92365
#> ECA4514   15.07471   15.09213   15.35504   15.38680
```
