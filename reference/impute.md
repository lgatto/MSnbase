# Quantitative proteomics data imputation

The `impute` method performs data imputation on `MSnSet` instances using
a variety of methods.

Users should proceed with care when imputing data and take precautions
to assure that the imputation produce valid results, in particular with
naive imputations such as replacing missing values with 0.

See
[`MsCoreUtils::impute_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
for details on the different imputation methods available and
strategies.

## Usage

``` r
# S4 method for class 'MSnSet'
impute(object, method, ...)
```

## Arguments

- object:

  An `MSnSet` object with missing values to be imputed.

- method:

  `character(1)` defining the imputation method. See
  [`MsCoreUtils::imputeMethods()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
  for available ones. See
  [`MsCoreUtils::impute_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
  for details.

- ...:

  Additional parameters passed to the inner imputation function. See
  [`MsCoreUtils::impute_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
  for details.

## Examples

``` r

data(naset)

## table of missing values along the rows
table(fData(naset)$nNA)
#> 
#>   0   1   2   3   4   8   9  10 
#> 301 247  91  13   2  23  10   2 

## table of missing values along the columns
pData(naset)$nNA
#>  [1] 34 45 56 39 47 52 49 61 41 42 55 45 51 43 57 53

## non-random missing values
notna <- which(!fData(naset)$randna)
length(notna)
#> [1] 35
notna
#>  [1]   6  20  79  88 130 187 227 231 238 264 275 317 324 363 373 382 409 437 445
#> [20] 453 456 474 484 485 492 514 516 546 568 580 594 631 648 664 671

impute(naset, method = "min")
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 689 features, 16 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: M1F1A M1F4A ... M2F11B (16 total)
#>   varLabels: nNA
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: AT1G09210 AT1G21750 ... AT4G39080 (689 total)
#>   fvarLabels: nNA randna
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> Data imputation using min Fri Apr 10 14:44:29 2026 
#>  MSnbase version: 1.15.6 

if (require("imputeLCMD")) {
    impute(naset, method = "QRILC")
    impute(naset, method = "MinDet")
}
#> Loading required package: imputeLCMD
#> Loading required package: tmvtnorm
#> Loading required package: mvtnorm
#> Loading required package: Matrix
#> 
#> Attaching package: ‘Matrix’
#> The following object is masked from ‘package:S4Vectors’:
#> 
#>     expand
#> Loading required package: gmm
#> Loading required package: sandwich
#> 
#> Attaching package: ‘sandwich’
#> The following object is masked from ‘package:generics’:
#> 
#>     estfun
#> Loading required package: norm
#> This package has some major limitations
#> (for example, it does not work reliably when
#> the number of variables exceeds 30),
#> and has been superseded by the norm2 package.
#> Loading required package: pcaMethods
#> 
#> Attaching package: ‘pcaMethods’
#> The following object is masked from ‘package:stats’:
#> 
#>     loadings
#> Loading required package: impute
#> Imputing along margin 2 (samples/columns).
#> Imputing along margin 2 (samples/columns).
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 689 features, 16 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: M1F1A M1F4A ... M2F11B (16 total)
#>   varLabels: nNA
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: AT1G09210 AT1G21750 ... AT4G39080 (689 total)
#>   fvarLabels: nNA randna
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> Data imputation using MinDet Fri Apr 10 14:44:30 2026 
#>  MSnbase version: 1.15.6 

if (require("norm"))
    impute(naset, method = "MLE")
#> Imputing along margin 2 (samples/columns).
#> Warning: NAs introduced by coercion to integer range
#> Iterations of EM: 
#> 1...2...
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 689 features, 16 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: M1F1A M1F4A ... M2F11B (16 total)
#>   varLabels: nNA
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: AT1G09210 AT1G21750 ... AT4G39080 (689 total)
#>   fvarLabels: nNA randna
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> Data imputation using MLE Fri Apr 10 14:44:30 2026 
#>  MSnbase version: 1.15.6 

impute(naset, "mixed",
       randna = fData(naset)$randna,
       mar = "knn", mnar = "QRILC")
#> Imputing along margin 1 (features/rows).
#> Imputing along margin 1 (features/rows).
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 689 features, 16 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: M1F1A M1F4A ... M2F11B (16 total)
#>   varLabels: nNA
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: AT1G09210 AT1G21750 ... AT4G39080 (689 total)
#>   fvarLabels: nNA randna
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> Data imputation using mixed Fri Apr 10 14:44:30 2026 
#>  MSnbase version: 1.15.6 


## neighbour averaging
x <- naset[1:4, 1:6]

exprs(x)[1, 1] <- NA ## min value
exprs(x)[2, 3] <- NA ## average
exprs(x)[3, 1:2] <- NA ## min value and average
## 4th row: no imputation
exprs(x)
#>              M1F1A    M1F4A   M1F7A  M1F11A    M1F2B    M1F5B
#> AT1G09210       NA 0.275500 0.21600 0.18525 0.465667 0.199667
#> AT1G21750 0.332000 0.279667      NA 0.16600 0.451500 0.200375
#> AT1G51760       NA       NA 0.16825 0.18825 0.459750 0.214500
#> AT1G56340 0.336733       NA      NA      NA 0.487167 0.201833

exprs(impute(x, "nbavg"))
#> Assuming values are ordered.
#> Imputing along margin 1 (features/rows).
#>              M1F1A    M1F4A     M1F7A  M1F11A    M1F2B    M1F5B
#> AT1G09210 0.166000 0.275500 0.2160000 0.18525 0.465667 0.199667
#> AT1G21750 0.332000 0.279667 0.2228335 0.16600 0.451500 0.200375
#> AT1G51760 0.166000 0.167125 0.1682500 0.18825 0.459750 0.214500
#> AT1G56340 0.336733       NA        NA      NA 0.487167 0.201833
```
