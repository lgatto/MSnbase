# Import mgf files as 'MSnExp' instances.

Reads a mgf file and generates an
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
object.

## Usage

``` r
readMgfData(filename, pdata = NULL, centroided = TRUE, smoothed = FALSE,
verbose = isMSnbaseVerbose(), cache = 1)
```

## Arguments

- filename:

  character vector with file name to be read.

- pdata:

  an object of class
  `"`[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`"`.

- smoothed:

  `Logical` indicating whether spectra already smoothed or not. Default
  is 'FALSE'. Used to initialise
  `"`[`MSnProcess`](https://lgatto.github.io/MSnbase/reference/MSnProcess-class.md)`"`
  object in `processingData` slot.

- centroided:

  `Logical` indicating whether spectra are centroided or not. Default is
  'TRUE'. Used to initialise
  `"`[`MSnProcess`](https://lgatto.github.io/MSnbase/reference/MSnProcess-class.md)`"`
  object in `processingData` slot.

- cache:

  Numeric indicating caching level. Default is 1. Under development.

- verbose:

  verbosity flag.

## Value

An instance of

## Author

Guangchuang Yu and Laurent Gatto

## Details

Note that when reading an mgf file, the original order of the spectra is
lost. Thus, if the data was originally written to mgf from an `MSnExp`
object using `writeMgfData`, although the feature names will be
identical, the spectra are not as a result of the reordering. See
example below.

## See also

[`writeMgfData`](https://lgatto.github.io/MSnbase/reference/writeMgfData-methods.md)
method to write the content of
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
or
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
objects to mgf files. Raw data files can also be read with the
[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)
function.

## Examples

``` r
  data(itraqdata)
  writeMgfData(itraqdata, con="itraqdata.mgf", COM="MSnbase itraqdata")
  itraqdata2 <- readMgfData("itraqdata.mgf")
  ## note that the order of the spectra is altered
  ## and precision of some values (precursorMz for instance)
  match(signif(precursorMz(itraqdata2),4),signif(precursorMz(itraqdata),4))
#>  [1]  1 10 11 12 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31
#> [26] 32 33 34 35 36 37 38 39  4 40 41 42 43 44 45 46 47 48 49  5 50 51 52 53 54
#> [51] 55  6  7  8  9
  ## [1]  1 10 11 12 13 14 15 16 17 18 ...
  ## ... but all the precursors are there
  all.equal(sort(precursorMz(itraqdata2)),
            sort(precursorMz(itraqdata)),
            check.attributes=FALSE,
            tolerance=10e-5)
#> [1] TRUE
  ## is TRUE
  all.equal(as.data.frame(itraqdata2[[1]]),as.data.frame(itraqdata[[1]]))
#> [1] TRUE
  ## is TRUE
  all.equal(as.data.frame(itraqdata2[[3]]),as.data.frame(itraqdata[[11]]))
#> [1] TRUE
  ## is TRUE
  f <- dir(system.file(package="MSnbase",dir="extdata"),
           full.name=TRUE,
           pattern="test.mgf")
  (x <- readMgfData(f))
#> MSn experiment data ("MSnExp")
#> Object size in memory: 0.01 Mb
#> - - - Spectra data - - -
#>  MS level(s): 2 
#>  Number of spectra: 3 
#>  MSn retention times: 17:08 - 18:47 minutes
#> - - - Processing information - - -
#> Data loaded: Wed Mar 18 16:57:27 2026 
#>  MSnbase version: 2.37.1 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: 1
#>   varLabels: sampleNames fileNumbers
#>   varMetadata: labelDescription
#> Loaded from:
#>   test.mgf 
#> protocolData: none
#> featureData
#>   featureNames: X1 X2 X3
#>   fvarLabels: TITLE PEPMASS ... SCANS (5 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
  x[[2]]
#> Object of class "Spectrum2"
#>  Precursor: 787.8283 
#>  Retention time: 18:37 
#>  Charge: 2 
#>  MSn level: 2 
#>  Peaks count: 21 
#>  Total ion count: 124558.3 
  precursorMz(x[[2]])
#> [1] 787.8283
  precursorIntensity(x[[2]])
#> [1] 880650.4
  precursorMz(x[[1]])
#> [1] 816.3383
  precursorIntensity(x[[1]]) ## was not in test.mgf
#> [1] 0
  scanIndex(x)
#>   X1   X2   X3 
#> 2162 2406 2432 
```
