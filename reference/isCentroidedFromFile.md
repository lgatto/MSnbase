# Get mode from mzML data file

The function extracts the mode (profile or centroided) from the raw mass
spectrometry file by parsing the mzML file directly. If the object `x`
stems from any other type of file, `NA`s are returned.

## Usage

``` r
isCentroidedFromFile(x)
```

## Arguments

- x:

  An object of class
  [OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md).

## Value

A named `logical` vector of the same length as `x`.

## Details

This function is much faster than
[`isCentroided()`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md),
which estimates mode from the data, but is limited to data stemming from
mzML files which are still available in their original location (and
accessed with `fileNames(x)`).

## Author

Laurent Gatto

## Examples

``` r
library("msdata")
f <- proteomics(full.names = TRUE,
                pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
x <- readMSData(f, mode = "onDisk")
table(isCentroidedFromFile(x), msLevel(x))
#>        
#>           1   2
#>   FALSE  58   0
#>   TRUE    0 451
```
