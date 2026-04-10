# Imports mass-spectrometry raw data files as 'MSnExp' instances.

Reads as set of XML-based mass-spectrometry data files and generates an
[MSnExp](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
object. This function uses the functionality provided by the `mzR`
package to access data and meta data in `mzData`, `mzXML` and `mzML`.

## Usage

``` r
readMSData(
  files,
  pdata = NULL,
  msLevel. = NULL,
  verbose = isMSnbaseVerbose(),
  centroided. = NA,
  smoothed. = NA,
  cache. = 1L,
  mode = c("inMemory", "onDisk")
)
```

## Arguments

- files:

  A `character` with file names to be read and parsed.

- pdata:

  An object of class Biobase::AnnotatedDataFrame or `NULL` (default).

- msLevel.:

  MS level spectra to be read. In `inMemory` mode, use `1` for MS1
  spectra or any larger numeric for MSn spectra. Default is `2` for
  `InMemory` mode. `onDisk` mode supports multiple levels and will, by
  default, read all the data.

- verbose:

  Verbosity flag. Default is to use
  [`isMSnbaseVerbose()`](https://lgatto.github.io/MSnbase/reference/MSnbaseOptions.md).

- centroided.:

  A `logical`, indicating whether spectra are centroided or not. Default
  is `NA` in which case the information is extracted from the raw file
  (for mzML or mzXML files). In `onDisk`, it can also be set for
  different MS levels by a vector of logicals, where the first element
  is for MS1, the second element is for MS2, ... See
  [OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
  for an example.

- smoothed.:

  A `logical` indicating whether spectra already smoothed or not.
  Default is `NA`.

- cache.:

  Numeric indicating caching level. Default is 0 for MS1 and 1 MS2 (or
  higher). Only relevant for `inMemory` mode.

- mode:

  On of `"inMemory"` (default) or `"onDisk"`. The former loads the raw
  data in memory, while the latter only generates the object and the raw
  data is accessed on disk when needed. See the *benchmarking* vignette
  for memory and speed implications.

## Value

An [MSnExp](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
object for `inMemory` mode and a
[OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
object for `onDisk` mode.

## Details

When using the `inMemory` mode, the whole MS data is read from file and
kept in memory as
[Spectrum](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
objects within the
[MSnExp](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)'es
`assayData` slot.

To reduce the memory footpring especially for large MS1 data sets it is
also possible to read only selected information from the MS files and
fetch the actual spectrum data (i.e. the M/Z and intensity values) only
on demand from the original data files. This can be achieved by setting
`mode = "onDisk"`. The function returns then an
[OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
object instead of a
[MSnExp](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
object.

## Note

`readMSData` uses `normalizePath` to replace relative with absolute file
paths.

## See also

[`readMgfData()`](https://lgatto.github.io/MSnbase/reference/readMgfData.md)
to read `mgf` peak lists.

## Author

Laurent Gatto

## Examples

``` r
file <- dir(system.file(package = "MSnbase", dir = "extdata"),
            full.name = TRUE,
            pattern = "mzXML$")
mem <- readMSData(file, mode = "inMemory")
mem
#> MSn experiment data ("MSnExp")
#> Object size in memory: 0.18 Mb
#> - - - Spectra data - - -
#>  MS level(s): 2 
#>  Number of spectra: 5 
#>  MSn retention times: 25:01 - 25:02 minutes
#> - - - Processing information - - -
#> Data loaded: Fri Apr 10 15:50:47 2026 
#>  MSnbase version: 2.37.3 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: dummyiTRAQ.mzXML
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   dummyiTRAQ.mzXML 
#> protocolData: none
#> featureData
#>   featureNames: F1.S1 F1.S2 ... F1.S5 (5 total)
#>   fvarLabels: spectrum
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
dsk <- readMSData(file, mode = "onDisk")
dsk
#> MSn experiment data ("OnDiskMSnExp")
#> Object size in memory: 0.03 Mb
#> - - - Spectra data - - -
#>  MS level(s): 2 
#>  Number of spectra: 5 
#>  MSn retention times: 25:01 - 25:02 minutes
#> - - - Processing information - - -
#> Data loaded [Fri Apr 10 15:50:49 2026] 
#>  MSnbase version: 2.37.3 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: dummyiTRAQ.mzXML
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   dummyiTRAQ.mzXML 
#> protocolData: none
#> featureData
#>   featureNames: F1.S1 F1.S2 ... F1.S5 (5 total)
#>   fvarLabels: fileIdx spIdx ... spectrum (36 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
```
