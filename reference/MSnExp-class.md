# The 'MSnExp' Class for MS Data And Meta-Data

The `MSnExp` class encapsulates data and meta-data for mass spectrometry
experiments, as described in the `slots` section. Several data files
(currently in `mzXML`) can be loaded together with the function
[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md).

This class extends the virtual
`"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`
class.

In version 1.19.12, the `polarity` slot had been added to the
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
class (previously in
`"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`).
Hence, `"MSnExp"` objects created prior to this change will not be valid
anymore, since all MS2 spectra will be missing the `polarity` slot.
Object can be appropriately updated using the `updateObject` method.

The feature variables in the feature data slot will depend on the file.
See also the documentation in the `mzR` package that parses the raw data
files and produces these data.

## Objects from the Class

Objects can be created by calls of the form `new("MSnExp",...)`.
However, it is preferred to use the
[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)
function that will read raw mass spectrometry data to generate a valid
`"MSnExp"` instance.

## Slots

- `assayData`::

  Object of class `"environment"` containing the MS spectra (see
  `"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`
  and
  `"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`).
  Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `phenoData`::

  Object of class `"AnnotatedDataFrame"` containing
  experimenter-supplied variables describing sample (i.e the individual
  tags for an labelled MS experiment) See `phenoData` for more details.
  Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `featureData`::

  Object of class `"AnnotatedDataFrame"` containing variables describing
  features (spectra in our case), e.g. identificaiton data, peptide
  sequence, identification score,... (inherited from `"eSet"`). See
  `featureData` for more details. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `experimentData`::

  Object of class
  `"`[`MIAPE`](https://lgatto.github.io/MSnbase/reference/MIAPE-class.md)`"`,
  containing details of experimental methods. See `experimentData` for
  more details. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `protocolData`::

  Object of class `"AnnotatedDataFrame"` containing equipment-generated
  variables (inherited from `"eSet"`). See `protocolData` for more
  details. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `processingData`::

  Object of class
  `"`[`MSnProcess`](https://lgatto.github.io/MSnbase/reference/MSnProcess-class.md)`"`
  that records all processing. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `.__classVersion__`::

  Object of class `"Versions"` describing the versions of R, the Biobase
  package,
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`
  and `MSnExp` of the current instance. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.
  Intended for developer use and debugging (inherited from `"eSet"`).

## Extends

Class
`"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`,
directly. Class `"VersionedBiobase"`, by class "pSet", distance 2. Class
`"Versioned"`, by class "pSet", distance 3.

## Methods

See the
`"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`
class for documentation on accessors inherited from `pSet`, subsetting
and general attribute accession.

- bin:

  `signature(object = "MSnExp")`: Bins spectra. See
  [`bin`](https://lgatto.github.io/MSnbase/reference/bin-methods.md)
  documentation for more details and examples.

- clean:

  `signature(object = "MSnExp")`: Removes unused 0 intensity data
  points. See
  [`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)
  documentation for more details and examples.

- compareSpectra:

  `signature(x = "Spectrum", y = "missing")`: Compares spectra. See
  [`compareSpectra`](https://lgatto.github.io/MSnbase/reference/compareSpectra-methods.md)
  documentation for more details and examples.

- extractPrecSpectra:

  `signature(object = "MSnExp", prec = "numeric")`: extracts spectra
  with precursor MZ value equal to `prec` and returns an object of class
  'MSnExp'. See
  [`extractPrecSpectra`](https://lgatto.github.io/MSnbase/reference/extractPrecSpectra-methods.md)
  documentation for more details and examples.

- pickPeaks:

  `signature(object = "MSnExp")`: Performs the peak picking to generate
  centroided spectra. Parameter `msLevel.` allows to restrict peak
  picking to spectra of certain MS level(s). See
  [`pickPeaks`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md)
  documentation for more details and examples.

- estimateNoise:

  `signature(object = "MSnExp")`: Estimates the noise in all profile
  spectra of `object`. See
  [`estimateNoise`](https://lgatto.github.io/MSnbase/reference/estimateNoise-method.md)
  documentation for more details and examples.

- plot:

  `signature(x = "MSnExp", y = "missing")`: Plots the `MSnExp` instance.
  See
  [`plot.MSnExp`](https://lgatto.github.io/MSnbase/reference/plot-methods.md)
  documentation for more details.

- plot2d:

  `signature(object = "MSnExp", ...)`: Plots retention time against
  precursor MZ for `MSnExp` instances. See
  [`plot2d`](https://lgatto.github.io/MSnbase/reference/plot2d-methods.md)
  documentation for more details.

- plotDensity:

  `signature(object = "MSnExp", ...)`: Plots the density of parameters
  of interest. instances. See
  [`plotDensity`](https://lgatto.github.io/MSnbase/reference/plotDensity-methods.md)
  documentation for more details.

- plotMzDelta:

  `signature(object = "MSnExp", ...)`: Plots a histogram of the m/z
  difference betwee all of the highest peaks of all MS2 spectra of an
  experiment. See
  [`plotMzDelta`](https://lgatto.github.io/MSnbase/reference/plotMzDelta-methods.md)
  documentation for more details.

- quantify:

  `signature(object = "MSnExp")`: Performs quantification for all the
  MS2 spectra of the `MSnExp` instance. See
  [`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)
  documentation for more details. Also for `OnDiskMSnExp` objects.

- removePeaks:

  `signature(object = "MSnExp")`: Removes peaks lower that a threshold
  `t`. See
  [`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
  documentation for more details and examples.

- removeReporters:

  `signature(object = "MSnExp", ...)`: Removes reporter ion peaks from
  all MS2 spectra of an experiment. See
  [`removeReporters`](https://lgatto.github.io/MSnbase/reference/removeReporters-methods.md)
  documentation for more details and examples.

- smooth:

  `signature(x = "MSnExp")`: Smooths spectra. See
  [`smooth`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md)
  documentation for more details and examples.

- addIdentificationData:

  `signature(object = "MSnExp", ...)`: Adds identification data to an
  experiment. See
  [`addIdentificationData`](https://lgatto.github.io/MSnbase/reference/addIdentificationData-methods.md)
  documentation for more details and examples.

- removeNoId:

  `signature(object = "MSnExp", fcol = "pepseq", keep = NULL)`: Removes
  non-identified features. See
  [`removeNoId`](https://lgatto.github.io/MSnbase/reference/removeNoId-methods.md)
  documentation for more details and examples.

- removeMultipleAssignment:

  `signature(object = "MSnExp", fcol = "nprot")`: Removes protein groups
  (or feature belong to protein groups) with more than one member. The
  latter is defined by extracting a feature variable (default is
  `"nprot"`). Also removes non-identified features.

- idSummary:

  `signature(object = "MSnExp", ...)`: Prints a summary that lists the
  percentage of identified features per file (called `coverage`).

- show:

  `signature(object = "MSnExp")`: Displays object content as text.

- isolationWindow:

  `signature(object = "MSnExp", ...)`: Returns the isolation window
  offsets for the MS2 spectra. See `isolationWindow` in the `mzR`
  package for details.

- trimMz:

  `signature(object = "MSnExp")`: Trims the MZ range of all the spectra
  of the `MSnExp` instance. See
  [`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
  documentation for more details and examples.

- `isCentroided(object, k = 0.025, qtl = 0.9, verbose = TRUE)`:

  A heuristic assessing if the spectra in the `object` are in profile or
  centroided mode. The function takes the `qtl`th quantile top peaks,
  then calculates the difference between adjacent M/Z value and returns
  `TRUE` if the first quartile is greater than `k`. (See
  `MSnbase:::.isCentroided` for the code.) If `verbose` (default), a
  table indicating mode for all MS levels is printed.

  The function has been tuned to work for MS1 and MS2 spectra and data
  centroided using different peak picking algorithms, but false
  positives can occur. See
  <https://github.com/lgatto/MSnbase/issues/131> for details. For whole
  experiments, where all MS1 and MS2 spectra are expected to be in the
  same, albeit possibly different modes, it is advised to assign the
  majority result for MS1 and MS2 spectra, rather than results for
  individual spectra. See an example below.

- as:

  `signature(object = "MSnExp", "data.frame")`: Coerces the `MSnExp`
  object to a four-column `data.frame` with columns `"file"` (file index
  in `object`), `"rt"` (retention time), `"mz"` (m/z values) and `"i"`
  (intensity values).

- as:

  `signature(object = "MSnExp", "MSpectra")`: Coerces the `MSnExp`
  object to a
  [`MSpectra`](https://lgatto.github.io/MSnbase/reference/MSpectra.md)
  object with all feature annotations added as metadata columns
  (`mcols`).

Filtering and subsetting functions:

- filterRt:

  `signature(object = "MSnExp", rt = "numeric", msLevel. = "numeric")`:
  Retains MS spectra of level `msLevel.` with a retention times within
  `rt[1]` and `rt[2]`.

- filterMsLevel:

  `signature(object = "MSnExp", msLevel. = "numeric")`: Retains MS
  spectra of level `msLevel.`.

- filterPolarity:

  `signature(object = "MSnExp", polarity. = "numeric")`: Retains MS
  spectra of polarity `polarity.`.

- filterMz:

  `signature(object = "MSnExp", mz = "numeric", msLevel. = "numeric")`.
  See
  [`filterMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
  for details.

- filterFile:

  `signature(object = "MSnExp", file)`: Retains MS data of files
  matching the file index or file name provided with parameter `file`.

- filterAcquisitionNum:

- filterEmptySpectra:

  `signature(object = "MSnExp")`: Remove empty spectra from `object`
  (see `isEmpty`).

- filterPrecursorScan:

  `signature(object = "MSnExp", acquisitionNum = "numeric")`: Retain
  parent (e.g. MS1) and children scans (e.g. MS2) of `acquisitionNum`.
  See
  [`OnDiskMSnExp`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
  for an example.

- splitByFile:

  `signature(object = "MSnExp", f = "factor")`: split a `MSnExp` object
  by file into a `list` of `MSnExp` objects given the grouping in
  `factor` `f`.

- filterPrecursorMz:

  `signature(object = "MSnExp", mz, ppm = 10)`: retain spectra with a
  precursor m/z equal or similar to the one defined with parameter `mz`.
  Parameter `ppm` allows to define an accepted difference between the
  provided m/z and the spectrum's m/z.

- filterIsolationWindow:

  `signature(object = "MSnExp", mz)`: retain spectra with isolation
  windows that contain (which m/z range contain) the specified m/z.

## References

Information about the mzXML format as well converters from vendor
specific formats to mzXML:
<http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML>.

## Author

Laurent Gatto

## See also

`"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`
and
[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)
for loading `mzXML`, `mzData` or `mzML` files to generate an instance of
`MSnExp`.

The
`"`[`OnDiskMSnExp`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)`"`
manual page contains further details and examples.

[`chromatogram`](https://lgatto.github.io/MSnbase/reference/chromatogram-MSnExp-method.md)
to extract chromatographic data from a `MSnExp` or `OnDiskMSnExp`
object.

[`write`](https://rdrr.io/r/base/write.html) for the function to write
the data to mzML or mzXML file(s).

## Examples

``` r
mzxmlfile <- dir(system.file("extdata",package="MSnbase"),
                 pattern="mzXML",full.names=TRUE)
msnexp <- readMSData(mzxmlfile)
msnexp
#> MSn experiment data ("MSnExp")
#> Object size in memory: 0.18 Mb
#> - - - Spectra data - - -
#>  MS level(s): 2 
#>  Number of spectra: 5 
#>  MSn retention times: 25:01 - 25:02 minutes
#> - - - Processing information - - -
#> Data loaded: Fri Apr 10 14:42:31 2026 
#>  MSnbase version: 2.37.2 
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
```
