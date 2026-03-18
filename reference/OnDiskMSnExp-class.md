# The `OnDiskMSnExp` Class for MS Data And Meta-Data

Like the
[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
class, the `OnDiskMSnExp` class encapsulates data and meta-data for mass
spectrometry experiments, but does, in contrast to the former, not keep
the spectrum data in memory, but fetches the M/Z and intensity values on
demand from the raw files. This results in some instances to a reduced
performance, has however the advantage of a much smaller memory
footprint.

## Details

The `OnDiskMSnExp` object stores many spectrum related information into
the `featureData`, thus, some calls, like `rtime` to retrieve the
retention time of the individual scans does not require the raw data to
be read. Only M/Z and intensity values are loaded on-the-fly from the
original files. Extraction of values for individual scans is, for mzML
files, very fast. Extraction of the full data (all spectra) are
performed in a per-file parallel processing strategy.

Data manipulations related to spectras' M/Z or intensity values (e.g.
[`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
or
[`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md))
are (for `OnDiskMSnExp` objects) not applied immediately, but are stored
for later execution into the `spectraProcessingQueue`. The manipulations
are performed *on-the-fly* upon data retrieval. Other manipulations,
like removal of individual spectra are applied directly, since the
corresponding data is available in the object's `featureData` slot.

## Objects from the Class

Objects can be created by calls of the form `new("OnDiskMSnExp",...)`.
However, it is preferred to use the
[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)
function with argument `backend="disk"` that will read raw mass
spectrometry data to generate a valid `"OnDiskMSnExp"` instance.

## Slots

- `backend`::

  Character string specifying the used backend.

- `spectraProcessingQueue`::

  `list` of
  [`ProcessingStep`](https://lgatto.github.io/MSnbase/reference/ProcessingStep-class.md)
  objects defining the functions to be applied *on-the-fly* to the
  spectra data (M/Z and intensity duplets).

- `assayData`::

  Object of class `"environment"` that is however empty, as no spectrum
  data is stored. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `phenoData`::

  Object of class
  `"`[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`"`
  containing experimenter-supplied variables describing sample (i.e the
  individual tags for an labelled MS experiment) See
  [`phenoData`](https://rdrr.io/pkg/Biobase/man/phenoData.html) for more
  details. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `featureData`::

  Object of class
  `"`[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`"`
  containing variables describing features (spectra in our case). See
  [`featureData`](https://rdrr.io/pkg/Biobase/man/featureData.html) for
  more details. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `experimentData`::

  Object of class
  `"`[`MIAPE`](https://lgatto.github.io/MSnbase/reference/MIAPE-class.md)`"`,
  containing details of experimental methods. See
  [`experimentData`](https://rdrr.io/pkg/Biobase/man/abstract.html) for
  more details. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `protocolData`::

  Object of class
  `"`[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`"`
  containing equipment-generated variables (inherited from
  `"`[`eSet`](https://rdrr.io/pkg/Biobase/man/class.eSet.html)`"`). See
  [`protocolData`](https://rdrr.io/pkg/Biobase/man/protocolData.html)
  for more details. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `processingData`::

  Object of class
  `"`[`MSnProcess`](https://lgatto.github.io/MSnbase/reference/MSnProcess-class.md)`"`
  that records all processing. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.

- `.__classVersion__`::

  Object of class
  `"`[`Versions`](https://rdrr.io/pkg/Biobase/man/class.Versions.html)`"`
  describing the versions of R, the Biobase package,
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`
  and `MSnExp` of the current instance. Slot is inherited from
  `"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`.
  Intended for developer use and debugging (inherited from
  `"`[`eSet`](https://rdrr.io/pkg/Biobase/man/class.eSet.html)`"`).

## Extends

Class
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`,
directly. Class
`"`[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`"`,
by class "MSnExp", distance 3. Class
`"`[`VersionedBiobase`](https://rdrr.io/pkg/Biobase/man/class.VersionedBiobase.html)`"`,
by class "pSet", distance 4. Class
`"`[`Versioned`](https://rdrr.io/pkg/Biobase/man/class.Versioned.html)`"`,
by class "pSet", distance 5.

## Getter/setter methods

(in alphabetical order) See also methods for
[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
or [`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)
objects.

- \[:

  `object[i]`:subset the `OnDiskMSnExp` by spectra. `i` can be a
  `numeric` or `logical` vector specifying to which spectra the data set
  should be reduced (with `i` being the index of the spectrum in the
  object's `featureData`).

  The method returns a `OnDiskMSnExp` object with the data sub-set.

- \[\[:

  `object[[i]]`: extract s single spectrum from the `OnDiskMSnExp`
  object `object`. Argument `i` can be either numeric or character
  specifying the index or the name of the spectrum in the object (i.e.
  in the `featureData`). The relevant information will be extracted from
  the corresponding raw data file.

  The method returns a `Spectrum1` object.

- acquisitionNum:

  `acquisitionNum(signature(object="OnDiskMSnExp"))`: get the
  acquisition number of each spectrum in each individual file. The
  relevant information is extracted from the object's `featureData`
  slot.

  Returns a numeric vector with names corresponding to the spectrum
  names.

- assayData:

  `assayData(signature(object = "OnDiskMSnExp"))`: Extract the full
  data, i.e. read all spectra from the original files, apply all
  processing steps from the `spectraProcessingQueue` slot and return the
  data. Due to the required processing time accessing the full data
  should be avoided wherever possible.

  Returns an `environment`.

- centroided,centroided\<-:

  `centroided(signature(object="OnDiskMSnExp", msLevel, = "numeric"))`:
  whether individual spectra are centroided or uncentroided. The
  relevant information is extracted from the object's `featureData`
  slot. Returns a logical vector with names corresponding to the
  spectrum names. Use `centroided(object) <- value` to update the
  information, with value being a logical vector of length equal to the
  number of spectra in the experiment.

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
  individual spectra.

  See also
  [`isCentroidedFromFile`](https://lgatto.github.io/MSnbase/reference/isCentroidedFromFile.md)
  that accessed the mode directly from the raw data file.

- fromFile:

  `fromFile(signature(object = "OnDiskMSnExp"))`: get the index of the
  file (in `fileNames(object)`) from which the spectra were read. The
  relevant information is extracted from the object's `featureData`
  slot.

  Returns a numeric vector with names corresponding to the spectrum
  names.

- intensity:

  `intensity(signature(object="OnDiskMSnExp"))`: return the intensities
  from each spectrum in the data set. Intensities are first read from
  the raw files followed by an optional processing (depending on the
  processing steps defined in the `spectraProcessingQueue`). To reduce
  the amount of required memory, this is performed on a per-file basis.
  The `BPPARAM` argument allows to specify how and if parallel
  processing should be used. Information from individual files will be
  processed in parallel (one process per original file).

  The method returns a `list` of numeric intensity values. Each list
  element represents the intensities from one spectrum.

- ionCount:

  `ionCount(signature(object="OnDiskMSnExp", BPPARAM=bpparam()))`:
  extract the ion count (i.e. sum of intensity values) for each spectrum
  in the data set. The relevant data has to be extracted from the raw
  files (with eventually applying processing steps). The `BPPARAM`
  argument can be used to define how and if parallel processing should
  be used. Information from individual files will be processed in
  parallel (one process per original file).

  Returns a numeric vector with names corresponding to the spectrum
  names.

- isolationWindowLowerMz:

  `isolationWindowLowerMz(object = "OnDiskMSnExp")`: return the lower
  m/z boundary for the isolation window.

  Returns a numeric vector of length equal to the number of spectra with
  the lower m/z value of the isolation window or `NA` if not specified
  in the original file.

- isolationWindowUpperMz:

  `isolationWindowUpperMz(object = "OnDiskMSnExp")`: return the upper
  m/z boundary for the isolation window.

  Returns a numeric vector of length equal to the number of spectra with
  the upper m/z value of the isolation window or `NA` if not specified
  in the original file.

- length:

  `length(signature(object="OnDiskMSnExp"))`: Returns the number of
  spectra of the current experiment.

- msLevel:

  `msLevel(signature(object = "OnDiskMSnExp"))`: extract the MS level
  from the spectra. The relevant information is extracted from the
  object's `featureData` slot.

  Returns a numeric vector with names corresponding to the spectrum
  names.

- mz:

  `mz(signature(object="OnDiskMSnExp"))`: return the M/Z values from
  each spectrum in the data set. M/Z values are first read from the raw
  files followed by an optional processing (depending on the processing
  steps defined in the `spectraProcessingQueue`). To reduce the amount
  of required memory, this is performed on a per-file basis. The
  `BPPARAM` argument allows to specify how and if parallel processing
  should be used. Information from individual files will be processed in
  parallel (one process per original file).

  The method returns a `list` of numeric M/Z values. Each list element
  represents the values from one spectrum.

- peaksCount:

  `peaksCount(signature(object="OnDiskMSnExp", scans="numeric"), BPPARAM=bpparam())`:
  extrac the peaks count from each spectrum in the object. Depending on
  the eventually present `ProcessingStep` objects in the
  `spectraProcessingQueue` raw data will be loaded to calculate the
  peaks count. If no steps are present, the data is extracted from the
  `featureData`. Optional argument `scans` allows to specify the index
  of specific spectra from which the count should be returned. The
  `BPPARAM` argument can be used to define how and if parallel
  processing should be used. Information from individual files will be
  processed in parallel (one process per original file).

  Returns a numeric vector with names corresponding to the spectrum
  names.

- polarity:

  `polarity(signature(object="OnDiskMSnExp"))`: returns a numeric vector
  with the polarity of the individual spectra in the data set. The
  relevant information is extracted from the `featureData`.

- rtime:

  `rtime(signature(object="OnDiskMSnExp"))`: extrac the retention time
  of the individual spectra in the data set (from the `featureData`).

  Returns a numeric vector with names corresponding to the spectrum
  names.

- scanIndex:

  `scanIndex(signature(object="OnDiskMSnExp"))`: get the spectra scan
  indices within the respective file. The relevant information is
  extracted from the object's `featureData` slot. Returns a numeric
  vector of indices with names corresponding to the spectrum names.

- smoothed,smoothed\<-:

  `smoothed(signature(object="OnDiskMSnExp", msLevel. = "numeric"))`:
  whether individual spectra are smoothed or unsmoothed. The relevant
  information is extracted from the object's `featureData` slot. Returns
  a logical vector with names corresponding to the spectrum names. Use
  `smoothed(object) <- value` to update the information, with value
  being a logical vector of length equal to the number of spectra in the
  experiment.

- spectra:

  `spectra(signature(object="OnDiskMSnExp"), BPPARAM=bpparam())`:
  extract spectrum data from the individual files. This causes the
  spectrum data to be read from the original raw files. After that all
  processing steps defined in the `spectraProcessingQueue` are applied
  to it. The results are then returned as a `list` of
  [`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)
  objects.

  The `BPPARAM` argument can be used to define how and if parallel
  processing should be used. Information from individual files will be
  processed in parallel (one process per file). Note: extraction of
  selected spectra results in a considerable processing speed and should
  thus be preferred over whole data extraction.

  Returns a `list` of
  [`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)
  objects with names corresponding to the spectrum names.

- tic:

  `tic(signature(object="OnDiskMSnExp"), initial = TRUE, BPPARAM = bpparam())`:
  get the total ion current (TIC) of each spectrum in the data set. If
  `initial = TRUE`, the information is extracted from the object's
  `featureData` and represents the tic provided in the header of the
  original raw data files. For `initial = FALSE`, the TIC is calculated
  from the actual intensity values in each spectrum after applying all
  data manipulation methods (if any).

  See also https://github.com/lgatto/MSnbase/issues/332 for more
  details.

  `BPPARAM` parameter: see `spectra` method above.

  Returns a numeric vector with names corresponding to the spectrum
  names.

- bpi:

  `bpi(signature(object="OnDiskMSnExp"), initial = TRUE, BPPARAM = bpparam())`:
  get the base peak intensity (BPI), i.e. the maximum intensity from
  each spectrum in the data set. If `initial = TRUE`, the information is
  extracted from the object's `featureData` and represents the bpi
  provided in the header of the original raw data files. For
  `initial = FALSE`, the BPI is calculated from the actual intensity
  values in each spectrum after applying all eventual data manipulation
  methods.

  See also https://github.com/lgatto/MSnbase/issues/332 for more
  details.

  `BPPARAM` parameter: see `spectra` method above.

  Returns a numeric vector with names corresponding to the spectrum
  names.

- featureNames:

  `tic(signature(object="OnDiskMSnExp"))`: return a `character` of
  length `length(object)` containing the feature names. A replacement
  method is also available.

- spectrapply:

  `spectrapply(signature(object = "OnDiskMSnExp"), FUN = NULL, BPPARAM = bpparam(), ...)`:
  applies the function `FUN` to each spectrum passing additional
  parameters in `...` to that function and return its results. For
  `FUN = NULL` it returns the list of spectra (same as a call to
  `spectra`). Parameter `BPPARAM` allows to specify how and if parallel
  processing should be enabled.

  Returns a list with the result for each of spectrum.

## Data manipulation methods

(in alphabetical order) See also methods for
[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
or [`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)
objects. In contrast to the same-named methods for
[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md) or
[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
classes, the actual data manipulation is not performed immediately, but
only on-demand, e.g. when intensity or M/Z values are loaded.

- clean:

  `clean(signature(object="OnDiskMSnExp"), all=TRUE, verbose=TRUE)`: add
  an *clean* processing step to the lazy processing queue of the
  `OnDiskMSnExp` object. The `clean` command will only be executed when
  spectra information (including M/Z and intensity values) is requested
  from the `OnDiskMSnExp` object. Optional arguments to the methods are
  `all=TRUE` and `verbose=TRUE`.

  The method returns an `OnDiskMSnExp` object.

  For more details see documentation of the
  [`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)
  method.

- normalize:

  `normalize(signature(object="OnDiskMSnExp"), method=c("max", "sum"), ...)`:
  add a `normalize` processing step to the lazy processing queue of the
  returned `OnDiskMSnExp` object.

  The method returns an `OnDiskMSnExp` object.

  For more details see documentation of the
  [`normalize`](https://lgatto.github.io/MSnbase/reference/normalise-methods.md)
  method.

- removePeaks:

  `removePeaks(signature(object="OnDiskMSnExp"), t="min", verbose=TRUE)`:
  add a `removePeaks` processing step to the lazy processing queue of
  the returned `OnDiskMSnExp` object.

  The method returns an `OnDiskMSnExp` object.

  For more details see documentation of the
  [`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
  method.

- trimMz:

  `trimMz(signature(object="OnDiskMSnExp", mzlim="numeric"),...)`: add a
  `trimMz` processing step to the lazy processing queue of the returned
  `OnDiskMSnExp` object.

  The method returns an `OnDiskMSnExp` object.

  For more details see documentation of the
  [`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
  method.

## Other methods and functions

- validateOnDiskMSnExp:

  `validateOnDiskMSnExp(signature(object = "OnDiskMSnExp"))`: validates
  an `OnDiskMSnExp` object and all of its spectra. In addition to the
  *standard* `validObject` method, this method reads also all spectra
  from the original files, applies eventual processing steps and
  evaluates their validity.

- `as(from, "MSnExp")`:

  Converts the `OnDiskMSnExp` object `from`, to an in-memory `MSnExp`.
  Also available as an S3 method `as.MSnExp()`.

## Author

Johannes Rainer \<johannes.rainer@eurac.edu\>

## See also

[`pSet`](https://lgatto.github.io/MSnbase/reference/pSet-class.md),
[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md),
[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)

## Examples

``` r
## Get some example mzML files
library(msdata)
#> 
#>  IMPORTANT: msdata will be deprecated and replaced by the new
#>  MsDataHub package -- see https://bioconductor.org/packages/MsDataHub.
#>  Please open an issue at https://github.com/rformassspectrometry/MsDataHub/issues
#>  to inform us of any data from msdata that would need to be moved to MsDataHub.
mzfiles <- c(system.file("microtofq/MM14.mzML", package="msdata"),
       system.file("microtofq/MM8.mzML", package="msdata"))
## Read the data as an OnDiskMSnExp
odmse <- readMSData(mzfiles, msLevel=1, centroided = TRUE)

## Get the length of data, i.e. the total number of spectra.
length(odmse)
#> [1] 310

## Get the MS level
head(msLevel(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       1       1       1       1       1       1 

## Get the featureData, use fData to return as a data.frame
head(fData(odmse))
#>         spectrum
#> F1.S001        1
#> F1.S002        2
#> F1.S003        3
#> F1.S004        4
#> F1.S005        5
#> F1.S006        6

## Get to know from which file the spectra are
head(fromFile(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       1       1       1       1       1       1 

## And the file names:
fileNames(odmse)
#> [1] "/__w/_temp/Library/msdata/microtofq/MM14.mzML"
#> [2] "/__w/_temp/Library/msdata/microtofq/MM8.mzML" 

## Scan index and acquisitionNum
head(scanIndex(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       1       2       3       4       5       6 
head(acquisitionNum(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       1       2       3       4       5       6 

## Extract the spectra; the data is retrieved from the raw files.
head(spectra(odmse))
#> $F1.S001
#> Object of class "Spectrum1"
#>  Retention time: 4:30 
#>  MSn level: 1 
#>  Total ion count: 1378 
#>  Polarity: 1 
#> 
#> $F1.S002
#> Object of class "Spectrum1"
#>  Retention time: 4:31 
#>  MSn level: 1 
#>  Total ion count: 1356 
#>  Polarity: 1 
#> 
#> $F1.S003
#> Object of class "Spectrum1"
#>  Retention time: 4:31 
#>  MSn level: 1 
#>  Total ion count: 1404 
#>  Polarity: 1 
#> 
#> $F1.S004
#> Object of class "Spectrum1"
#>  Retention time: 4:31 
#>  MSn level: 1 
#>  Total ion count: 1496 
#>  Polarity: 1 
#> 
#> $F1.S005
#> Object of class "Spectrum1"
#>  Retention time: 4:32 
#>  MSn level: 1 
#>  Total ion count: 1525 
#>  Polarity: 1 
#> 
#> $F1.S006
#> Object of class "Spectrum1"
#>  Retention time: 4:32 
#>  MSn level: 1 
#>  Total ion count: 1498 
#>  Polarity: 1 
#> 

## Extracting individual spectra or a subset is much faster.
spectra(odmse[1:50])
#> $F1.S001
#> Object of class "Spectrum1"
#>  Retention time: 4:30 
#>  MSn level: 1 
#>  Total ion count: 1378 
#>  Polarity: 1 
#> 
#> $F1.S002
#> Object of class "Spectrum1"
#>  Retention time: 4:31 
#>  MSn level: 1 
#>  Total ion count: 1356 
#>  Polarity: 1 
#> 
#> $F1.S003
#> Object of class "Spectrum1"
#>  Retention time: 4:31 
#>  MSn level: 1 
#>  Total ion count: 1404 
#>  Polarity: 1 
#> 
#> $F1.S004
#> Object of class "Spectrum1"
#>  Retention time: 4:31 
#>  MSn level: 1 
#>  Total ion count: 1496 
#>  Polarity: 1 
#> 
#> $F1.S005
#> Object of class "Spectrum1"
#>  Retention time: 4:32 
#>  MSn level: 1 
#>  Total ion count: 1525 
#>  Polarity: 1 
#> 
#> $F1.S006
#> Object of class "Spectrum1"
#>  Retention time: 4:32 
#>  MSn level: 1 
#>  Total ion count: 1498 
#>  Polarity: 1 
#> 
#> $F1.S007
#> Object of class "Spectrum1"
#>  Retention time: 4:32 
#>  MSn level: 1 
#>  Total ion count: 1484 
#>  Polarity: 1 
#> 
#> $F1.S008
#> Object of class "Spectrum1"
#>  Retention time: 4:33 
#>  MSn level: 1 
#>  Total ion count: 1532 
#>  Polarity: 1 
#> 
#> $F1.S009
#> Object of class "Spectrum1"
#>  Retention time: 4:33 
#>  MSn level: 1 
#>  Total ion count: 1532 
#>  Polarity: 1 
#> 
#> $F1.S010
#> Object of class "Spectrum1"
#>  Retention time: 4:33 
#>  MSn level: 1 
#>  Total ion count: 1513 
#>  Polarity: 1 
#> 
#> $F1.S011
#> Object of class "Spectrum1"
#>  Retention time: 4:34 
#>  MSn level: 1 
#>  Total ion count: 1563 
#>  Polarity: 1 
#> 
#> $F1.S012
#> Object of class "Spectrum1"
#>  Retention time: 4:34 
#>  MSn level: 1 
#>  Total ion count: 1522 
#>  Polarity: 1 
#> 
#> $F1.S013
#> Object of class "Spectrum1"
#>  Retention time: 4:34 
#>  MSn level: 1 
#>  Total ion count: 1571 
#>  Polarity: 1 
#> 
#> $F1.S014
#> Object of class "Spectrum1"
#>  Retention time: 4:35 
#>  MSn level: 1 
#>  Total ion count: 1553 
#>  Polarity: 1 
#> 
#> $F1.S015
#> Object of class "Spectrum1"
#>  Retention time: 4:35 
#>  MSn level: 1 
#>  Total ion count: 1606 
#>  Polarity: 1 
#> 
#> $F1.S016
#> Object of class "Spectrum1"
#>  Retention time: 4:35 
#>  MSn level: 1 
#>  Total ion count: 1583 
#>  Polarity: 1 
#> 
#> $F1.S017
#> Object of class "Spectrum1"
#>  Retention time: 4:36 
#>  MSn level: 1 
#>  Total ion count: 1567 
#>  Polarity: 1 
#> 
#> $F1.S018
#> Object of class "Spectrum1"
#>  Retention time: 4:36 
#>  MSn level: 1 
#>  Total ion count: 1662 
#>  Polarity: 1 
#> 
#> $F1.S019
#> Object of class "Spectrum1"
#>  Retention time: 4:36 
#>  MSn level: 1 
#>  Total ion count: 1628 
#>  Polarity: 1 
#> 
#> $F1.S020
#> Object of class "Spectrum1"
#>  Retention time: 4:37 
#>  MSn level: 1 
#>  Total ion count: 1604 
#>  Polarity: 1 
#> 
#> $F1.S021
#> Object of class "Spectrum1"
#>  Retention time: 4:37 
#>  MSn level: 1 
#>  Total ion count: 1647 
#>  Polarity: 1 
#> 
#> $F1.S022
#> Object of class "Spectrum1"
#>  Retention time: 4:37 
#>  MSn level: 1 
#>  Total ion count: 1594 
#>  Polarity: 1 
#> 
#> $F1.S023
#> Object of class "Spectrum1"
#>  Retention time: 4:38 
#>  MSn level: 1 
#>  Total ion count: 1643 
#>  Polarity: 1 
#> 
#> $F1.S024
#> Object of class "Spectrum1"
#>  Retention time: 4:38 
#>  MSn level: 1 
#>  Total ion count: 1673 
#>  Polarity: 1 
#> 
#> $F1.S025
#> Object of class "Spectrum1"
#>  Retention time: 4:38 
#>  MSn level: 1 
#>  Total ion count: 1695 
#>  Polarity: 1 
#> 
#> $F1.S026
#> Object of class "Spectrum1"
#>  Retention time: 4:39 
#>  MSn level: 1 
#>  Total ion count: 1595 
#>  Polarity: 1 
#> 
#> $F1.S027
#> Object of class "Spectrum1"
#>  Retention time: 4:39 
#>  MSn level: 1 
#>  Total ion count: 1592 
#>  Polarity: 1 
#> 
#> $F1.S028
#> Object of class "Spectrum1"
#>  Retention time: 4:39 
#>  MSn level: 1 
#>  Total ion count: 1614 
#>  Polarity: 1 
#> 
#> $F1.S029
#> Object of class "Spectrum1"
#>  Retention time: 4:40 
#>  MSn level: 1 
#>  Total ion count: 1558 
#>  Polarity: 1 
#> 
#> $F1.S030
#> Object of class "Spectrum1"
#>  Retention time: 4:40 
#>  MSn level: 1 
#>  Total ion count: 1539 
#>  Polarity: 1 
#> 
#> $F1.S031
#> Object of class "Spectrum1"
#>  Retention time: 4:40 
#>  MSn level: 1 
#>  Total ion count: 1532 
#>  Polarity: 1 
#> 
#> $F1.S032
#> Object of class "Spectrum1"
#>  Retention time: 4:41 
#>  MSn level: 1 
#>  Total ion count: 1546 
#>  Polarity: 1 
#> 
#> $F1.S033
#> Object of class "Spectrum1"
#>  Retention time: 4:41 
#>  MSn level: 1 
#>  Total ion count: 1506 
#>  Polarity: 1 
#> 
#> $F1.S034
#> Object of class "Spectrum1"
#>  Retention time: 4:41 
#>  MSn level: 1 
#>  Total ion count: 1428 
#>  Polarity: 1 
#> 
#> $F1.S035
#> Object of class "Spectrum1"
#>  Retention time: 4:42 
#>  MSn level: 1 
#>  Total ion count: 1448 
#>  Polarity: 1 
#> 
#> $F1.S036
#> Object of class "Spectrum1"
#>  Retention time: 4:42 
#>  MSn level: 1 
#>  Total ion count: 1431 
#>  Polarity: 1 
#> 
#> $F1.S037
#> Object of class "Spectrum1"
#>  Retention time: 4:42 
#>  MSn level: 1 
#>  Total ion count: 1393 
#>  Polarity: 1 
#> 
#> $F1.S038
#> Object of class "Spectrum1"
#>  Retention time: 4:43 
#>  MSn level: 1 
#>  Total ion count: 1409 
#>  Polarity: 1 
#> 
#> $F1.S039
#> Object of class "Spectrum1"
#>  Retention time: 4:43 
#>  MSn level: 1 
#>  Total ion count: 1410 
#>  Polarity: 1 
#> 
#> $F1.S040
#> Object of class "Spectrum1"
#>  Retention time: 4:43 
#>  MSn level: 1 
#>  Total ion count: 1342 
#>  Polarity: 1 
#> 
#> $F1.S041
#> Object of class "Spectrum1"
#>  Retention time: 4:44 
#>  MSn level: 1 
#>  Total ion count: 1366 
#>  Polarity: 1 
#> 
#> $F1.S042
#> Object of class "Spectrum1"
#>  Retention time: 4:44 
#>  MSn level: 1 
#>  Total ion count: 1351 
#>  Polarity: 1 
#> 
#> $F1.S043
#> Object of class "Spectrum1"
#>  Retention time: 4:44 
#>  MSn level: 1 
#>  Total ion count: 1367 
#>  Polarity: 1 
#> 
#> $F1.S044
#> Object of class "Spectrum1"
#>  Retention time: 4:45 
#>  MSn level: 1 
#>  Total ion count: 1342 
#>  Polarity: 1 
#> 
#> $F1.S045
#> Object of class "Spectrum1"
#>  Retention time: 4:45 
#>  MSn level: 1 
#>  Total ion count: 1365 
#>  Polarity: 1 
#> 
#> $F1.S046
#> Object of class "Spectrum1"
#>  Retention time: 4:45 
#>  MSn level: 1 
#>  Total ion count: 1379 
#>  Polarity: 1 
#> 
#> $F1.S047
#> Object of class "Spectrum1"
#>  Retention time: 4:46 
#>  MSn level: 1 
#>  Total ion count: 1325 
#>  Polarity: 1 
#> 
#> $F1.S048
#> Object of class "Spectrum1"
#>  Retention time: 4:46 
#>  MSn level: 1 
#>  Total ion count: 1349 
#>  Polarity: 1 
#> 
#> $F1.S049
#> Object of class "Spectrum1"
#>  Retention time: 4:46 
#>  MSn level: 1 
#>  Total ion count: 1383 
#>  Polarity: 1 
#> 
#> $F1.S050
#> Object of class "Spectrum1"
#>  Retention time: 4:47 
#>  MSn level: 1 
#>  Total ion count: 1342 
#>  Polarity: 1 
#> 

## Alternatively, we could also subset the whole object by spectra and/or samples:
subs <- odmse[rtime(odmse) >= 2 & rtime(odmse) <= 20, ]
fileNames(subs)
#> [1] "/__w/_temp/Library/msdata/microtofq/MM8.mzML"
rtime(subs)
#>   F2.S006   F2.S007   F2.S008   F2.S009   F2.S010   F2.S011   F2.S012   F2.S013 
#>  2.169000  2.506002  2.842998  3.178998  3.516000  3.852000  4.188000  4.525002 
#>   F2.S014   F2.S015   F2.S016   F2.S017   F2.S018   F2.S019   F2.S020   F2.S021 
#>  4.861002  5.197998  5.533998  5.871000  6.207000  6.544020  6.880020  7.216980 
#>   F2.S022   F2.S023   F2.S024   F2.S025   F2.S026   F2.S027   F2.S028   F2.S029 
#>  7.554000  7.891020  8.227020  8.563980  8.899980  9.237000  9.573000  9.910020 
#>   F2.S030   F2.S031   F2.S032   F2.S033   F2.S034   F2.S035   F2.S036   F2.S037 
#> 10.246020 10.582020 10.918980 11.254980 11.592000 11.928000 12.265020 12.601020 
#>   F2.S038   F2.S039   F2.S040   F2.S041   F2.S042   F2.S043   F2.S044   F2.S045 
#> 12.937020 13.273980 13.609980 13.947000 14.283000 14.620020 14.956020 15.292020 
#>   F2.S046   F2.S047   F2.S048   F2.S049   F2.S050   F2.S051   F2.S052   F2.S053 
#> 15.628980 15.964980 16.302000 16.638000 16.975020 17.311020 17.649000 17.986020 
#>   F2.S054   F2.S055   F2.S056   F2.S057   F2.S058 
#> 18.322980 18.660000 18.996000 19.333020 19.669020 

## Extract intensities and M/Z values per spectrum; the methods return a list,
## each element representing the values for one spectrum.
ints <- intensity(odmse)
mzs <- mz(odmse)

## Return a data.frame with mz and intensity pairs for each spectrum from the
## object
res <- spectrapply(odmse, FUN = as, Class = "data.frame")

## Calling removePeaks, i.e. setting intensity values below a certain threshold to 0.
## Unlike the name suggests, this is not actually removing peaks. Such peaks with a 0
## intensity are then removed by the "clean" step.
## Also, the manipulations are not applied directly, but put into the "lazy"
## processing queue.
odmse <- removePeaks(odmse, t=10000)
odmse <- clean(odmse)

## The processing steps are only applied when actual raw data is extracted.
spectra(odmse[1:2])
#> $F1.S001
#> Object of class "Spectrum1"
#>  Retention time: 4:30 
#>  MSn level: 1 
#>  Total ion count: 0 
#>  Polarity: 1 
#> 
#> $F1.S002
#> Object of class "Spectrum1"
#>  Retention time: 4:31 
#>  MSn level: 1 
#>  Total ion count: 0 
#>  Polarity: 1 
#> 

## Get the polarity of the spectra.
head(polarity(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       1       1       1       1       1       1 

## Get the retention time of all spectra
head(rtime(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#> 270.334 270.671 271.007 271.343 271.680 272.016 

## Get the intensities after removePeaks and clean
intsAfter <- intensity(odmse)

head(lengths(ints))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>    1378    1356    1404    1496    1525    1498 
head(lengths(intsAfter))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       0       0       0       3       6       9 

## The same for the M/Z values
mzsAfter <- intensity(odmse)
head(lengths(mzs))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>    1378    1356    1404    1496    1525    1498 
head(lengths(mzsAfter))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       0       0       0       3       6       9 


## Centroided or profile mode
f <- msdata::proteomics(full.names = TRUE,
      pattern = "MS3TMT11.mzML")
odmse <- readMSData(f, mode = "onDisk")
validObject(odmse)
#> [1] TRUE
odmse[[1]]
#> Object of class "Spectrum1"
#>  Retention time: 45:27 
#>  MSn level: 1 
#>  Total ion count: 10768 
#>  Polarity: 1 

table(isCentroidedFromFile(odmse), msLevel(odmse))
#>        
#>           1   2   3
#>   FALSE  30   0   0
#>   TRUE    0 482 482

## centroided status could be set manually
centroided(odmse, msLevel = 1) <- FALSE
centroided(odmse, msLevel = 2) <- TRUE
centroided(odmse, msLevel = 3) <- TRUE

## or when reading the data
odmse2 <- readMSData(f, centroided = c(FALSE, TRUE, TRUE),
         mode = "onDisk")
table(centroided(odmse), msLevel(odmse))
#>        
#>           1   2   3
#>   FALSE  30   0   0
#>   TRUE    0 482 482

## Filtering precursor scans

head(acquisitionNum(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>   21945   21946   21947   21948   21949   21950 
head(msLevel(odmse))
#> F1.S001 F1.S002 F1.S003 F1.S004 F1.S005 F1.S006 
#>       1       2       2       3       2       2 

## Extract all spectra stemming from the first MS1 spectrum
(from1 <- filterPrecursorScan(odmse, 21945))
#> MSn experiment data ("OnDiskMSnExp")
#> Object size in memory: 0.05 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 2 3 
#>  Number of spectra: 35 
#>  MSn retention times: 45:27 - 45:30 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 18 16:56:19 2026] 
#> Filter: select parent/children scans for 21945 [Wed Mar 18 16:56:20 2026] 
#>  MSnbase version: 2.37.1 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: MS3TMT11.mzML
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   MS3TMT11.mzML 
#> protocolData: none
#> featureData
#>   featureNames: F1.S001 F1.S002 ... F1.S035 (35 total)
#>   fvarLabels: fileIdx spIdx ... spectrum (36 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
table(msLevel(from1))
#> 
#>  1  2  3 
#>  1 17 17 


## Extract the second sepctrum's parent (MS1) and children (MS3)
## spectra
(from2 <- filterPrecursorScan(odmse, 21946))
#> MSn experiment data ("OnDiskMSnExp")
#> Object size in memory: 0.03 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 2 3 
#>  Number of spectra: 3 
#>  MSn retention times: 45:27 - 45:27 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 18 16:56:19 2026] 
#> Filter: select parent/children scans for 21946 [Wed Mar 18 16:56:20 2026] 
#>  MSnbase version: 2.37.1 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: MS3TMT11.mzML
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   MS3TMT11.mzML 
#> protocolData: none
#> featureData
#>   featureNames: F1.S001 F1.S002 F1.S004
#>   fvarLabels: fileIdx spIdx ... spectrum (36 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
table(msLevel(from2))
#> 
#> 1 2 3 
#> 1 1 1 
```
