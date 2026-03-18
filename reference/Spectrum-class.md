# The "Spectrum" Class

Virtual container for spectrum data common to all different types of
spectra. A `Spectrum` object can not be directly instanciated. Use
`"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`
and
`"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`
instead.

In version 1.19.12, the `polarity` slot has been added to this class
(previously in
`"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`).

## Slots

- `msLevel`::

  Object of class `"integer"` indicating the MS level: 1 for MS1 level
  `Spectrum1` objects and 2 for MSMSM `Spectrum2` objects. Levels \> 2
  have not been tested and will be handled as MS2 spectra.

- `polarity`::

  Object of class `"integer"` indicating the polarity if the ion.

- `peaksCount`::

  Object of class `"integer"` indicating the number of MZ peaks.

- `rt`::

  Object of class `"numeric"` indicating the retention time (in seconds)
  for the current ions.

- `tic`::

  Object of class `"numeric"` indicating the total ion current, as
  reported in the original raw data file.

- `acquisitionNum`::

  Object of class `"integer"` corresponding to the acquisition number of
  the current spectrum.

- `scanIndex`::

  Object of class `"integer"` indicating the scan index of the current
  spectrum.

- `mz`::

  Object of class `"numeric"` of length equal to the peaks count (see
  `peaksCount` slot) indicating the MZ values that have been measured
  for the current ion.

- `intensity`::

  Object of class `"numeric"` of same length as `mz` indicating the
  intensity at which each `mz` datum has been measured.

- `centroided`::

  Object of class `"logical"` indicating if instance is centroided
  ('TRUE') of uncentroided ('FALSE'). Default is `NA`.

- `smoothed`::

  Object of class `"logical"` indicating if instance is smoothed
  ('TRUE') of unsmoothed ('FALSE'). Default is `NA`.

- `fromFile`::

  Object of class `"integer"` referencing the file the spectrum
  originates. The file names are stored in the `processingData` slot of
  the
  `"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
  or
  `"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
  instance that contains the current `"Spectrum"` instance.

- `.__classVersion__`::

  Object of class `"Versions"` indicating the version of the `Spectrum`
  class. Intended for developer use and debugging.

## Extends

Class
`"`[`Versioned`](https://rdrr.io/pkg/Biobase/man/class.Versioned.html)`"`,
directly.

## Methods

- `acquisitionNum(object)`:

  Returns the acquisition number of the spectrum as an integer.

- `scanIndex(object)`:

  Returns the scan index of the spectrum as an integer.

- `centroided(object)`:

  Indicates whether spectrum is centroided (`TRUE`), in profile mode
  (`FALSE`), or unkown (`NA`).

- `isCentroided(object, k=0.025, qtl=0.9)`:

  A heuristic assessing if a spectrum is in profile or centroided mode.
  The function takes the `qtl`th quantile top peaks, then calculates the
  difference between adjacent M/Z value and returns `TRUE` if the first
  quartile is greater than `k`. (See `MSnbase:::.isCentroided` for the
  code.) The function has been tuned to work for MS1 and MS2 spectra and
  data centroided using different peak picking algorithms, but false
  positives can occur. See
  <https://github.com/lgatto/MSnbase/issues/131> for details. It should
  however be safe to use is at the experiment level, assuming that all
  MS level have the same mode. See
  [`class?MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
  for an example.

- `smoothed(object)`:

  Indicates whether spectrum is smoothed (`TRUE`) or not (`FALSE`).

- `centroided(object) <- value`:

  Sets the `centroided` status of the spectrum object.

- `smoothed(object) <- value`:

  Sets the `smoothed` status of the spectrum object.

- `fromFile(object)`:

  Returns the index of the raw data file from which the current
  instances originates as an integer.

- `intensity(object)`:

  Returns an object of class `numeric` containing the intensities of the
  spectrum.

- `msLevel(object)`:

  Returns an MS level of the spectrum as an integer.

- `mz(object, ...)`:

  Returns an object of class `numeric` containing the MZ value of the
  spectrum peaks. Additional arguments are currently ignored.

- `peaksCount(object)`:

  Returns the number of peaks (possibly of 0 intensity) as an integer.

- `rtime(object, ...)`:

  Returns the retention time for the spectrum as an integer. Additional
  arguments are currently ignored.

- `ionCount(object)`:

  Returns the total ion count for the spectrum as a numeric.

- `tic(object, ...)`:

  Returns the total ion current for the spectrum as a numeric.
  Additional arguments are currently ignored. This is the total ion
  current as originally reported in the raw data file. To get the
  current total ion count, use `ionCount`.

- bin:

  `signature(object = "Spectrum")`: Bins Spectrum. See
  [`bin`](https://lgatto.github.io/MSnbase/reference/bin-methods.md)
  documentation for more details and examples.

- clean:

  `signature(object = "Spectrum")`: Removes unused 0 intensity data
  points. See
  [`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)
  documentation for more details and examples.

- compareSpectra:

  `signature(x = "Spectrum", y = "Spectrum")`: Compares spectra. See
  [`compareSpectra`](https://lgatto.github.io/MSnbase/reference/compareSpectra-methods.md)
  documentation for more details and examples.

- estimateNoise:

  `signature(object = "Spectrum")`: Estimates the noise in a profile
  spectrum. See
  [`estimateNoise`](https://lgatto.github.io/MSnbase/reference/estimateNoise-method.md)
  documentation for more details and examples.

- pickPeaks:

  `signature(object = "Spectrum")`: Performs the peak picking to
  generate a centroided spectrum. See
  [`pickPeaks`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md)
  documentation for more details and examples.

- plot:

  `signature(x = "Spectrum", y = "missing")`: Plots intensity against
  mz. See
  [`plot.Spectrum`](https://lgatto.github.io/MSnbase/reference/plot-methods.md)
  documentation for more details.

- plot:

  `signature(x = "Spectrum", y = "Spectrum")`: Plots two spectra
  above/below each other. See
  [`plot.Spectrum.Spectrum`](https://lgatto.github.io/MSnbase/reference/plotSpectrumSpectrum-methods.md)
  documentation for more details.

- plot:

  `signature(x = "Spectrum", y = "character")`: Plots an MS2 level
  spectrum and its highlight the fragmention peaks. See
  [`plot.Spectrum.character`](https://lgatto.github.io/MSnbase/reference/plot-methods.md)
  documentation for more details.

- quantify:

  `signature(object = "Spectrum")`: Quatifies defined peaks in the
  spectrum. See
  [`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)
  documentation for more details.

- removePeaks:

  `signature(object = "Spectrum")`: Remove peaks lower that a threshold
  `t`. See
  [`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
  documentation for more details and examples.

- smooth:

  `signature(x = "Spectrum")`: Smooths spectrum. See
  [`smooth`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md)
  documentation for more details and examples.

- show:

  `signature(object = "Spectrum")`: Displays object content as text.

- trimMz:

  `signature(object = "Spectrum")`: Trims the MZ range of all the
  spectra of the `MSnExp` instance. See
  [`trimMz`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
  documentation for more details and examples.

- isEmpty:

  `signature(x = "Spectrum")`: Checks if the `x` is an empty `Spectrum`.

- as:

  `signature(object = "Spectrum", "data.frame")`: Coerces the `Spectrum`
  object to a two-column `data.frame` containing intensities and MZ
  values.

## Author

Laurent Gatto

## Note

This is a virtual class and can not be instanciated directly.

## See also

Instaciable sub-classes
`"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`
and
`"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`
for MS1 and MS2 spectra.
