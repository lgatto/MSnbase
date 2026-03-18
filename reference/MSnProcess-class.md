# The "MSnProcess" Class

`MSnProcess` is a container for MSnExp and MSnSet processing
information. It records data files, processing steps, thresholds,
analysis methods and times that have been applied to MSnExp or MSnSet
instances.

## Slots

- `files`::

  Object of class `"character"` storing the raw data files used in
  experiment described by the `"MSnProcess"` instance.

- `processing`::

  Object of class `"character"` storing all the processing steps and
  times.

- `merged`::

  Object of class `"logical"` indicating whether spectra have been
  merged.

- `cleaned`::

  Object of class `"logical"` indicating whether spectra have been
  cleaned. See
  [`clean`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)
  for more details and examples.

- `removedPeaks`::

  Object of class `"character"` describing whether peaks have been
  removed and which threshold was used. See
  [`removePeaks`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
  for more details and examples.

- `smoothed`::

  Object of class `"logical"` indicating whether spectra have been
  smoothed.

- `trimmed`::

  Object of class `"numeric"` documenting if/how the data has been
  trimmed.

- `normalised`::

  Object of class `"logical"` describing whether and how data have been
  normalised.

- `MSnbaseVersion`::

  Object of class `"character"` indicating the version of MSnbase.

- `.__classVersion__`::

  Object of class `"Versions"` indicating the version of the
  `MSnProcess` instance. Intended for developer use and debugging.

## Extends

Class `"Versioned"`, directly.

## Methods

- fileNames:

  `signature(object = "MSnProcess")`: Returns the file names used in
  experiment described by the `"MSnProcess"` instance.

- show:

  `signature(object = "MSnProcess")`: Displays object content as text.

- combine:

  `signature(x = "MSnProcess", y = "MSnProcess")`: Combines multiple
  `MSnProcess` instances.

## Author

Laurent Gatto

## Note

This class is likely to be updated using an `AnnotatedDataFrame`.

## See also

See the
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
and
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
classes that actually use `MSnProcess` as a slot.

## Examples

``` r
showClass("MSnProcess")
#> Class "MSnProcess" [package "MSnbase"]
#> 
#> Slots:
#>                                                                               
#> Name:              files        processing            merged           cleaned
#> Class:         character         character           logical           logical
#>                                                                               
#> Name:       removedPeaks          smoothed           trimmed        normalised
#> Class:         character           logical           numeric           logical
#>                                           
#> Name:     MSnbaseVersion .__classVersion__
#> Class:         character          Versions
#> 
#> Extends: "Versioned"
```
