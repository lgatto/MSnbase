# The "Spectrum2" Class for MSn Spectra

`Spectrum2` extends the
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
class and introduces several MS2 specific attributes in addition to the
slots in
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`.
Since version 1.99.2, this class is used for any MS levels \> 1.
`Spectrum2` are not created directly but are contained in the
`assayData` slot of an
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`.

In version 1.19.12, the `polarity` slot had been added to the
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
class (previously in
`"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`).
Hence, `"Spectrum2"` objects created prior to this change will not be
valid anymore, since they will miss the `polarity` slots. Object can be
appropriately updated using the `updateObject` method.

## Slots

See the
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
class for inherited slots.

- `merged`::

  Object of class `"numeric"` indicating of how many combination the
  current spectrum is the result of.

- `precScanNum`::

  Object of class `"integer"` indicating the precursor MS scan index in
  the original input file. Accessed with the `precScanNum` or
  `precAcquisitionNum` methods.

- `precursorMz`::

  Object of class `"numeric"` providing the precursor ion MZ value.

- `precursorIntensity`::

  Object of class `"numeric"` providing the precursor ion intensity.

- `precursorCharge`::

  Object of class `"integer"` indicating the precursor ion charge.

- `collisionEnergy`::

  Object of class `"numeric"` indicating the collision energy used to
  fragment the parent ion.

## Extends

Class
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`,
directly. Class
`"`[`Versioned`](https://rdrr.io/pkg/Biobase/man/class.Versioned.html)`"`,
by class "Spectrum", distance 2.

## Methods

See
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
for additional accessors and methods for `Spectrum2` objects.

- `precursorMz(object)`:

  Returns the precursor MZ value as a numeric.

- `precursorMz(object)`:

  Returns the precursor scan number in the original data file as an
  integer.

- `precursorIntensity(object)`:

  Returns the precursor intensity as a numeric.

- `precursorCharge(object)`:

  Returns the precursor intensity as a integer.

- `collisionEnergy(object)`:

  Returns the collision energy as an numeric.

- `removeReporters(object, ...)`:

  Removes all reporter ion peaks. See
  [`removeReporters`](https://lgatto.github.io/MSnbase/reference/removeReporters-methods.md)
  documentation for more details and examples.

- `precAcquisitionNum`::

  Returns the precursor's acquisition number.

- `precScanNum`::

  See `precAcquisitionNum`.

- calculateFragments:

  `signature(sequence = "character", object = "Spectrum2")`: Calculates
  and matches the theoretical fragments of a peptide `sequence` with the
  ones observed in a spectrum. See
  [`calculateFragments`](https://lgatto.github.io/MSnbase/reference/calculateFragments-methods.md)
  documentation for more details and examples.

## Author

Laurent Gatto

## See also

Virtual super-class
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`,
`"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`
for MS1 spectra and
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
for a full experiment container.
