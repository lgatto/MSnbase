# The "Spectrum1" Class for MS1 Spectra

`Spectrum1` extends the
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
class and introduces an MS1 specific attribute in addition to the slots
in
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`.
`Spectrum1` instances are not created directly but are contained in the
`assayData` slot of an
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`.

## Slots

See the
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
class for inherited slots.

## Extends

Class
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`,
directly. Class
`"`[`Versioned`](https://rdrr.io/pkg/Biobase/man/class.Versioned.html)`"`,
by class "Spectrum", distance 2.

## Methods

See
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
for additional accessors and methods to process `Spectrum1` objects.

- `polarity(object)`:

  Returns the polarity of the spectrum as an integer.

## Author

Laurent Gatto

## See also

Virtual super-class
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`,
`"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`
for MS2 spectra and
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
for a full experiment container.
