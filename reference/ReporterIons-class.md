# The "ReporterIons" Class

The `ReporterIons` class allows to define a set of isobaric reporter
ions that are used for quantification in MSMS mode, e.g. iTRAQ (isobaric
tag for relative and absolute quantitation) or TMT (tandem mass tags).
`ReporterIons` instances can them be used when quantifying
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
data of plotting the reporters peaks based on in
`"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`
ojects.

Some reporter ions are provided with `MSnbase` an can be loaded with the
[`data`](https://rdrr.io/r/utils/data.html) function. These reporter
ions data sets are:

- `iTRAQ4`::

  `ReporterIon` object for the iTRAQ 4-plex set. Load with
  `data(iTRAQ4)`.

- `iTRAQ5`::

  `ReporterIon` object for the iTRAQ 4-plex set plus the isobaric tag.
  Load with `data(iTRAQ5)`.

- `TMT6`::

  `ReporterIon` object for the TMT 6-plex set. Load with `data(TMT6)`.

- `TMT7`::

  `ReporterIon` object for the TMT 6-plex set plus the isobaric tag.
  Load with `data(TMT6)`.

## Objects from the Class

Objects can be created by calls of the form `new("ReporterIons", ...)`.

## Slots

- `name`::

  Object of class `"character"` to identify the `ReporterIons` instance.

- `reporterNames`::

  Object of class `"character"` naming each individual reporter of the
  `ReporterIons` instance. If not provided explicitely, they are names
  by concatenating the `ReporterIons` name and the respective MZ values.

- `description`::

  Object of class `"character"` to describe the `ReporterIons` instance.

- `mz`::

  Object of class `"numeric"` providing the MZ values of the reporter
  ions.

- `col`::

  Object of class `"character"` providing colours to highlight the
  reporters on plots.

- `width`::

  Object of class `"numeric"` indicating the width around the individual
  reporter ions MZ values were to search for peaks. This is dependent on
  the mass spectrometer's resolution and is used for peak picking when
  quantifying the reporters. See
  [`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)
  for more details about quantification.

- `.__classVersion__`::

  Object of class `"Versions"` indicating the version of the
  `ReporterIons` instance. Intended for developer use and debugging.

## Extends

Class
`"`[`Versioned`](https://rdrr.io/pkg/Biobase/man/class.Versioned.html)`"`,
directly.

## Methods

- `show(object)`:

  Displays object content as text.

- `object[]`:

  Subsets one or several reporter ions of the `ReporterIons` object and
  returns a new instance of the same class.

- `length(object)`:

  Returns the number of reporter ions in the instance.

- `mz(object, ...)`:

  Returns the expected mz values of reporter ions. Additional arguments
  are currently ignored.

- `reporterColours(object)` or reporterColors(object):

  Returns the colours used to highlight the reporter ions.

- `reporterNames(object)`:

  Returns the name of the individual reporter ions. If not specified or
  is an incorrect number of names is provided at initialisation, the
  names are generated automatically by concatenating the instance name
  and the reporter's MZ values.

- `reporterNames(object) <- value`:

  Sets the reporter names to `value`, which must be a character of the
  same length as the number of reporter ions.

- `width(object)`:

  Returns the widths in which the reporter ion peaks are expected.

- `names(object)`:

  Returns the name of the `ReporterIons` object.

- `description(object)`:

  Returns the description of the `ReporterIons` object.

## References

Ross PL, Huang YN, Marchese JN, Williamson B, Parker K, Hattan S,
Khainovski N, Pillai S, Dey S, Daniels S, Purkayastha S, Juhasz P,
Martin S, Bartlet-Jones M, He F, Jacobson A, Pappin DJ. "Multiplexed
protein quantitation in Saccharomyces cerevisiae using amine-reactive
isobaric tagging reagents." *Mol Cell Proteomics*, 2004
Dec;3(12):1154-69. Epub 2004 Sep 22. PubMed PMID: 15385600.

Thompson A, Schäfer J, Kuhn K, Kienle S, Schwarz J, Schmidt G, Neumann
T, Johnstone R, Mohammed AK, Hamon C. "Tandem mass tags: a novel
quantification strategy for comparative analysis of complex protein
mixtures by MS/MS." *Anal Chem.* 2003 Apr 15;75(8):1895-904. *Erratum*
in: *Anal Chem.* 2006 Jun 15;78(12):4235. Mohammed, A Karim A \[added\]
and *Anal Chem.* 2003 Sep 15;75(18):4942. Johnstone, R \[added\]. PubMed
PMID: 12713048.

## Author

Laurent Gatto

## See also

[`TMT6`](https://lgatto.github.io/MSnbase/reference/TMT6.md) or
[`iTRAQ4`](https://lgatto.github.io/MSnbase/reference/iTRAQ4.md) for
readily available examples.

## Examples

``` r
## Code used for the iTRAQ4 set
ri <- new("ReporterIons",
          description="4-plex iTRAQ",
          name="iTRAQ4",
          reporterNames=c("iTRAQ4.114","iTRAQ4.115",
                          "iTRAQ4.116","iTRAQ4.117"),
          mz=c(114.1,115.1,116.1,117.1),
          col=c("red","green","blue","yellow"),
          width=0.05)
ri
#> Object of class "ReporterIons"
#> iTRAQ4: '4-plex iTRAQ' with 4 reporter ions
#>  - [iTRAQ4.114] 114.1 +/- 0.05 (red)
#>  - [iTRAQ4.115] 115.1 +/- 0.05 (green)
#>  - [iTRAQ4.116] 116.1 +/- 0.05 (blue)
#>  - [iTRAQ4.117] 117.1 +/- 0.05 (yellow)
reporterNames(ri)
#> [1] "iTRAQ4.114" "iTRAQ4.115" "iTRAQ4.116" "iTRAQ4.117"
ri[1:2]
#> Object of class "ReporterIons"
#> iTRAQ4[1:2]: 'subset of 4-plex iTRAQ' with 2 reporter ions
#>  - [iTRAQ4.114] 114.1 +/- 0.05 (red)
#>  - [iTRAQ4.115] 115.1 +/- 0.05 (green)
```
