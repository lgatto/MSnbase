# Write an experiment or spectrum to an mgf file

Methods `writeMgfData` write individual
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
instances of whole
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
experiments to a file in Mascot Generic Format (mgf) (see
[http://www.matrixscience.com/help/data_file_help.html](http://www.matrixscience.com/help/data_file_help.md)
for more details). Function `readMgfData` read spectra from and mgf file
and creates an
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
object.

## Arguments

- object:

  An instance of class
  `"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
  or
  `"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`.

- con:

  A valid `connection` or a `character` string with the name of the file
  to save the object. In case of the latter, a `file` connection is
  created. If not specified, 'spectrum.mgf' or 'experiment.mgf' are used
  depending on the class of `object`. Note that existing files are
  overwritted.

- COM:

  Optional character vector with the value for the 'COM' field.

- TITLE:

  Optional character vector with the value for the spectrum 'TITLE'
  field. Not applicable for experiments.

## Methods

- `signature(object = "MSnExp")`:

  Writes the full exeriment to an mgf file.

- `signature(object = "Spectrum")`:

  Writes an individual spectrum to an mgf file.

## See also

[`readMgfData`](https://lgatto.github.io/MSnbase/reference/readMgfData.md)
function to read data from and mgf file.

## Details

Note that when reading an mgf file, the original order of the spectra is
lost. Thus, if the data was originally written to mgf from an `MSnExp`
object using `writeMgfData`, although the feature names will be
identical, the spectra are not as a result of the reordering. See
example below.

## Examples

``` r
data(itraqdata)

f <- tempfile()

writeMgfData(itraqdata, con = f)

itraqdata2 <- readMgfData(f)

## note that the order of the spectra and precision of some values
## (precursorMz for instance) are altered
match(signif(precursorMz(itraqdata2),4),
      signif(precursorMz(itraqdata),4))
#>  [1]  1 10 11 12 13 14 15 16 17 18 19  2 20 21 22 23 24 25 26 27 28 29  3 30 31
#> [26] 32 33 34 35 36 37 38 39  4 40 41 42 43 44 45 46 47 48 49  5 50 51 52 53 54
#> [51] 55  6  7  8  9

## [1]  1 10 11 12 13 14 15 16 17 18 ...
## ... but all the precursors are there
all.equal(sort(precursorMz(itraqdata2)),
          sort(precursorMz(itraqdata)),
          check.attributes = FALSE,
          tolerance = 10e-5)
#> [1] TRUE

all.equal(as.data.frame(itraqdata2[[1]]),
          as.data.frame(itraqdata[[1]]))
#> [1] TRUE

all.equal(as.data.frame(itraqdata2[[3]]),
          as.data.frame(itraqdata[[11]]))
#> [1] TRUE

all(featureNames(itraqdata2) == featureNames(itraqdata))
#> [1] TRUE
```
