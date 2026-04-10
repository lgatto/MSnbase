# Normalisation of `MSnExp`, `MSnSet` and `Spectrum` objects

The `normalise` method (also available as `normalize`) performs basic
normalisation on spectra intensities of single spectra
(`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
or
`"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`
objects), whole experiments
(`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
objects) or quantified expression data
(`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
objects).

Raw spectra and experiments are normalised using `max` or `sum` only.
For MSMS spectra could be normalised to their `precursor` additionally.
Each peak intensity is divided by the highest intensity in the spectrum,
the sum of intensities or the intensity of the precursor. These methods
aim at facilitating relative peaks heights between different spectra.

The `method` parameter for
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
can be one of `sum`, `max`, `quantiles`, `center.mean`, `center.median`,
`.median`, `quantiles.robust` or `vsn`. For `sum` and `max`, each
feature's reporter intensity is divided by the maximum or the sum
respectively. These two methods are applied along the features (rows).

`center.mean` and `center.median` translate the respective sample
(column) intensities according to the column mean or median.
`diff.median` translates all samples (columns) so that they all match
the grand median. Using `quantiles` or `quantiles.robust` applies
(robust) quantile normalisation, as implemented in `normalize.quantiles`
and `normalize.quantiles.robust` of the `preprocessCore` package. `vsn`
uses the `vsn2` function from the `vsn` package. Note that the latter
also glog-transforms the intensities. See respective manuals for more
details and function arguments.

A `scale` method, mimicking the base `scale` method exists for
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
instances. See `?base::`[`scale`](https://rdrr.io/r/base/scale.html) for
details.

## Arguments

- object:

  An object of class
  `"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`,
  `"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`,
  `"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
  or
  `"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`.

- method:

  A character vector of length one that describes how to normalise the
  object. See description for details.

- ...:

  Additional arguments passed to the normalisation function.

## Methods

The `normalise` methods:

- `signature(object = "MSnSet", method = "character")`:

  Normalises the `object` reporter ions intensities using `method`.

- `signature(object = "MSnExp", method = "character")`:

  Normalises the `object` peak intensities using `method`.

- `signature(object = "Spectrum", method = "character")`:

  Normalises the `object` peak intensities using `method`.

- `signature(object = "Spectrum2", method = "character", precursorIntensity)`:

  Normalises the `object` peak intensities using `method`. If
  `method == "precursor"`, `precursorIntensity` allows to specify the
  intensity of the precursor manually.

The `scale` method:

- `signature(x = "MSnSet", center = "logical", scale = "logical")`:

  See `?base::`[`scale`](https://rdrr.io/r/base/scale.html).

## Examples

``` r
## quantifying full experiment
data(msnset)
msnset.nrm <- normalise(msnset, "quantiles")
msnset.nrm
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 55 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#>   varLabels: mz reporters
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: X1 X10 ... X9 (55 total)
#>   fvarLabels: spectrum ProteinAccession ... collision.energy (15 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation: No annotation 
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> iTRAQ4 quantification by trapezoidation: Wed Apr  1 21:41:53 2015 
#> Normalised (quantiles): Fri Apr 10 14:44:47 2026 
#>  MSnbase version: 1.1.22 
```
