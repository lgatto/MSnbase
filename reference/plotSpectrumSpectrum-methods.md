# Plotting a 'Spectrum' vs another 'Spectrum' object.

These method plot mass spectra MZ values against the intensities as line
plots. The first spectrum is plotted in the upper panel and the other in
upside down in the lower panel. Common peaks are drawn in a slightly
darker colour. If a peptide sequence is provided it automatically
calculates and labels the fragments.

## Arguments

- x:

  Object of class
  `"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
  .

- y:

  Object of class
  `"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
  .

- ...:

  Further arguments passed to internal functions.

## Methods

- `signature(x = "Spectrum", y = "Spectrum", ...)`:

  Plots two spectra against each other. Common peaks are drawn in a
  slightly darker colour. The `...` arguments are passed to the internal
  functions. Currently `tolerance`, `relative`, `sequences` and most of
  the [`plot.default`](https://rdrr.io/r/graphics/plot.default.html)
  arguments (like `xlim`, `ylim`, `main`, `xlab`, `ylab`, ...) are
  supported. You could change the `tolerance` (default `25e-6`) and
  decide whether this tolerance should be applied relative (default
  `relative = TRUE`) or absolute (`relative = FALSE`) to find and colour
  common peaks. Use a `character` vector of length 2 to provide
  `sequences` which would be used to calculate and draw the
  corresponding fragments. If `sequences` are given the `type` argument
  (default: `type=c("b", "y")` specify the fragment types which should
  calculated. Also it is possible to allow some `modifications`.
  Therefore you have to apply a named `character` vector for
  `modifications` where the name corresponds to the one-letter-code of
  the modified amino acid (default: Carbamidomethyl
  `modifications=c(C=57.02146)`). Additional you can specifiy the type
  of `neutralLoss` (default:
  [`PSMatch::defaultNeutralLoss()`](https://rdrr.io/pkg/PSMatch/man/calculateFragments.html)).
  See
  [`calculateFragments`](https://lgatto.github.io/MSnbase/reference/calculateFragments-methods.md)
  for details.

  There are a lot of graphical arguments available to control the
  representation of the peaks and fragments. Use `peaks.pch` to set the
  character on top of the peaks (default: `peaks.pch=19`). In a similar
  way you can set the line width `peaks.lwd=1` and the magnification
  `peaks.cex=0.5` of the peaks. The size of the fragment/legend labels
  could be set using `fragments.cex=0.75` or `legend.cex` respectively.
  See [`par`](https://rdrr.io/r/graphics/par.html) for details about
  graphical parameters in general.

## Author

Sebastian Gibb \<mail@sebastiangibb.de\>

## See also

More spectrum plotting available in
[`plot.Spectrum`](https://lgatto.github.io/MSnbase/reference/plot-methods.md).

More details about fragment calculation:
[`calculateFragments`](https://lgatto.github.io/MSnbase/reference/calculateFragments-methods.md).

## Examples

``` r
## find path to a mzXML file
file <- dir(system.file(package = "MSnbase", dir = "extdata"),
            full.name = TRUE, pattern = "mzXML$")

## create basic MSnExp
msexp <- readMSData(file, centroided.=FALSE)

## centroid them
msexp <- pickPeaks(msexp)
#> Error: unable to find an inherited method for function ‘pickPeaks’ for signature ‘object = "MSnExp"’

## plot the first against the second spectrum
plot(msexp[[1]], msexp[[2]])
#> Your spectrum is not centroided.
#> Your spectrum is not centroided.


## add sequence information
plot(msexp[[1]], msexp[[2]], sequences=c("VESITARHGEVLQLRPK",
                                         "IDGQWVTHQWLKK"))
#> Your spectrum is not centroided.
#> Your spectrum is not centroided.



itraqdata2 <- pickPeaks(itraqdata)
#> Error: unable to find an inherited method for function ‘pickPeaks’ for signature ‘object = "MSnExp"’
(k <- which(fData(itraqdata2)[, "PeptideSequence"] == "TAGIQIVADDLTVTNPK"))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': error in evaluating the argument 'object' in selecting a method for function 'fData': object 'itraqdata2' not found
mzk <- precursorMz(itraqdata2)[k]
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'precursorMz': object 'itraqdata2' not found
zk <- precursorCharge(itraqdata2)[k]
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'precursorCharge': object 'itraqdata2' not found
mzk * zk
#> Error: object 'mzk' not found
plot(itraqdata2[[k[1]]], itraqdata2[[k[2]]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'plot': object 'itraqdata2' not found
```
