# Quantifies 'MSnExp' and 'Spectrum' objects

This method quantifies individual
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
objects or full
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
experiments. Current, MS2-level isobar tagging using iTRAQ and TMT (or
any arbitrary peaks of interest, see
`"`[`ReporterIons`](https://lgatto.github.io/MSnbase/reference/ReporterIons-class.md)`"`)
and MS2-level label-free quantitation (spectral counting, spectral index
or spectral abundance factor) are available.

Isobaric tag peaks of single spectra or complete experiments can be
quantified using appropriate `methods`. Label-free quantitation is
available only for `MSnExp` experiments.

Since version 1.13.5, parallel quantitation is supported by the
`BiocParallel` package and controlled by the `BPPARAM` argument.

## Arguments

- object:

  An instance of class
  `"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`
  (isobaric tagging only) or
  `"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`.

- method:

  Peak quantitation method. For isobaric tags, one of, possibly
  abreviated `"trapezoidation"`, `"max"`, or `"sum"`. These methods
  return respectively the area under the peak(s), the maximum of the
  peak(s) or the sum of all intensities of the peak(s).

  For label-free quantitation, one of `"SI"` (spectral index), `"SIgi"`
  (global intensity spectral index), `"SIn"` (normalised spectral
  index), `"SAF"` (spectral abundance factor) or `"NSAF"` (normalised
  spectral abundance factor).

  Finally, the simple `"count"` method counts the occurrence of the
  respective spectra (at this stage all 1s) that can then be used as
  input to
  [`combineFeatures`](https://lgatto.github.io/MSnbase/reference/combineFeatures.md)
  to implement spectra counting.

- reporters:

  An instance of class
  `"`[`ReporterIons`](https://lgatto.github.io/MSnbase/reference/ReporterIons-class.md)`"`
  that defines the peak(s) to be quantified. For isobaric tagging only.

- strict:

  For isobaric tagging only. If strict is `FALSE` (default), the
  quantitation is performed using data points along the entire width of
  a peak. If strict is set to `TRUE`, once the apex(es) is/are
  identified, only data points within apex +/- width of reporter (see
  `"`[`ReporterIons`](https://lgatto.github.io/MSnbase/reference/ReporterIons-class.md)`"`)
  are used for quantitation.

- BPPARAM:

  Support for parallel processing using the `BiocParallel`
  infrastructure. When missing (default), the default registered
  `BiocParallelParam` parameters are applied using
  [`bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
  Alternatively, one can pass a valid `BiocParallelParam` parameter
  instance: `SnowParam`, `MulticoreParam`, `DoparParam`, ... see the
  `BiocParallel` package for details.

- parallel:

  Deprecated. Please see `BPPARAM`.

- qual:

  Should the `qual` slot be populated. Default is `TRUE`.

- pepseq:

  A `character` giving the peptide sequence column in the feature data.
  Default is `"sequence"`.

- verbose:

  Verbose of the output (only for `MSnExp` objects).

- ...:

  Further arguments passed to the quantitation functions.

## Details

`"`[`ReporterIons`](https://lgatto.github.io/MSnbase/reference/ReporterIons-class.md)`"`
define specific MZ at which peaks are expected and a window around that
MZ value. A peak of interest is searched for in that window. Since
version 1.1.2, warnings are not thrown anymore in case no data is found
in that region or if the peak extends outside the window. This can be
checked manually after quantitation, by inspecting the quantitation data
(using the `exprs` accessor) for `NA` values or by comaring the
`lowerMz` and `upperMz` columns in the
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
`qual` slot against the respective expected `mz(reporters)` +/-
`width(reporters)`.

Once the range of the curve is found, quantification is performed. If no
data points are found in the expected region, `NA` is returned for the
reporter peak MZ.

Note that for label-free, spectra that have not been identified (the
corresponding fields in the feature data are populated with `NA` values)
or that have been uniquely assigned to a protein (the `nprot` feature
data is greater that 1) are removed prior to quantitation. The latter
does not apply for `method = "count"` but can be applied manually with
[`removeMultipleAssignment`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md).

## Methods

- `signature(object = "MSnExp", method = "character", reporters = "ReporterIons", verbose = "logical", ...)`:

  For isobaric tagging, quantifies peaks defined in `reporters` using
  `method` in all spectra of the `MSnExp` object. If verbose is set to
  `TRUE`, a progress bar will be displayed.

  For label-free quantitation, the respective quantitation methods and
  normalisations are applied to the spectra. These methods require two
  additional arguments (`...`), namely the protein accession of
  identifiers (`fcol`, with detault value `"DatabaseAccess"`) and the
  protein lengths (`plength`, with default value `"DBseqLength"`). These
  values are available of the identification data had been collated
  using
  [`addIdentificationData`](https://lgatto.github.io/MSnbase/reference/addIdentificationData-methods.md).

  An object of class
  `"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
  is returned containing the quantified feature expression and all meta
  data inherited from the `MSnExp` `object` argument.

- `signature(object = "Spectrum", method = "character", reporters = "ReporterIons")`:

  Quantifies peaks defined in `reporters` using `method` in the
  `Spectrum` object (isobaric tagging only).

  A list of length 2 will be returned. The first element, named
  `peakQuant`, is a 'numeric' of length equal to `length(reporters)`
  with quantitation of the reporter peaks using `method`.

  The second element, names `curveStats`, is a 'data.frame' of dimension
  `length(reporters)` times 7 giving, for each reporter curve
  parameters: maximum intensity (`'maxInt'`), number of maxima
  (`'nMaxInt'`), number of data points defined the curve
  (`'baseLength'`), lower and upper MZ values for the curve (`'lowerMz'`
  and `'upperMz'`), reporter (`'reporter'`) and precursor MZ value
  (`'precursor'`) when available.

## References

For details about the spectral index (SI), see Griffin NM, Yu J, Long F,
Oh P, Shore S, Li Y, Koziol JA, Schnitzer JE. *Label-free, normalized
quantification of complex mass spectrometry data for proteomic
analysis*. Nat Biotechnol. 2010 Jan;28(1):83-9. doi: 10.1038/nbt.1592.
PMID: 20010810; PubMed Central PMCID: PMC2805705.

For details about the spectra abundance factor, see Paoletti AC, Parmely
TJ, Tomomori-Sato C, Sato S, Zhu D, Conaway RC, Conaway JW, Florens L,
Washburn MP. *Quantitative proteomic analysis of distinct mammalian
Mediator complexes using normalized spectral abundance factors*. PNAS.
2006 Dec 12;103(50):18928-33. PMID: 17138671; PubMed Central PMCID:
PMC1672612.

## Author

Laurent Gatto and Sebastian Gibb

## Examples

``` r

## Quantifying a full experiment using iTRAQ4-plex tagging
data(itraqdata)
msnset <- quantify(itraqdata, method = "trap", reporters = iTRAQ4)
msnset
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
#> Updated from version 0.3.0 to 0.3.1 [Fri Jul  8 20:23:25 2016] 
#> iTRAQ4 quantification by trapezoidation: Fri Apr 10 15:50:44 2026 
#>  MSnbase version: 1.1.22 

## specifying a custom parallel framework
## bp <- MulticoreParam(2L) # on Linux/OSX
## bp <- SnowParam(2L) # on Windows
## quantify(itraqdata[1:10], method = "trap", iTRAQ4, BPPARAM = bp)

## Checking for non-quantified peaks
sum(is.na(exprs(msnset)))
#> [1] 1

## Quantifying a single spectrum
qty <- quantify(itraqdata[[1]], method = "trap", iTRAQ4[1])
qty$peakQuant
#> iTRAQ4.114 
#>   1347.616 
qty$curveStats
#>        maxInt nMaxInt baseLength  lowerMz  upperMz reporter precursor
#> [1,] 197455.8       1         11 114.1012 114.1196 114.1112  520.7833


## Label-free quantitation
## Raw (mzXML) and identification (mzid) files
quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")

msexp <- readMSData(quantFile)
msexp <- addIdentificationData(msexp, identFile)
fData(msexp)$DatabaseAccess
#> [1] "ECA0984" "ECA1028" NA        NA        "ECA0510"

si <- quantify(msexp, method = "SIn")
processingData(si)
#> - - - Processing information - - -
#> Data loaded: Fri Apr 10 15:50:45 2026 
#> Filtered 2 unidentified peptides out [Fri Apr 10 15:50:45 2026] 
#> Quantitation by total ion current [Fri Apr 10 15:50:45 2026] 
#> Combined 3 features into 3 using sum: Fri Apr 10 15:50:45 2026 
#> Quantification by SIn [Fri Apr 10 15:50:45 2026] 
#>  MSnbase version: 2.37.3 
exprs(si)
#>         dummyiTRAQ.mzXML
#> ECA0510     0.0006553518
#> ECA0984     0.0035384487
#> ECA1028     0.0002684726

saf <- quantify(msexp, method = "NSAF")
processingData(saf)
#> - - - Processing information - - -
#> Data loaded: Fri Apr 10 15:50:45 2026 
#> Filtered 2 unidentified peptides out [Fri Apr 10 15:50:46 2026] 
#> Filtered 0 unidentified peptides out [Fri Apr 10 15:50:46 2026] 
#> Quantitation by count [Fri Apr 10 15:50:46 2026] 
#> Combined 3 features into 3 using user-defined function: Fri Apr 10 15:50:46 2026 
#> Quantification by NSAF [Fri Apr 10 15:50:46 2026] 
#>  MSnbase version: 2.37.3 
exprs(saf)
#>         dummyiTRAQ.mzXML
#> ECA0510        0.4306167
#> ECA0984        0.3094475
#> ECA1028        0.2599359
```
