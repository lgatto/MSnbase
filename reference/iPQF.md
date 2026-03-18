# iPQF: iTRAQ (and TMT) Protein Quantification based on Features

The iPQF spectra-to-protein summarisation method integrates peptide
spectra characteristics and quantitative values for protein quantitation
estimation. Spectra features, such as charge state, sequence length,
identification score and others, contain valuable information concerning
quantification accuracy. The iPQF algorithm assigns weights to spectra
according to their overall feature reliability and computes a weighted
mean to estimate protein quantities. See also
[`combineFeatures`](https://lgatto.github.io/MSnbase/reference/combineFeatures.md)
for a more general overview of feature aggregation and examples.

## Usage

``` r
iPQF(
  object,
  groupBy,
  low.support.filter = FALSE,
  ratio.calc = "sum",
  method.combine = FALSE,
  feature.weight = c(7, 6, 4, 3, 2, 1, 5)^2
)
```

## Arguments

- object:

  An instance of class `MSnSet` containing absolute ion intensities.

- groupBy:

  Vector defining spectra to protein matching. Generally, this is a
  feature variable such as `fData(object)$accession`.

- low.support.filter:

  A `logical` specifying if proteins being supported by only 1-2 peptide
  spectra should be filtered out. Default is `FALSE`.

- ratio.calc:

  Either `"none"` (don't calculate any ratios), `"sum"` (default), or a
  specific channel (one of `sampleNames(object)`) defining how to
  calculate relative peptides intensities.

- method.combine:

  A `logical` defining whether to further use median polish to combine
  features.

- feature.weight:

  Vector `"numeric"` giving weight to the different features. Default is
  the squared order of the features redundant -unique-distance metric,
  charge state, ion intensity, sequence length, identification score,
  modification state, and mass based on a robustness analysis.

## Value

A `matrix` with estimated protein ratios.

## References

iPQF: a new peptide-to-protein summarization method using peptide
spectra characteristics to improve protein quantification. Fischer M,
Renard BY. Bioinformatics. 2016 Apr 1;32(7):1040-7.
doi:10.1093/bioinformatics/btv675. Epub 2015 Nov 20. PubMed
PMID:26589272.

## Author

Martina Fischer

## Examples

``` r
data(msnset2)
head(exprs(msnset2))
#>       X114.ions  X115.ions  X116.ions X117.ions
#> 11695 7.0384398  7.6349574  7.2886949 6.7348692
#> 11696 9.4098886 10.1034528 10.0697640 9.4273709
#> 11697 4.5859834  4.6926775  4.5016924 4.2924254
#> 11698 3.7262916  3.8913142  3.7928041 3.3721019
#> 11699 1.2229379  1.4005805  1.3940034 1.2957861
#> 11700 0.8303542  0.7397313  0.7686411 0.7402981
prot <- combineFeatures(msnset2,
                        groupBy = fData(msnset2)$accession,
                        method = "iPQF")
#> The following 1 proteins are only supported by 1 or 2 peptide spectra,
#> hence, protein quantification is not reliable and can only be calculated
#> by the 'mean' in these cases, corresponding protein accessions are:
#>   O95678
head(exprs(prot))
#>        X114.ions X115.ions X116.ions X117.ions
#> O95678 0.2404726 0.2682764 0.2584247 0.2328263
#> P01766 0.2610278 0.2467206 0.2544715 0.2377801
#> P01776 0.2678859 0.2591250 0.2423396 0.2306495
#> P02749 0.2640340 0.2523566 0.2510357 0.2325736
#> P02763 0.2503318 0.2524583 0.2501628 0.2470472
#> P07225 0.2533961 0.2506013 0.2504353 0.2455673
```
