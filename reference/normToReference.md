# Combine peptides into proteins.

This function combines peptides into their proteins by normalising the
intensity values to a reference run/sample for each protein.

## Usage

``` r
normToReference(
  x,
  group,
  reference = .referenceFractionValues(x = x, group = group)
)
```

## Arguments

- x:

  `matrix`, [`exprs`](https://rdrr.io/pkg/Biobase/man/exprs.html) matrix
  of an
  [MSnSet](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
  object.

- group:

  `double` or `factor`, grouping variable, i.e. protein accession; has
  to be of length equal `nrow(x)`.

- reference:

  `double`, vector of reference values, has to be of the same length as
  `group` and `nrow(x)`.

## Value

a matrix with one row per protein.

## Details

This function is not intented to be used directly (that's why it is not
exported via `NAMESPACE`). Instead the user should use
[`combineFeatures`](https://lgatto.github.io/MSnbase/reference/combineFeatures.md).

The algorithm is described in Nikolovski et al., briefly it works as
follows:

1.  Find reference run (column) for each protein (grouped rows). We use
    the run (column) with the lowest number of `NA`. If multiple
    candidates are available we use the one with the highest intensity.
    This step is skipped if the user use his own `reference` vector.

2.  For each protein (grouped rows) and each run (column):

    1.  Find peptides (grouped rows) shared by the current run (column)
        and the reference run (column).

    2.  Sum the shared peptides (grouped rows) for the current run
        (column) and the reference run (column).

    3.  The ratio of the shared peptides (grouped rows) of the current
        run (column) and the reference run (column) is the new intensity
        for the current protein for the current run.

## References

Nikolovski N, Shliaha PV, Gatto L, Dupree P, Lilley KS. Label-free
protein quantification for plant Golgi protein localization and
abundance. Plant Physiol. 2014 Oct;166(2):1033-43. DOI:
10.1104/pp.114.245589. PubMed PMID: 25122472.

## See also

[`combineFeatures`](https://lgatto.github.io/MSnbase/reference/combineFeatures.md)

## Author

Sebastian Gibb <mail@sebastiangibb.de>, Pavel Shliaha

## Examples

``` r
library("MSnbase")
data(msnset)

# choose the reference run automatically
combineFeatures(msnset, groupBy=fData(msnset)$ProteinAccession)
#> Your data contains missing values. Please read the relevant section in
#> the combineFeatures manual page for details on the effects of missing
#> values on data aggregation.
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 40 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#>   varLabels: mz reporters
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: BSA ECA0172 ... ENO (40 total)
#>   fvarLabels: spectrum ProteinAccession ... CV.iTRAQ4.117 (19 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> iTRAQ4 quantification by trapezoidation: Wed Apr  1 21:41:53 2015 
#> Combined 55 features into 40 using mean: Wed Mar 18 17:49:23 2026 
#>  MSnbase version: 2.37.2 

# use a user-given reference
combineFeatures(msnset, groupBy=fData(msnset)$ProteinAccession,
 reference=rep(2, 55))
#> Your data contains missing values. Please read the relevant section in
#> the combineFeatures manual page for details on the effects of missing
#> values on data aggregation.
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 40 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#>   varLabels: mz reporters
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: BSA ECA0172 ... ENO (40 total)
#>   fvarLabels: spectrum ProteinAccession ... CV.iTRAQ4.117 (19 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> iTRAQ4 quantification by trapezoidation: Wed Apr  1 21:41:53 2015 
#> Combined 55 features into 40 using mean: Wed Mar 18 17:49:23 2026 
#>  MSnbase version: 2.37.2 
```
