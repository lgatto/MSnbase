# Example `MSnExp` and `MSnSet` data sets

`itraqdata` is and example data sets is an iTRAQ 4-plex experiment that
has been run on an Orbitrap Velos instrument. It includes identification
data in the feature data slot obtain from the Mascot search engine. It
is a subset of an spike-in experiment where proteins have spiked in an
*Erwinia* background, as described in

Karp et al. (2010), *Addressing accuracy and precision issues in iTRAQ
quantitation*, Mol Cell Proteomics. 2010 Sep;9(9):1885-97. Epub 2010 Apr
10. (PMID 20382981).

The spiked-in proteins in `itradata` are BSA and ENO and are present in
relative abundances 1, 2.5, 5, 10 and 10, 5, 2.5, 1 in the 114, 115, 116
and 117 reporter tags.

The `msnset` object is produced by running the `quantify` method on the
`itraqdata` experimental data, as detailed in the
[`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)
example. This example data set is used in the MSnbase-demo vignette,
available with `vignette("MSnbase-demo", package="MSnbase")`.

The `msnset2` object is another example iTRAQ4 data that is used to
demonstrate features of the package, in particular the `iPQF` feature
aggregation method, described in
[`iPQF`](https://lgatto.github.io/MSnbase/reference/iPQF.md). It
corresponds to 11 proteins with spectra measurements from the original
data set described by Breitwieser et al. (2011) *General statistical
modeling of data from protein relative expression isobaric tags*. J.
Proteome Res., 10, 2758-2766.

## Usage

``` r
itraqdata
```

## Examples

``` r
data(itraqdata)
itraqdata
#> MSn experiment data ("MSnExp")
#> Object size in memory: 1.9 Mb
#> - - - Spectra data - - -
#>  MS level(s): 2 
#>  Number of spectra: 55 
#>  MSn retention times: 19:09 - 50:18 minutes
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> Updated from version 0.3.0 to 0.3.1 [Fri Jul  8 20:23:25 2016] 
#>  MSnbase version: 1.1.22 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: 1
#>   varLabels: sampleNames sampleNumbers
#>   varMetadata: labelDescription
#> Loaded from:
#>   dummyiTRAQ.mzXML 
#> protocolData: none
#> featureData
#>   featureNames: X1 X10 ... X9 (55 total)
#>   fvarLabels: spectrum ProteinAccession ProteinDescription
#>     PeptideSequence
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'

## created by
## msnset <- quantify(itraqdata, method = "trap", reporters = iTRAQ4)
data(msnset)
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
#> iTRAQ4 quantification by trapezoidation: Wed Apr  1 21:41:53 2015 
#>  MSnbase version: 1.1.22 

data(msnset2)
msnset2
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 614 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData: none
#> featureData
#>   featureNames: 11695 11696 ... 390 (614 total)
#>   fvarLabels: accession sequence ... rudist.org (9 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> Subset [13758,4][614,4] Mon Jun 15 21:03:40 2015 
#>  MSnbase version: 1.14.1 
```
