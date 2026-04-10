# Read an 'mzTab' file

This function can be used to create a
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
by reading and parsing an `mzTab` file. The metadata section is always
used to populate the `MSnSet`'s `experimentData` slot.

## Usage

``` r
readMzTabData_v0.9(file, what = c("PRT", "PEP"), verbose = isMSnbaseVerbose())
```

## Arguments

- file:

  A `character` with the `mzTab` file to be read in.

- what:

  One of `"PRT"` or `"PEP"`, defining which of protein of peptide
  section should be parse. The metadata section, when available, is
  always used to populate the `experimentData` slot.

- verbose:

  Produce verbose output.

## Value

An instance of class `MSnSet`.

## See also

[`writeMzTabData`](https://lgatto.github.io/MSnbase/reference/writeMzTabData.md)
to save an
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
as an `mzTab` file.

## Author

Laurent Gatto

## Examples

``` r
testfile <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/legacy/jmztab-1.0/examples/mztab_itraq_example.txt"

prot <- readMzTabData_v0.9(testfile, "PRT")
#> Warning: Version 0.9 is deprecated. Please see '?readMzTabData' and '?MzTab' for details.

prot
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 2 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: sub[1] sub[2] sub[3] sub[4]
#>   varLabels: abundance
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: 1 2
#>   fvarLabels: accession unit_id ... protein_abundance_std_error_sub[4]
#>     (26 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> mzTab read: Fri Apr 10 14:45:09 2026 
#>  MSnbase version: 2.37.2 

pep <- readMzTabData_v0.9(testfile, "PEP")
#> Warning: Version 0.9 is deprecated. Please see '?readMzTabData' and '?MzTab' for details.

pep
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 6 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: sub[1] sub[2] sub[3] sub[4]
#>   varLabels: abundance
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: 1 2 ... 6 (6 total)
#>   fvarLabels: sequence accession ... peptide_abundance_std_error_sub[4]
#>     (23 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#> mzTab read: Fri Apr 10 14:45:10 2026 
#>  MSnbase version: 2.37.2 
```
