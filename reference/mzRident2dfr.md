# Coerce identification data to a `data.frame`

A function to convert the identification data contained in an `mzRident`
object to a `data.frame`. Each row represents a scan, which can however
be repeated several times if the PSM matches multiple proteins and/or
contains two or more modifications. To reduce the `data.frame` so that
rows/scans are unique and use semicolon-separated values to combine
information pertaining a scan, use
[`reduce`](https://lgatto.github.io/MSnbase/reference/reduce-data.frame-method.md).

## Arguments

- from:

  An object of class `mzRident` defined in the `mzR` package.

## Value

A `data.frame`

## Details

See also the *Tandem MS identification data* section in the
*MSnbase-demo* vignette.

## Author

Laurent Gatto

## Examples

``` r
## find path to a mzIdentML file
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")
library("mzR")
x <- openIDfile(identFile)
x
#> Identification file handle.
#> Filename:  dummyiTRAQ.mzid 
#> Number of psms:  4 
as(x, "data.frame")
#>            sequence spectrumID chargeState rank passThreshold
#> 1     IDGQWVTHQWLKK     scan=2           3    1          TRUE
#> 2 VESITARHGEVLQLRPK     scan=1           3    1          TRUE
#> 3 IKPQAVIETLHRLTEGK     scan=1           3    2          TRUE
#> 4           LVILLFR     scan=5           2    1          TRUE
#>   experimentalMassToCharge calculatedMassToCharge peptideRef modNum isDecoy
#> 1                 546.9586               546.9633       Pep1      0   FALSE
#> 2                 645.3741               645.0375       Pep2      0   FALSE
#> 3                 645.3741               645.0458       Pep3      0   FALSE
#> 4                 437.8040               437.2997       Pep4      0   FALSE
#>   post pre start end DatabaseAccess DBseqLength DatabaseSeq
#> 1    A   K    50  62        ECA1028         275            
#> 2    A   R   170 186        ECA0984         231            
#> 3    A   K   372 388        ECA3829         572            
#> 4    L   K    22  28        ECA0510         166            
#>                                                          DatabaseDescription
#> 1 ECA1028 2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-succinyltransferase
#> 2                                        ECA0984 DNA mismatch repair protein
#> 3                    ECA3829 acetolactate synthase isozyme III large subunit
#> 4           ECA0510 putative capsular polysacharide biosynthesis transferase
#>   scan.number.s. acquisitionNum     spectrumFile          idFile MS.GF.RawScore
#> 1              2              2 dummyiTRAQ.mzXML dummyiTRAQ.mzid            -30
#> 2              1              1 dummyiTRAQ.mzXML dummyiTRAQ.mzid            -39
#> 3              1              1 dummyiTRAQ.mzXML dummyiTRAQ.mzid            -39
#> 4              5              5 dummyiTRAQ.mzXML dummyiTRAQ.mzid            -42
#>   MS.GF.DeNovoScore MS.GF.SpecEValue MS.GF.EValue modPeptideRef modName modMass
#> 1                39     9.399048e-06     13.46615          <NA>    <NA>      NA
#> 2                77     5.527468e-05     79.36958          <NA>    <NA>      NA
#> 3                77     5.527468e-05     79.36958          <NA>    <NA>      NA
#> 4                 5     2.577830e-04    366.38422          <NA>    <NA>      NA
#>   modLocation subOriginalResidue subReplacementResidue subLocation
#> 1          NA               <NA>                  <NA>          NA
#> 2          NA               <NA>                  <NA>          NA
#> 3          NA               <NA>                  <NA>          NA
#> 4          NA               <NA>                  <NA>          NA
```
