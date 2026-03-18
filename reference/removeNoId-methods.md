# Removes non-identified features

The method removes non-identifed features in `MSnExp` and `MSnSet`
instances using relevant information from the `feaureData` slot of a
user-provide filtering vector of logicals.

## Methods

- `signature(object = "MSnExp", fcol = "pepseq", keep = NULL)`:

  Removes the feature from `object` that have a feature `fcol` (default
  is `"pepseq"`) equal to `NA`. Alternatively, one can also manually
  define `keep`, a vector of logical, defining the feature to be
  retained.

- `signature(object = "MSnSet", fcol = "pepseq", keep = NULL)`:

  As above of `MSnSet` instances.

## Author

Laurent Gatto

## See also

[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
and
[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md).

## Examples

``` r
  quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "mzXML$")
  identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
  msexp <- readMSData(quantFile)
  msexp <- addIdentificationData(msexp, identFile)
  fData(msexp)$sequence
#> [1] "VESITARHGEVLQLRPK" "IDGQWVTHQWLKK"     NA                 
#> [4] NA                  "LVILLFR"          
  length(msexp)
#> [1] 5

  ## using default fcol
  msexp2 <- removeNoId(msexp)
  length(msexp2)
#> [1] 3
  fData(msexp2)$sequence
#> [1] "VESITARHGEVLQLRPK" "IDGQWVTHQWLKK"     "LVILLFR"          

  ## using keep
  print(fvarLabels(msexp))
#>  [1] "spectrum"                 "acquisition.number"      
#>  [3] "sequence"                 "chargeState"             
#>  [5] "rank"                     "passThreshold"           
#>  [7] "experimentalMassToCharge" "calculatedMassToCharge"  
#>  [9] "peptideRef"               "modNum"                  
#> [11] "isDecoy"                  "post"                    
#> [13] "pre"                      "start"                   
#> [15] "end"                      "DatabaseAccess"          
#> [17] "DBseqLength"              "DatabaseSeq"             
#> [19] "DatabaseDescription"      "scan.number.s."          
#> [21] "idFile"                   "MS.GF.RawScore"          
#> [23] "MS.GF.DeNovoScore"        "MS.GF.SpecEValue"        
#> [25] "MS.GF.EValue"             "modPeptideRef"           
#> [27] "modName"                  "modMass"                 
#> [29] "modLocation"              "subOriginalResidue"      
#> [31] "subReplacementResidue"    "subLocation"             
#> [33] "nprot"                    "npep.prot"               
#> [35] "npsm.prot"                "npsm.pep"                
  (k <- fData(msexp)$'MS.GF.EValue' > 75)
#> [1]  TRUE FALSE    NA    NA  TRUE
  k[is.na(k)] <- FALSE
  k
#> [1]  TRUE FALSE FALSE FALSE  TRUE
  msexp3 <- removeNoId(msexp, keep = k)
  length(msexp3)
#> [1] 2
  fData(msexp3)$sequence
#> [1] "VESITARHGEVLQLRPK" "LVILLFR"          
```
