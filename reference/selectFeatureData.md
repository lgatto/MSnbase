# Select feature variables of interest

Select feature variables to be retained.

`requiredFvarLabels` returns a `character` vector with the required
feature data variable names (`fvarLabels`, i.e. the column names in the
`fData` `data.frame`) for the specified object.

## Usage

``` r
selectFeatureData(object, graphics = TRUE, fcol)

requiredFvarLabels(x = c("OnDiskMSnExp", "MSnExp", "MSnSet"))
```

## Arguments

- object:

  An `MSnSet`, `MSnExp` or `OnDiskMSnExp`.

- graphics:

  A `logical` (default is `TRUE`) indicating whether a shiny application
  should be used if available. Otherwise, a text menu is used. Ignored
  if `k` is not missing.

- fcol:

  A `numeric`, `logical` or `character` of valid feature variables to be
  passed directly.

- x:

  `character(1)` specifying the class name for which the required
  feature data variable names should be returned.

## Value

For `selectFeatureData`: updated object containing only selected feature
variables.

For `requiredFvarLabels`: `character` with the required feature variable
names.

## Author

Laurent Gatto

## Examples

``` r

library("pRolocdata")
data(hyperLOPIT2015)
## 5 first feature variables
x <- selectFeatureData(hyperLOPIT2015, fcol = 1:5)
fvarLabels(x)
#> [1] "entry.name"          "protein.description" "peptides.rep1"      
#> [4] "peptides.rep2"       "psms.rep1"          
if (FALSE) { # \dontrun{
## select via GUI
x <- selectFeatureData(hyperLOPIT2015)
fvarLabels(x)
} # }

## Subset the feature data of an OnDiskMSnExp object to the minimal
## required columns
f <- system.file("microtofq/MM14.mzML", package = "msdata")
od <- readMSData(f, mode = "onDisk")

## what columns do we have?
fvarLabels(od)
#>  [1] "fileIdx"                    "spIdx"                     
#>  [3] "smoothed"                   "seqNum"                    
#>  [5] "acquisitionNum"             "msLevel"                   
#>  [7] "polarity"                   "originalPeaksCount"        
#>  [9] "totIonCurrent"              "retentionTime"             
#> [11] "basePeakMZ"                 "basePeakIntensity"         
#> [13] "collisionEnergy"            "electronBeamEnergy"        
#> [15] "ionisationEnergy"           "lowMZ"                     
#> [17] "highMZ"                     "precursorScanNum"          
#> [19] "precursorMZ"                "precursorCharge"           
#> [21] "precursorIntensity"         "mergedScan"                
#> [23] "mergedResultScanNum"        "mergedResultStartScanNum"  
#> [25] "mergedResultEndScanNum"     "injectionTime"             
#> [27] "filterString"               "spectrumId"                
#> [29] "centroided"                 "ionMobilityDriftTime"      
#> [31] "isolationWindowTargetMZ"    "isolationWindowLowerOffset"
#> [33] "isolationWindowUpperOffset" "scanWindowLowerLimit"      
#> [35] "scanWindowUpperLimit"       "spectrum"                  

## Reduce the feature data data.frame to the required columns only
od <- selectFeatureData(od, fcol = requiredFvarLabels(class(od)))
fvarLabels(od)
#> [1] "fileIdx"          "spIdx"            "acquisitionNum"   "retentionTime"   
#> [5] "msLevel"          "precursorScanNum"
```
