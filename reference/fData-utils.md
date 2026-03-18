# Expand or merge feature variables

The `expandFeatureVars` and `mergeFeatureVars` respectively expand and
merge groups of feature variables. Using these functions, a set of
columns in a feature data can be merged into a single new
data.frame-column variables and a data.frame-column can be expanded into
single feature columns. The original feature variables are removed.

## Usage

``` r
expandFeatureVars(x, fcol, prefix)

mergeFeatureVars(x, fcol, fcol2)
```

## Arguments

- x:

  An object of class `MSnSet`.

- fcol:

  A [`character()`](https://rdrr.io/r/base/character.html) of feature
  variables to expand (for `expandFeatureVars`) or merge (for
  `mergeFeatureVars`).

- prefix:

  A `character(1)` to use as prefix to the new feature variables. If
  missing (default), then `fcol` is used instead. If `NULL`, then no
  prefix is used.

- fcol2:

  A `character(1)` defining the name of the new feature variable.

## Value

An `MSnSet` for expanded (merged) feature variables.

## Author

Laurent Gatto

## Examples

``` r
library("pRolocdata")
data(hyperLOPIT2015)
fvarLabels(hyperLOPIT2015)
#>  [1] "entry.name"                "protein.description"      
#>  [3] "peptides.rep1"             "peptides.rep2"            
#>  [5] "psms.rep1"                 "psms.rep2"                
#>  [7] "phenodisco.input"          "phenodisco.output"        
#>  [9] "curated.phenodisco.output" "markers"                  
#> [11] "svm.classification"        "svm.score"                
#> [13] "svm.top.quartile"          "final.assignment"         
#> [15] "first.evidence"            "curated.organelles"       
#> [17] "cytoskeletal.components"   "trafficking.proteins"     
#> [19] "protein.complexes"         "signalling.cascades"      
#> [21] "oct4.interactome"          "nanog.interactome"        
#> [23] "sox2.interactome"          "cell.surface.proteins"    
#> [25] "markers2015"               "TAGM"                     
## Let's merge all svm prediction feature variables
(k <- grep("^svm", fvarLabels(hyperLOPIT2015), value = TRUE))
#> [1] "svm.classification" "svm.score"          "svm.top.quartile"  
hl <- mergeFeatureVars(hyperLOPIT2015, fcol = k, fcol2 = "SVM")
fvarLabels(hl)
#>  [1] "entry.name"                "protein.description"      
#>  [3] "peptides.rep1"             "peptides.rep2"            
#>  [5] "psms.rep1"                 "psms.rep2"                
#>  [7] "phenodisco.input"          "phenodisco.output"        
#>  [9] "curated.phenodisco.output" "markers"                  
#> [11] "final.assignment"          "first.evidence"           
#> [13] "curated.organelles"        "cytoskeletal.components"  
#> [15] "trafficking.proteins"      "protein.complexes"        
#> [17] "signalling.cascades"       "oct4.interactome"         
#> [19] "nanog.interactome"         "sox2.interactome"         
#> [21] "cell.surface.proteins"     "markers2015"              
#> [23] "TAGM"                      "SVM"                      
head(fData(hl)$SVM)
#>                             svm.classification svm.score   svm.top.quartile
#> Q9JHU4                              Peroxisome     0.303            unknown
#> Q9QXS1-3                          60S Ribosome     0.223            unknown
#> Q9ERU9                     Nucleus - Chromatin     0.737            unknown
#> P26039                      Actin cytoskeleton     1.000 Actin cytoskeleton
#> Q8BTM8   Endoplasmic reticulum/Golgi apparatus     0.218            unknown
#> A2ARV4                                Lysosome     0.681            unknown

## Let's expand the new SVM into individual columns
hl2 <- expandFeatureVars(hl, "SVM")
fvarLabels(hl2)
#>  [1] "entry.name"                "protein.description"      
#>  [3] "peptides.rep1"             "peptides.rep2"            
#>  [5] "psms.rep1"                 "psms.rep2"                
#>  [7] "phenodisco.input"          "phenodisco.output"        
#>  [9] "curated.phenodisco.output" "markers"                  
#> [11] "final.assignment"          "first.evidence"           
#> [13] "curated.organelles"        "cytoskeletal.components"  
#> [15] "trafficking.proteins"      "protein.complexes"        
#> [17] "signalling.cascades"       "oct4.interactome"         
#> [19] "nanog.interactome"         "sox2.interactome"         
#> [21] "cell.surface.proteins"     "markers2015"              
#> [23] "TAGM"                      "SVM.svm.classification"   
#> [25] "SVM.svm.score"             "SVM.svm.top.quartile"     
## We can set the prefix manually
hl2 <- expandFeatureVars(hl, "SVM", prefix = "Expanded")
fvarLabels(hl2)
#>  [1] "entry.name"                  "protein.description"        
#>  [3] "peptides.rep1"               "peptides.rep2"              
#>  [5] "psms.rep1"                   "psms.rep2"                  
#>  [7] "phenodisco.input"            "phenodisco.output"          
#>  [9] "curated.phenodisco.output"   "markers"                    
#> [11] "final.assignment"            "first.evidence"             
#> [13] "curated.organelles"          "cytoskeletal.components"    
#> [15] "trafficking.proteins"        "protein.complexes"          
#> [17] "signalling.cascades"         "oct4.interactome"           
#> [19] "nanog.interactome"           "sox2.interactome"           
#> [21] "cell.surface.proteins"       "markers2015"                
#> [23] "TAGM"                        "Expanded.svm.classification"
#> [25] "Expanded.svm.score"          "Expanded.svm.top.quartile"  
## If we don't want any prefix
hl2 <- expandFeatureVars(hl, "SVM", prefix = NULL)
fvarLabels(hl2)
#>  [1] "entry.name"                "protein.description"      
#>  [3] "peptides.rep1"             "peptides.rep2"            
#>  [5] "psms.rep1"                 "psms.rep2"                
#>  [7] "phenodisco.input"          "phenodisco.output"        
#>  [9] "curated.phenodisco.output" "markers"                  
#> [11] "final.assignment"          "first.evidence"           
#> [13] "curated.organelles"        "cytoskeletal.components"  
#> [15] "trafficking.proteins"      "protein.complexes"        
#> [17] "signalling.cascades"       "oct4.interactome"         
#> [19] "nanog.interactome"         "sox2.interactome"         
#> [21] "cell.surface.proteins"     "markers2015"              
#> [23] "TAGM"                      "svm.classification"       
#> [25] "svm.score"                 "svm.top.quartile"         
```
