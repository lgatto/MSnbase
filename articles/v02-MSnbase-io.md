# MSnbase IO capabilities

Abstract

This vignette describes *MSnbase*’s input and output capabilities.

## Foreword

This software is free and open-source software. If you use it, please
support the project by citing it in publications:

> Gatto L, Lilley KS. MSnbase-an R/Bioconductor package for isobaric
> tagged mass spectrometry data visualization, processing and
> quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. doi:
> [10.1093/bioinformatics/btr645](https://doi.org/10.1093/bioinformatics/btr645).
> PMID: [22113085](https://www.ncbi.nlm.nih.gov/pubmed/22113085).

> *`MSnbase`, efficient and elegant R-based processing and visualisation
> of raw mass spectrometry data*. Laurent Gatto, Sebastian Gibb,
> Johannes Rainer. bioRxiv 2020.04.29.067868; doi:
> <https://doi.org/10.1101/2020.04.29.067868>

## Questions and bugs

For bugs, typos, suggestions or other questions, please file an issue in
our tracking system (<https://github.com/lgatto/MSnbase/issues>)
providing as much information as possible, a reproducible example and
the output of
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

If you don’t have a GitHub account or wish to reach a broader audience
for general questions about proteomics analysis using R, you may want to
use the Bioconductor support site: <https://support.bioconductor.org/>.

## Overview

*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)*’s aims are
to facilitate the reproducible analysis of mass spectrometry data within
the R environment, from raw data import and processing, feature
quantification, quantification and statistical analysis of the results
(Gatto and Lilley 2012). Data import functions for several formats are
provided and intermediate or final results can also be saved or
exported. These capabilities are presented below.

## Data input

##### Raw data

Data stored in one of the published `XML`-based formats. i.e. `mzXML`
(Pedrioli et al. 2004), `mzData` (Orchard et al. 2007) or `mzML`
(Martens et al. 2010), can be imported with the `readMSData` method,
which makes use of the
*[mzR](https://bioconductor.org/packages/3.23/mzR)* package to create
`MSnExp` objects. The files can be in profile or centroided mode. See
[`?readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)
for details.

Data from `mzML` files containing chromatographic data (e.g. generated
in SRM/MRM experiments) can be imported with the `readSRMData` function
that returns the chromatographic data as a `MChromatograms` object. See
[`?readSRMData`](https://lgatto.github.io/MSnbase/reference/readSRMData.md)
for more details.

##### Peak lists

Peak lists in the `mgf` format[^1] can be imported using the
`readMgfData`. In this case, the peak data has generally been
pre-processed by other software. See
[`?readMgfData`](https://lgatto.github.io/MSnbase/reference/readMgfData.md)
for details.

##### Quantitation data

Third party software can be used to generate quantitative data and
exported as a spreadsheet (generally comma or tab separated format).
This data as well as any additional meta-data can be imported with the
`readMSnSet` function. See
[`?readMSnSet`](https://lgatto.github.io/MSnbase/reference/readMSnSet.md)
for details.

*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* also
supports the `mzTab` format[^2], a light-weight, tab-delimited file
format for proteomics data developed within the Proteomics Standards
Initiative (PSI). `mzTab` files can be read into R with `readMzTabData`
to create and `MSnSet` instance.

![MSnbase input capabilities. The white and red boxes represent R
functions/methods and objects respectively. The blue boxes represent
different disk storage formats.](Figures/MSnbase-io-in.png)

*MSnbase* input capabilities. The white and red boxes represent R
functions/methods and objects respectively. The blue boxes represent
different disk storage formats.

## Data output

##### RData files

R objects can most easily be stored on disk with the `save` function. It
creates compressed binary images of the data representation that can
later be read back from the file with the `load` function.

##### mzML/mzXML files

`MSnExp` and `OnDiskMSnExp` files can be written to MS data files in
`mzML` or `mzXML` files with the `writeMSData` method. See
[`?writeMSData`](https://lgatto.github.io/MSnbase/reference/writeMSData.md)
for details.

##### Peak lists

`MSnExp` instances as well as individual spectra can be written as `mgf`
files with the `writeMgfData` method. Note that the meta-data in the
original R object can not be included in the file. See
[`?writeMgfData`](https://lgatto.github.io/MSnbase/reference/writeMgfData-methods.md)
for details.

##### Quantitation data

Quantitation data can be exported to spreadsheet files with the
`write.exprs` method. Feature meta-data can be appended to the feature
intensity values. See
[`?writeMgfData`](https://lgatto.github.io/MSnbase/reference/writeMgfData-methods.md)
for details.

**Deprecated** `MSnSet` instances can also be exported to `mzTab` files
using the `writeMzTabData` function.

![MSnbase output capabilities. The white and red boxes represent R
functions/methods and objects respectively. The blue boxes represent
different disk storage formats.](Figures/MSnbase-io-out.png)

*MSnbase* output capabilities. The white and red boxes represent R
functions/methods and objects respectively. The blue boxes represent
different disk storage formats.

## Creating `MSnSet` from text spread sheets

This section describes the generation of `MSnSet` objects using data
available in a text-based spreadsheet. This entry point into R and
*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* allows to
import data processed by any of the third party mass-spectrometry
processing software available and proceed with data exploration,
normalisation and statistical analysis using functions available in and
the numerous Bioconductor packages.

### A complete work flow

The following section describes a work flow that uses three input files
to create the `MSnSet`. These files respectively describe the
quantitative expression data, the sample meta-data and the feature
meta-data. It is taken from the
*[pRoloc](https://bioconductor.org/packages/3.23/pRoloc)* tutorial and
uses example files from the
*[pRolocdat](https://bioconductor.org/packages/3.23/pRolocdat)* package.

We start by describing the `csv` to be used as input using the
`read.csv` function.

``` r

## The original data for replicate 1, available
## from the pRolocdata package
f0 <- dir(system.file("extdata", package = "pRolocdata"),
          full.names = TRUE,
          pattern = "pr800866n_si_004-rep1.csv")
csv <- read.csv(f0)
```

The three first lines of the original spreadsheet, containing the data
for replicate one, are illustrated below (using the function `head`). It
contains 888 rows (proteins) and 16 columns, including protein
identifiers, database accession numbers, gene symbols, reporter ion
quantitation values, information related to protein identification, …

``` r

head(csv, n=3)
```

    ##   Protein.ID        FBgn Flybase.Symbol No..peptide.IDs Mascot.score
    ## 1    CG10060 FBgn0001104    G-ialpha65A               3       179.86
    ## 2    CG10067 FBgn0000044         Act57B               5       222.40
    ## 3    CG10077 FBgn0035720        CG10077               5       219.65
    ##   No..peptides.quantified area.114 area.115 area.116 area.117
    ## 1                       1 0.379000 0.281000 0.225000 0.114000
    ## 2                       9 0.420000 0.209667 0.206111 0.163889
    ## 3                       3 0.187333 0.167333 0.169667 0.476000
    ##   PLS.DA.classification Peptide.sequence Precursor.ion.mass
    ## 1                    PM                                    
    ## 2                    PM                                    
    ## 3                                                          
    ##   Precursor.ion.charge pd.2013 pd.markers
    ## 1                           PM    unknown
    ## 2                           PM    unknown
    ## 3                      unknown    unknown

Below read in turn the spread sheets that contain the quantitation data
(`exprsFile.csv`), feature meta-data (`fdataFile.csv`) and sample
meta-data (`pdataFile.csv`).

``` r

## The quantitation data, from the original data
f1 <- dir(system.file("extdata", package = "pRolocdata"),
          full.names = TRUE, pattern = "exprsFile.csv")
exprsCsv <- read.csv(f1)
## Feature meta-data, from the original data
f2 <- dir(system.file("extdata", package = "pRolocdata"),
          full.names = TRUE, pattern = "fdataFile.csv")
fdataCsv <- read.csv(f2)
## Sample meta-data, a new file
f3 <- dir(system.file("extdata", package = "pRolocdata"),
          full.names = TRUE, pattern = "pdataFile.csv")
pdataCsv <- read.csv(f3)
```

`exprsFile.csv` contains the quantitation (expression) data for the 888
proteins and 4 reporter tags.

``` r

head(exprsCsv, n = 3)
```

    ##          FBgn     X114     X115     X116     X117
    ## 1 FBgn0001104 0.379000 0.281000 0.225000 0.114000
    ## 2 FBgn0000044 0.420000 0.209667 0.206111 0.163889
    ## 3 FBgn0035720 0.187333 0.167333 0.169667 0.476000

`fdataFile.csv` contains meta-data for the 888 features (here proteins).

``` r

head(fdataCsv, n = 3)
```

    ##          FBgn ProteinID FlybaseSymbol NoPeptideIDs MascotScore
    ## 1 FBgn0001104   CG10060   G-ialpha65A            3      179.86
    ## 2 FBgn0000044   CG10067        Act57B            5      222.40
    ## 3 FBgn0035720   CG10077       CG10077            5      219.65
    ##   NoPeptidesQuantified PLSDA
    ## 1                    1    PM
    ## 2                    9    PM
    ## 3                    3

`pdataFile.csv` contains samples (here fractions) meta-data. This simple
file has been created manually.

``` r

pdataCsv
```

    ##   sampleNames Fractions
    ## 1        X114       4/5
    ## 2        X115     12/13
    ## 3        X116        19
    ## 4        X117        21

The self-contained `MSnSet` can now easily be generated using the
`readMSnSet` constructor, providing the respective `csv` file names
shown above and specifying that the data is comma-separated (with
`sep = ","`). Below, we call that object `res` and display its content.

``` r

library("MSnbase")
res <- readMSnSet(exprsFile = f1,
                  featureDataFile = f2,
                  phenoDataFile = f3,
                  sep = ",")
res
```

    ## MSnSet (storageMode: lockedEnvironment)
    ## assayData: 888 features, 4 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: X114 X115 X116 X117
    ##   varLabels: Fractions
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: FBgn0001104 FBgn0000044 ... FBgn0001215 (888 total)
    ##   fvarLabels: ProteinID FlybaseSymbol ... PLSDA (6 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:  
    ## - - - Processing information - - -
    ## Quantitation data loaded: Wed Mar 18 17:50:47 2026  using readMSnSet. 
    ##  MSnbase version: 2.37.2

#### The `MSnSet` class

Although there are additional specific sub-containers for additional
meta-data (for instance to make the object MIAPE compliant), the feature
(the sub-container, or slot `featureData`) and sample (the `phenoData`
slot) are the most important ones. They need to meet the following
validity requirements (see figure below):

- the number of row in the expression/quantitation data and feature data
  must be equal and the row names must match exactly, and

- the number of columns in the expression/quantitation data and number
  of row in the sample meta-data must be equal and the column/row names
  must match exactly.

A detailed description of the `MSnSet` class is available by typing
[`?MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
in the R console.

![Dimension requirements for the respective expression, feature and
sample meta-data slots.](Figures/msnset.png)

Dimension requirements for the respective expression, feature and sample
meta-data slots.

The individual parts of this data object can be accessed with their
respective accessor methods:

- the quantitation data can be retrieved with `exprs(res)`,
- the feature meta-data with `fData(res)` and
- the sample meta-data with `pData(res)`.

### A shorter work flow

The `readMSnSet2` function provides a simplified import workforce. It
takes a single spreadsheet as input (default is `csv`) and extract the
columns identified by `ecol` to create the expression data, while the
others are used as feature meta-data. `ecol` can be a `character` with
the respective column labels or a numeric with their indices. In the
former case, it is important to make sure that the names match exactly.
Special characters like `'-'` or `'('` will be transformed by R into
`'.'` when the `csv` file is read in. Optionally, one can also specify a
column to be used as feature names. Note that these must be unique to
guarantee the final object validity.

``` r

ecol <- paste("area", 114:117, sep = ".")
fname <- "Protein.ID"
eset <- readMSnSet2(f0, ecol, fname)
eset
```

    ## MSnSet (storageMode: lockedEnvironment)
    ## assayData: 888 features, 4 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData: none
    ## featureData
    ##   featureNames: CG10060 CG10067 ... CG9983 (888 total)
    ##   fvarLabels: Protein.ID FBgn ... pd.markers (12 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:  
    ## - - - Processing information - - -
    ##  MSnbase version: 2.37.2

The `ecol` columns can also be queried interactively from R using the
`getEcols` and `grepEcols` function. The former return a character with
all column names, given a splitting character, i.e. the separation value
of the spreadsheet (typically `","` for `csv`, `"\t"` for `tsv`, …). The
latter can be used to grep a pattern of interest to obtain the relevant
column indices.

``` r

getEcols(f0, ",")
```

    ##  [1] "\"Protein ID\""              "\"FBgn\""                   
    ##  [3] "\"Flybase Symbol\""          "\"No. peptide IDs\""        
    ##  [5] "\"Mascot score\""            "\"No. peptides quantified\""
    ##  [7] "\"area 114\""                "\"area 115\""               
    ##  [9] "\"area 116\""                "\"area 117\""               
    ## [11] "\"PLS-DA classification\""   "\"Peptide sequence\""       
    ## [13] "\"Precursor ion mass\""      "\"Precursor ion charge\""   
    ## [15] "\"pd.2013\""                 "\"pd.markers\""

``` r

grepEcols(f0, "area", ",")
```

    ## [1]  7  8  9 10

``` r

e <- grepEcols(f0, "area", ",")
readMSnSet2(f0, e)
```

    ## MSnSet (storageMode: lockedEnvironment)
    ## assayData: 888 features, 4 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData: none
    ## featureData
    ##   featureNames: 1 2 ... 888 (888 total)
    ##   fvarLabels: Protein.ID FBgn ... pd.markers (12 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:  
    ## - - - Processing information - - -
    ##  MSnbase version: 2.37.2

The `phenoData` slot can now be updated accordingly using the
replacement functions `phenoData<-` or `pData<-` (see
[`?MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
for details).

## Session information

``` r

sessionInfo()
```

    ## R Under development (unstable) (2026-03-15 r89629)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] pRolocdata_1.49.0   MSnbase_2.37.2      ProtGenerics_1.43.0
    ##  [4] S4Vectors_0.49.0    mzR_2.45.0          Rcpp_1.1.1         
    ##  [7] Biobase_2.71.0      BiocGenerics_0.57.0 generics_0.1.4     
    ## [10] BiocStyle_2.39.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rlang_1.1.7                 magrittr_2.0.4             
    ##  [3] clue_0.3-67                 otel_0.2.0                 
    ##  [5] matrixStats_1.5.0           compiler_4.6.0             
    ##  [7] systemfonts_1.3.2           vctrs_0.7.1                
    ##  [9] reshape2_1.4.5              stringr_1.6.0              
    ## [11] pkgconfig_2.0.3             MetaboCoreUtils_1.19.2     
    ## [13] fastmap_1.2.0               XVector_0.51.0             
    ## [15] rmarkdown_2.30              preprocessCore_1.73.0      
    ## [17] ragg_1.5.1                  purrr_1.2.1                
    ## [19] xfun_0.56                   MultiAssayExperiment_1.37.2
    ## [21] cachem_1.1.0                jsonlite_2.0.0             
    ## [23] DelayedArray_0.37.0         BiocParallel_1.45.0        
    ## [25] parallel_4.6.0              cluster_2.1.8.2            
    ## [27] R6_2.6.1                    bslib_0.10.0               
    ## [29] stringi_1.8.7               RColorBrewer_1.1-3         
    ## [31] limma_3.67.0                GenomicRanges_1.63.1       
    ## [33] jquerylib_0.1.4             iterators_1.0.14           
    ## [35] Seqinfo_1.1.0               bookdown_0.46              
    ## [37] SummarizedExperiment_1.41.1 knitr_1.51                 
    ## [39] IRanges_2.45.0              Matrix_1.7-4               
    ## [41] igraph_2.2.2                tidyselect_1.2.1           
    ## [43] abind_1.4-8                 yaml_2.3.12                
    ## [45] doParallel_1.0.17           codetools_0.2-20           
    ## [47] affy_1.89.0                 lattice_0.22-9             
    ## [49] tibble_3.3.1                plyr_1.8.9                 
    ## [51] S7_0.2.1                    evaluate_1.0.5             
    ## [53] desc_1.4.3                  Spectra_1.21.5             
    ## [55] pillar_1.11.1               affyio_1.81.0              
    ## [57] BiocManager_1.30.27         MatrixGenerics_1.23.0      
    ## [59] foreach_1.5.2               MALDIquant_1.22.3          
    ## [61] ncdf4_1.24                  ggplot2_4.0.2              
    ## [63] scales_1.4.0                glue_1.8.0                 
    ## [65] lazyeval_0.2.2              tools_4.6.0                
    ## [67] mzID_1.49.0                 data.table_1.18.2.1        
    ## [69] QFeatures_1.21.0            vsn_3.79.5                 
    ## [71] fs_1.6.7                    XML_3.99-0.22              
    ## [73] grid_4.6.0                  impute_1.85.0              
    ## [75] tidyr_1.3.2                 MsCoreUtils_1.23.6         
    ## [77] PSMatch_1.15.1              cli_3.6.5                  
    ## [79] textshaping_1.0.5           S4Arrays_1.11.1            
    ## [81] dplyr_1.2.0                 AnnotationFilter_1.35.0    
    ## [83] pcaMethods_2.3.0            gtable_0.3.6               
    ## [85] sass_0.4.10                 digest_0.6.39              
    ## [87] SparseArray_1.11.11         htmlwidgets_1.6.4          
    ## [89] farver_2.1.2                htmltools_0.5.9            
    ## [91] pkgdown_2.2.0.9000          lifecycle_1.0.5            
    ## [93] statmod_1.5.1               MASS_7.3-65

## References

Gatto, Laurent, and Kathryn S Lilley. 2012. “MSnbase – an R/Bioconductor
Package for Isobaric Tagged Mass Spectrometry Data Visualization,
Processing and Quantitation.” *Bioinformatics* 28 (2): 288–89.
<https://doi.org/10.1093/bioinformatics/btr645>.

Martens, Lennart, Matthew Chambers, Marc Sturm, et al. 2010. “mzML - a
Community Standard for Mass Spectrometry Data.” *Molecular & Cellular
Proteomics : MCP*, ahead of print.
<https://doi.org/10.1074/mcp.R110.000133>.

Orchard, Sandra, Luisa Montechi-Palazzi, Eric W Deutsch, et al. 2007.
“Five Years of Progress in the Standardization of Proteomics Data 4th
Annual Spring Workshop of the HUPO-Proteomics Standards Initiative April
23-25, 2007 Ecole Nationale Supérieure (ENS), Lyon, France.”
*Proteomics* 7 (19): 3436–40. <https://doi.org/10.1002/pmic.200700658>.

Pedrioli, Patrick G A, Jimmy K Eng, Robert Hubley, et al. 2004. “A
Common Open Representation of Mass Spectrometry Data and Its Application
to Proteomics Research.” *Nat. Biotechnol.* 22 (11): 1459–66.
<https://doi.org/10.1038/nbt1031>.

[^1]: [http://www.matrixscience.com/help/data_file_help.html](http://www.matrixscience.com/help/data_file_help.md)

[^2]: <https://github.com/HUPO-PSI/mzTab>
