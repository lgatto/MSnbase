# A short introduction to \*MSnbase\* development

Abstract

This vignette describes the classes implemented in package. It is
intended as a starting point for developers or users who would like to
learn more or further develop/extend mass spectrometry and proteomics
data structures.

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

**NB** This document is going to be updated based on current major
development plans in `MSnbase`.

## Introduction

This document is not a replacement for the individual manual pages, that
document the slots of the
*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* classes. It
is a centralised high-level description of the package design.

*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* aims at
being compatible with the
*[Biobase](https://bioconductor.org/packages/3.23/Biobase)*
infrastructure (Gentleman et al. 2004). Many meta data structures that
are used in `eSet` and associated classes are also used here. As such,
knowledge of the *Biobase development and the new eSet* vignette would
be beneficial; the vignette can directly be accessed with
[`vignette("BiobaseDevelopment", package="Biobase")`](https://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/BiobaseDevelopment.html).

The initial goal is to use the
*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)*
infrastructure for MS2 labelled (iTRAQ (Ross et al. 2004) and TMT
(Thompson et al. 2003)) and label-free (spectral counting, index and
abundance) quantitation - see the documentation for the `quantify`
function for details. The infrastructure is currently extended to
support a wider range of technologies, including metabolomics.

## Coding style

`MSnbase` follows the [Bioconductor style
guide](https://bioconductor.org/developers/how-to/coding-style/). In
particular

- Do not use `.` when naming symbols.
- A leading `.` can be used for hidden/local functions or variables.
- Snake case should be restricted to internal functions. For
  consistency, we favour camel case for public functions.
- Class names should start with a capital and each class should posses a
  constructor with identical name. Running the constructor without any
  input should produce a valid empty object.
- Use `##` to start full-line comments.
- For roxygen headers `#'` is preferred, although `##'` is tolerated.
- Use spaces between `=` in function arguments or class definition:
  `f(a = 1, b = 2)`.
- Always use a space after a comma: `a, b, c`.
- Always use spaces around binary operators: `a + b`.
- Lines should be kept shorter than 80 characters. For example the
  following code isn’t accepted

``` r

# no wrap at 80
someVeryLongVariableName <- someVeryLongFunctionName(withSomeEvenLongerFunctionArgumentA = 1, withSomeEvenLongerFunctionArgumentB = 2)
```

and should be wrapped as shown below:

``` r

# alternative 1
someVeryLongVariableName <-
    someVeryLongFunctionName(withSomeEvenLongerFunctionArgumentA = 1,
                             withSomeEvenLongerFunctionArgumentB = 2)

# alternative 2
someVeryLongVariableName <- someVeryLongFunctionName(
    withSomeEvenLongerFunctionArgumentA = 1,
    withSomeEvenLongerFunctionArgumentB = 2)
```

## *[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* classes

All classes have a `.__classVersion__` slot, of class `Versioned` from
the *[Biobase](https://bioconductor.org/packages/3.23/Biobase)* package.
This slot documents the class version for any instance to be used for
debugging and object update purposes. Any change in a class
implementation should trigger a version change.

### `pSet`: a virtual class for raw mass spectrometry data and meta data

This virtual class is the main container for mass spectrometry data, i.e
spectra, and meta data. It is based on the `eSet` implementation for
genomic data. The main difference with `eSet` is that the `assayData`
slot is an environment containing any number of `Spectrum` instances
(see the [`Spectrum` section](#Spectrum)).

One new slot is introduced, namely `processingData`, that contains one
`MSnProcess` instance (see the [`MSnProcess` section](#MSnProcess)). and
the `experimentData` slot is now expected to contain `MIAPE` data. The
`annotation` slot has not been implemented, as no prior feature
annotation is known in shotgun proteomics.

``` r

getClass("pSet")
```

    ## Virtual Class "pSet" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                
    ## Name:           assayData          phenoData        featureData
    ## Class:        environment AnnotatedDataFrame AnnotatedDataFrame
    ##                                                                
    ## Name:      experimentData       protocolData     processingData
    ## Class:              MIAxE AnnotatedDataFrame         MSnProcess
    ##                                             
    ## Name:              .cache  .__classVersion__
    ## Class:        environment           Versions
    ## 
    ## Extends: "Versioned"
    ## 
    ## Known Subclasses: 
    ## Class "MSnExp", directly
    ## Class "OnDiskMSnExp", by class "MSnExp", distance 2, with explicit coerce

### `MSnExp`: a class for MS experiments

`MSnExp` extends `pSet` to store MS experiments. It does not add any new
slots to `pSet`. Accessors and setters are all inherited from `pSet` and
new ones should be implemented for `pSet`. Methods that manipulate
actual data in experiments are implemented for `MSnExp` objects.

``` r

getClass("MSnExp")
```

    ## Class "MSnExp" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                
    ## Name:           assayData          phenoData        featureData
    ## Class:        environment AnnotatedDataFrame AnnotatedDataFrame
    ##                                                                
    ## Name:      experimentData       protocolData     processingData
    ## Class:              MIAxE AnnotatedDataFrame         MSnProcess
    ##                                             
    ## Name:              .cache  .__classVersion__
    ## Class:        environment           Versions
    ## 
    ## Extends: 
    ## Class "pSet", directly
    ## Class "Versioned", by class "pSet", distance 2
    ## 
    ## Known Subclasses: 
    ## Class "OnDiskMSnExp", directly, with explicit coerce

### `OnDiskMSnExp`: a on-disk implementation of the `MSnExp` class

The `OnDiskMSnExp` class extends `MSnExp` and inherits all of its
functionality but is aimed to use as little memory as possible based on
a balance between memory demand and performance. Most of the
spectrum-specific data, like retention time, polarity, total ion current
are stored within the object’s `featureData` slot. The actual M/Z and
intensity values from the individual spectra are, in contrast to
`MSnExp` objects, not kept in memory (in the `assayData` slot), but are
fetched from the original files on-demand. Because mzML files are
indexed, using the *[mzR](https://bioconductor.org/packages/3.23/mzR)*
package to read the relevant spectrum data is fast and only moderately
slower than for in-memory `MSnExp`[^1].

To keep track of data manipulation steps that are applied to spectrum
data (such as performed by methods `removePeaks` or `clean`) a *lazy
execution* framework was implemented. Methods that manipulate or subset
a spectrum’s M/Z or intensity values can not be applied directly to a
`OnDiskMSnExp` object, since the relevant data is not kept in memory.
Thus, any call to a processing method that changes or subset M/Z or
intensity values are added as `ProcessingStep` items to the object’s
`spectraProcessingQueue`. When the spectrum data is then queried from an
`OnDiskMSnExp`, the spectra are read in from the file and all these
processing steps are applied on-the-fly to the spectrum data before
being returned to the user.

The operations involving extracting or manipulating spectrum data are
applied on a per-file basis, which enables parallel processing. Thus,
all corresponding method implementations for `OnDiskMSnExp` objects have
an argument `BPPARAM` and users can set a `PARALLEL_THRESH` option
flag[^2] that enables to define how and when parallel processing should
be performed (using the
*[BiocParallel](https://bioconductor.org/packages/3.23/BiocParallel)*
package).

Note that all data manipulations that are not applied to M/Z or
intensity values of a spectrum (e.g. sub-setting by retention time etc)
are very fast as they operate directly to the object’s `featureData`
slot.

``` r

getClass("OnDiskMSnExp")
```

    ## Class "OnDiskMSnExp" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                            
    ## Name:  spectraProcessingQueue                backend              assayData
    ## Class:                   list              character            environment
    ##                                                                            
    ## Name:               phenoData            featureData         experimentData
    ## Class:     AnnotatedDataFrame     AnnotatedDataFrame                  MIAxE
    ##                                                                            
    ## Name:            protocolData         processingData                 .cache
    ## Class:     AnnotatedDataFrame             MSnProcess            environment
    ##                              
    ## Name:       .__classVersion__
    ## Class:               Versions
    ## 
    ## Extends: 
    ## Class "MSnExp", directly
    ## Class "pSet", by class "MSnExp", distance 2
    ## Class "Versioned", by class "MSnExp", distance 3

The distinction between `MSnExp` and `OnDiskMSnExp` is often not
explicitly stated as it should not matter, from a user’s perspective,
which data structure they are working with, as both behave in equivalent
ways. Often, they are referred to as *in-memory* and *on-disk* `MSnExp`
implementations.

### `MSnSet`: a class for quantitative proteomics data

This class stores quantitation data and meta data after running
`quantify` on an `MSnExp` object or by creating an `MSnSet` instance
from an external file, as described in the *MSnbase-io* vignette and in
[`?readMSnSet`](https://lgatto.github.io/MSnbase/reference/readMSnSet.md),
`readMzTabData`, etc. The quantitative data is in form of a *n* by *p*
matrix, where *n* is the number of features/spectra originally in the
`MSnExp` used as parameter in `quantify` and *p* is the number of
reporter ions. If read from an external file, *n* corresponds to the
number of features (protein groups, proteins, peptides, spectra) in the
file and $`p`$ is the number of columns with quantitative data (samples)
in the file.

This prompted to keep a similar implementation as the `ExpressionSet`
class, while adding the proteomics-specific annotation slot introduced
in the `pSet` class, namely `processingData` for objects of class
`MSnProcess`.

``` r

getClass("MSnSet")
```

    ## Class "MSnSet" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                
    ## Name:      experimentData     processingData               qual
    ## Class:              MIAPE         MSnProcess         data.frame
    ##                                                                
    ## Name:           assayData          phenoData        featureData
    ## Class:          AssayData AnnotatedDataFrame AnnotatedDataFrame
    ##                                                                
    ## Name:          annotation       protocolData  .__classVersion__
    ## Class:          character AnnotatedDataFrame           Versions
    ## 
    ## Extends: 
    ## Class "eSet", directly
    ## Class "VersionedBiobase", by class "eSet", distance 2
    ## Class "Versioned", by class "eSet", distance 3

The `MSnSet` class extends the virtual `eSet` class to provide
compatibility for `ExpressionSet`-like behaviour. The experiment
meta-data in `experimentData` is also of class `MIAPE` . The
`annotation` slot, inherited from `eSet` is not used. As a result, it is
easy to convert `ExpressionSet` data from/to `MSnSet` objects with the
coersion method `as`.

``` r

data(msnset)
class(msnset)
```

    ## [1] "MSnSet"
    ## attr(,"package")
    ## [1] "MSnbase"

``` r

class(as(msnset, "ExpressionSet"))
```

    ## [1] "ExpressionSet"
    ## attr(,"package")
    ## [1] "Biobase"

``` r

data(sample.ExpressionSet)
class(sample.ExpressionSet)
```

    ## [1] "ExpressionSet"
    ## attr(,"package")
    ## [1] "Biobase"

``` r

class(as(sample.ExpressionSet, "MSnSet"))
```

    ## [1] "MSnSet"
    ## attr(,"package")
    ## [1] "MSnbase"

### `MSnProcess`: a class for logging processing meta data

This class aims at recording specific manipulations applied to `MSnExp`
or `MSnSet` instances. The `processing` slot is a `character` vector
that describes major processing. Most other slots are of class `logical`
that indicate whether the data has been centroided, smoothed, although
many of the functionality is not implemented yet. Any new processing
that is implemented should be documented and logged here.

It also documents the raw data file from which the data originates
(`files` slot) and the
*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* version that
was in use when the `MSnProcess` instance, and hence the
`MSnExp`/`MSnSet` objects, were originally created.

``` r

getClass("MSnProcess")
```

    ## Class "MSnProcess" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                               
    ## Name:              files        processing            merged           cleaned
    ## Class:         character         character           logical           logical
    ##                                                                               
    ## Name:       removedPeaks          smoothed           trimmed        normalised
    ## Class:         character           logical           numeric           logical
    ##                                           
    ## Name:     MSnbaseVersion .__classVersion__
    ## Class:         character          Versions
    ## 
    ## Extends: "Versioned"

### `MIAPE`: Minimum Information About a Proteomics Experiment

The Minimum Information About a Proteomics Experiment (Taylor et al.
2007; Taylor et al. 2008) `MIAPE` class describes the experiment,
including contact details, information about the mass spectrometer and
control and analysis software.

``` r

getClass("MIAPE")
```

    ## Class "MIAPE" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                         
    ## Name:                     title                      url
    ## Class:                character                character
    ##                                                         
    ## Name:                  abstract                pubMedIds
    ## Class:                character                character
    ##                                                         
    ## Name:                   samples            preprocessing
    ## Class:                     list                     list
    ##                                                         
    ## Name:                     other                dateStamp
    ## Class:                     list                character
    ##                                                         
    ## Name:                      name                      lab
    ## Class:                character                character
    ##                                                         
    ## Name:                   contact                    email
    ## Class:                character                character
    ##                                                         
    ## Name:           instrumentModel   instrumentManufacturer
    ## Class:                character                character
    ##                                                         
    ## Name:  instrumentCustomisations             softwareName
    ## Class:                character                character
    ##                                                         
    ## Name:           softwareVersion        switchingCriteria
    ## Class:                character                character
    ##                                                         
    ## Name:            isolationWidth            parameterFile
    ## Class:                  numeric                character
    ##                                                         
    ## Name:                 ionSource         ionSourceDetails
    ## Class:                character                character
    ##                                                         
    ## Name:                  analyser          analyserDetails
    ## Class:                character                character
    ##                                                         
    ## Name:              collisionGas        collisionPressure
    ## Class:                character                  numeric
    ##                                                         
    ## Name:           collisionEnergy             detectorType
    ## Class:                character                character
    ##                                                         
    ## Name:       detectorSensitivity        .__classVersion__
    ## Class:                character                 Versions
    ## 
    ## Extends: 
    ## Class "MIAxE", directly
    ## Class "Versioned", by class "MIAxE", distance 2

### `Spectrum` *et al.*: classes for MS spectra

`Spectrum` is a virtual class that defines common attributes to all
types of spectra. MS1 and MS2 specific attributes are defined in the
`Spectrum1` and `Spectrum2` classes, that directly extend `Spectrum`.

``` r

getClass("Spectrum")
```

    ## Virtual Class "Spectrum" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                               
    ## Name:            msLevel        peaksCount                rt    acquisitionNum
    ## Class:           integer           integer           numeric           integer
    ##                                                                               
    ## Name:          scanIndex               tic                mz         intensity
    ## Class:           integer           numeric           numeric           numeric
    ##                                                                               
    ## Name:           fromFile        centroided          smoothed          polarity
    ## Class:           integer           logical           logical           integer
    ##                         
    ## Name:  .__classVersion__
    ## Class:          Versions
    ## 
    ## Extends: "Versioned"
    ## 
    ## Known Subclasses: "Spectrum2", "Spectrum1"

``` r

getClass("Spectrum1")
```

    ## Class "Spectrum1" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                               
    ## Name:            msLevel        peaksCount                rt    acquisitionNum
    ## Class:           integer           integer           numeric           integer
    ##                                                                               
    ## Name:          scanIndex               tic                mz         intensity
    ## Class:           integer           numeric           numeric           numeric
    ##                                                                               
    ## Name:           fromFile        centroided          smoothed          polarity
    ## Class:           integer           logical           logical           integer
    ##                         
    ## Name:  .__classVersion__
    ## Class:          Versions
    ## 
    ## Extends: 
    ## Class "Spectrum", directly
    ## Class "Versioned", by class "Spectrum", distance 2

``` r

getClass("Spectrum2")
```

    ## Class "Spectrum2" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                
    ## Name:              merged        precScanNum        precursorMz
    ## Class:            numeric            integer            numeric
    ##                                                                
    ## Name:  precursorIntensity    precursorCharge    collisionEnergy
    ## Class:            numeric            integer            numeric
    ##                                                                
    ## Name:             msLevel         peaksCount                 rt
    ## Class:            integer            integer            numeric
    ##                                                                
    ## Name:      acquisitionNum          scanIndex                tic
    ## Class:            integer            integer            numeric
    ##                                                                
    ## Name:                  mz          intensity           fromFile
    ## Class:            numeric            numeric            integer
    ##                                                                
    ## Name:          centroided           smoothed           polarity
    ## Class:            logical            logical            integer
    ##                          
    ## Name:   .__classVersion__
    ## Class:           Versions
    ## 
    ## Extends: 
    ## Class "Spectrum", directly
    ## Class "Versioned", by class "Spectrum", distance 2

### `ReporterIons`: a class for isobaric tags

The iTRAQ and TMT (or any other peak of interest) are implemented
`ReporterIons` instances, that essentially defines an expected MZ
position for the peak and a width around this value as well a names for
the reporters.

``` r

getClass("ReporterIons")
```

    ## Class "ReporterIons" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                               
    ## Name:               name     reporterNames       description                mz
    ## Class:         character         character         character           numeric
    ##                                                             
    ## Name:                col             width .__classVersion__
    ## Class:         character           numeric          Versions
    ## 
    ## Extends: "Versioned"

### `Chromatogram` and `MChromatograms`: classes to handle chromatographic data

The `Chromatogram` class represents chromatographic MS data,
i.e. retention time and intensity duplets for one file/sample. The
`MChromatograms` class (Matrix of Chromatograms) allows to arrange
multiple `Chromatogram` instances in a two-dimensional grid, with
columns supposed to represent different samples and rows two-dimensional
areas in the plane spanned by the m/z and retention time dimensions from
which the intensities are extracted (e.g. an extracted ion chromatogram
for a specific ion). The `MChromatograms` class extends the base
`matrix` class. `MChromatograms` objects can be extracted from an
`MSnExp` or `OnDiskMSnExp` object using the `chromatogram` method.

``` r

getClass("Chromatogram")
```

    ## Class "Chromatogram" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                               
    ## Name:              rtime         intensity                mz          filterMz
    ## Class:           numeric           numeric           numeric           numeric
    ##                                                                               
    ## Name:        precursorMz         productMz          fromFile    aggregationFun
    ## Class:           numeric           numeric           integer         character
    ##                                           
    ## Name:            msLevel .__classVersion__
    ## Class:           integer          Versions
    ## 
    ## Extends: "Versioned"

``` r

getClass("MChromatograms")
```

    ## Class "MChromatograms" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                
    ## Name:               .Data          phenoData        featureData
    ## Class:             matrix AnnotatedDataFrame AnnotatedDataFrame
    ## 
    ## Extends: 
    ## Class "matrix", from data part
    ## Class "array", by class "matrix", distance 2
    ## Class "dataframeOrDataFrameOrmatrix", by class "matrix", distance 2
    ## Class "structure", by class "matrix", distance 3
    ## Class "matrix_OR_array_OR_table_OR_numeric", by class "matrix", distance 3
    ## Class "vector", by class "matrix", distance 4, with explicit coerce
    ## Class "vector_OR_factor", by class "matrix", distance 5, with explicit coerce
    ## Class "vector_OR_Vector", by class "matrix", distance 5, with explicit coerce

### Other classes

#### Lists of `MSnSet` instances

When several `MSnSet` instances are related to each other and should be
stored together as different objects, they can be grouped as a list into
and `MSnSetList` object. In addition to the actual `list` slot, this
class also has basic logging functionality and enables iteration over
the `MSnSet` instances using a dedicated `lapply` methods.

``` r

getClass("MSnSetList")
```

    ## Class "MSnSetList" [package "MSnbase"]
    ## 
    ## Slots:
    ##                                                                               
    ## Name:                  x               log       featureData .__classVersion__
    ## Class:              list              list         DataFrame          Versions
    ## 
    ## Extends: "Versioned"

## Miscellaneous

##### Unit tests

*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* implements
unit tests with the
*[testthat](https://CRAN.R-project.org/package=testthat)* package.

##### Processing methods

Methods that process raw data, i.e. spectra should be implemented for
`Spectrum` objects first and then `eapply`ed (or similar) to the
`assayData` slot of an `MSnExp` instance in the specific method.

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
    ## [1] MSnbase_2.37.1      ProtGenerics_1.43.0 S4Vectors_0.49.0   
    ## [4] mzR_2.45.0          Rcpp_1.1.1          Biobase_2.71.0     
    ## [7] BiocGenerics_0.57.0 generics_0.1.4      BiocStyle_2.39.0   
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

Gentleman, Robert C., Vincent J. Carey, Douglas M. Bates, et al. 2004.
“Bioconductor: Open Software Development for Computational Biology and
Bioinformatics.” *Genome Biol* 5 (10): –80.
<https://doi.org/10.1186/gb-2004-5-10-r80>.

Ross, Philip L., Yulin N. Huang, Jason N. Marchese, et al. 2004.
“Multiplexed Protein Quantitation in Saccharomyces Cerevisiae Using
Amine-Reactive Isobaric Tagging Reagents.” *Mol Cell Proteomics* 3 (12):
1154–69. <https://doi.org/10.1074/mcp.M400129-MCP200>.

Taylor, Chris F, Pierre-Alain Binz, Ruedi Aebersold, et al. 2008.
“Guidelines for Reporting the Use of Mass Spectrometry in Proteomics.”
*Nat. Biotechnol.* 26 (8): 860–61.
<https://doi.org/10.1038/nbt0808-860>.

Taylor, Chris F., Norman W. Paton, Kathryn S. Lilley, et al. 2007. “The
Minimum Information about a Proteomics Experiment (MIAPE).” *Nat
Biotechnol* 25 (8): 887–93. <https://doi.org/10.1038/nbt1329>.

Thompson, Andrew, Jürgen Schäfer, Karsten Kuhn, et al. 2003. “[Tandem
Mass Tags: A Novel Quantification Strategy for Comparative Analysis of
Complex Protein Mixtures by
MS/MS.](https://www.ncbi.nlm.nih.gov/pubmed/12713048)” *Anal. Chem.* 75
(8): 1895–904.

[^1]: The *benchmarking* vignette compares data size and operation speed
    of the two implementations.

[^2]: see
    [`?MSnbaseOptions`](https://lgatto.github.io/MSnbase/reference/MSnbaseOptions.md)
    for details.
