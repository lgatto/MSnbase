[![Build Status](https://travis-ci.org/lgatto/MSnbase.svg?branch=master)](https://travis-ci.org/lgatto/MSnbase) 
[![codecov.io](https://codecov.io/github/lgatto/MSnbase/coverage.svg?branch=master)](https://codecov.io/github/lgatto/MSnbase?branch=master)




# The `MSnbase` package

<img align = "right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/MSnbase/MSnbase.png" height="200">

[**MSnbase**](http://lgatto.github.io/MSnbase/) is an R/Bioconductor
package that provides infrastructure for plotting, manipulation and
processing mass spectrometry and proteomics data. The project was
started by [Laurent Gatto](http://lgatto.github.io/) in October 2010
(Mon Oct 4 23:35:23 2010, according to the git log) and has, since
then, benefited from
[various contributions](http://lgatto.github.io/msnbase-contribs/), in
particular [Sebastian Gibb](http://sebastiangibb.de/)
and [Johannes Rainer](https://github.com/jotsetung/).

The official page is the Bioconductor landing page
([release](http://www.bioconductor.org/packages/release/bioc/html/MSnbase.html) or
[devel](http://www.bioconductor.org/packages/devel/bioc/html/MSnbase.html) versions). The
[github page](https://github.com/lgatto/MSnbase) page is for active
development, issue tracking and forking/pulling purposes. 

To get an overview of the package, see the
[*MSnbase-demo*](https://lgatto.github.io/MSnbase/articles/MSnbase-demo.html) vignette. More
vignettes are available in the *Articles* tab. 

## Installation

To install *[MSnbase](http://bioconductor.org/packages/MSnbase)*


```r
library("BiocInstaller")
biocLite("MSnbase")
```

If you need the github version (not recommended unless you know what
you are doing), use


```r
biocLite("lgatto/MSnbase")
```

## Questions

General questions should be asked on
the [Bioconductor support forum](https://support.bioconductor.org/),
using `MSnbase` to tag the question. Feel also free to open a
GitHub [issue](https://github.com/lgatto/MSnbase/issues), in
particular for bug reports.

## Citation

To cite the `MSnbase` package in publications, please use:

> Gatto L, Lilley KS. *`MSnbase` - an R/Bioconductor package for
> isobaric tagged mass spectrometry data visualization, processing and
> quantitation*. Bioinformatics. 2012 Jan
> 15;28(2):288-9. doi:10.1093/bioinformatics/btr645. PubMed
> [PMID:22113085](http://www.ncbi.nlm.nih.gov/pubmed/22113085).

## Contributing

Contributions to the package are more than welcome. If you want to
contribute to this package, you should follow the same conventions as
the rest of the functions. Please do get in touch (preferable opening
a [github issue](https://github.com/lgatto/MSnbase/issues/)) to
discuss any suggestions. The
[`MSnbase` development vignette](http://lgatto.github.io/MSnbase/articles/MSnbase-development.html) gives
some background on the class infrastructure.

Please note that this project is released with a
[Contributor Code of Conduct](https://github.com/lgatto/MSnbase/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

