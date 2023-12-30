[![codecov.io](https://codecov.io/github/lgatto/MSnbase/coverage.svg?branch=master)](https://codecov.io/github/lgatto/MSnbase?branch=master)
[![R-CMD-check-bioc](https://github.com/lgatto/MSnbase/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lgatto/MSnbase/actions?query=workflow%3AR-CMD-check-bioc)


# The `MSnbase` package

<img align = "right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/MSnbase/MSnbase.png" height="200">

[**MSnbase**](https://lgatto.github.io/MSnbase/) is an R/Bioconductor
package that provides infrastructure for plotting, manipulation and
processing mass spectrometry and proteomics data. The project was
started by [Laurent Gatto](https://lgatto.github.io/) in October 2010
(Mon Oct 4 23:35:23 2010, according to the git log) and has, since
then, benefited from
[various contributions](https://lgatto.github.io/msnbase-contribs/), in
particular [Sebastian Gibb](https://sebastiangibb.de/)
and [Johannes Rainer](https://github.com/jorainer).

The official package page is the Bioconductor landing page
([release](https://www.bioconductor.org/packages/release/bioc/html/MSnbase.html) or
[devel](https://www.bioconductor.org/packages/devel/bioc/html/MSnbase.html) versions). The
[github page](https://github.com/lgatto/MSnbase) page is for active
development, issue tracking and forking/pulling purposes.

To get an overview of the package, see the
[*MSnbase-demo*](https://lgatto.github.io/MSnbase/articles/v01-MSnbase-demo.html)
vignette. More vignettes are available in the *Articles* tab.

## The R for Mass Spectrometry initiative

The aim of the *R for Mass Spectrometry* initiative is to provide
efficient, thoroughly documented, tested and flexible R software for
the analysis and interpretation of high throughput mass spectrometry
assays, including proteomics and metabolomics experiments. The project
formalises the longtime collaborative development efforts of its core
members under the R for Mass Spectrometry organisation to facilitate
dissemination and accessibility of their work.

If you are using MSnbase, consider switching to the **R for Mass
Spectrometry** packages, in particular,
[Spectra](https://rformassspectrometry.github.io/Spectra/) for raw
data, [PSMatch](https://rformassspectrometry.github.io/PSMatch/) for
identification data, and
[QFeatures](https://rformassspectrometry.github.io/QFeatures/) for
quantitative data. See https://RforMassSpectrometry.org for details.

## Installation

To install the package:


```r
install.packages("BiocManager")
BiocManager::install("MSnbase")
```

If you need the github version (not recommended unless you know what
you are doing), use


```r
BiocManager::install("lgatto/MSnbase")
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
> [PMID:22113085](https://www.ncbi.nlm.nih.gov/pubmed/22113085).


> *`MSnbase`, efficient and elegant R-based processing and
> visualisation of raw mass spectrometry data*. Laurent Gatto,
> Sebastian Gibb, Johannes Rainer. bioRxiv 2020.04.29.067868; doi:
> https://doi.org/10.1101/2020.04.29.067868


## Contributing

Contributions to the package are more than welcome. If you want to
contribute to this package, you should follow the same conventions as
the rest of the functions. Please do get in touch (preferable opening
a [github issue](https://github.com/lgatto/MSnbase/issues/)) to
discuss any suggestions. The
[`MSnbase` development vignette](https://lgatto.github.io/MSnbase/articles/v05-MSnbase-development.html) gives
some background on the class infrastructure.

Please note that this project is released with a
[Contributor Code of Conduct](https://github.com/lgatto/MSnbase/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.
