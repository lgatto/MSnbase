# MSnbase benchmarking

## Introduction

In this vignette, we will document various timings and benchmarkings of
the *[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* version
2, that focuses on *on-disk* data access (as opposed to *in-memory*).
More details about the new implementation are documented in the
respective classes manual pages and in

> *`MSnbase`, efficient and elegant R-based processing and visualisation
> of raw mass spectrometry data*. Laurent Gatto, Sebastian Gibb,
> Johannes Rainer. bioRxiv 2020.04.29.067868; doi:
> <https://doi.org/10.1101/2020.04.29.067868>

As a benchmarking dataset, we are going to use a subset of an TMT 6-plex
experiment acquired on an LTQ Orbitrap Velos, that is distributed with
the *[MsDataHub](https://bioconductor.org/packages/3.23/MsDataHub)*
package

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`"MsDataHub"`](https://rformassspectrometry.github.io/MsDataHub)`)`\
`f`` ``<-`` `[`TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz`](https://rformassspectrometry.github.io/MsDataHub/reference/PXD000001.html)`(``)`

    ## see ?MsDataHub and browseVignettes('MsDataHub') for documentation

    ## loading from cache

We need to load the
*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* package and
set the session-wide verbosity flag to `FALSE`.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`"MSnbase"`](https://lgatto.github.io/MSnbase)`)`\
[`setMSnbaseVerbose`](https://lgatto.github.io/MSnbase/reference/MSnbaseOptions.md)`(``FALSE``)`

## Benchmarking

### Reading data

We first read the data using the original behaviour `readMSData`
function by setting the `mode` argument to `"inMemory"` to generates an
in-memory representation of the MS2-level raw data and measure the time
needed for this operation.

\
[`system.time`](https://rdrr.io/r/base/system.time.html)`(``inmem`` ``<-`` `[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)`(``f``, msLevel. ``=`` ``2``,`\
`                                mode ``=`` ``"inMemory"``,`\
`                                centroided. ``=`` ``TRUE``)``)`

    ##    user  system elapsed 
    ##  41.062   0.371  41.262

Next, we use the `readMSData` function to generate an on-disk
representation of the same data by setting `mode = "onDisk"`.

\
[`system.time`](https://rdrr.io/r/base/system.time.html)`(``ondisk`` ``<-`` `[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.md)`(``f``, msLevel. ``=`` ``2``,`\
`                                  mode ``=`` ``"onDisk"``,`\
`                                  centroided. ``=`` ``TRUE``)``)`

    ##    user  system elapsed 
    ##   9.792   0.205   9.834

Creating the on-disk experiment is considerable faster and scales to
much bigger, multi-file data, both in terms of object creation time, but
also in terms of object size (see next section). We must of course make
sure that these two datasets are equivalent:

\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``inmem``, ``ondisk``)`

    ## [1] TRUE

### Data size

To compare the size occupied in memory of these two objects, we are
going to use the `object.size` function, which accounts for the data
(the spectra) in the `assayData` environment (as opposed to the
`object.size` function from the `utils` package).

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`object.size`](https://rdrr.io/r/utils/object.size.html)`(``inmem``)``, units ``=`` ``"MiB"``)`

    ## 0.5 MiB

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`object.size`](https://rdrr.io/r/utils/object.size.html)`(``ondisk``)``, units ``=`` ``"MiB"``)`

    ## 2.8 MiB

The difference is explained by the fact that for `ondisk`, the spectra
are not created and stored in memory; they are access on disk when
needed, such as for example for plotting:

\
[`plot`](https://lgatto.github.io/MSnbase/reference/plot-methods.md)`(``inmem``[[``200``]``]``, full ``=`` ``TRUE``)`\
[`plot`](https://lgatto.github.io/MSnbase/reference/plot-methods.md)`(``ondisk``[[``200``]``]``, full ``=`` ``TRUE``)`

![Plotting in-memory and on-disk
spectra](v04-benchmarking_files/figure-html/plot1-1.png)

Plotting in-memory and on-disk spectra

### Accessing spectra

The drawback of the on-disk representation is when the spectrum data has
to actually be accessed. To compare access time, we are going to use the
*[microbenchmark](https://CRAN.R-project.org/package=microbenchmark)*
and repeat access 10 times to compare access to all 6103 and a single
spectrum in-memory (i.e. pre-loaded and constructed) and on-disk
(i.e. on-the-fly access).

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`"microbenchmark"`](https://github.com/joshuaulrich/microbenchmark/)`)`\
`mb`` ``<-`` `[`microbenchmark`](https://rdrr.io/pkg/microbenchmark/man/microbenchmark.html)`(`[`spectra`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`(``inmem``)``,`\
`                     ``inmem``[[``200``]``]``,`\
`                     `[`spectra`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`(``ondisk``)``,`\
`                     ``ondisk``[[``200``]``]``,`\
`                     times ``=`` ``10``)`\
`mb`

    ## Unit: microseconds
    ##             expr         min          lq         mean      median          uq
    ##   spectra(inmem)     922.239    1434.270    1764.1386    1916.486    2024.046
    ##     inmem[[200]]      20.170      22.363      62.5161      57.170      81.742
    ##  spectra(ondisk) 4004717.471 4020451.792 4734164.1678 4037820.935 5795670.019
    ##    ondisk[[200]] 1597398.145 1609258.647 1619363.0977 1619704.890 1628712.366
    ##          max neval
    ##     2371.571    10
    ##      134.180    10
    ##  5816314.962    10
    ##  1640935.447    10

While it takes order or magnitudes more time to access the data
on-the-fly rather than a pre-generated spectrum, accessing all spectra
is only marginally slower than accessing all spectra, as most of the
time is spent preparing the file for access, which is done only once.

On-disk access performance will depend on the read throughput of the
disk. A comparison of the data import of the above file from an internal
solid state drive and from an USB3 connected hard disk showed only small
differences for the `onDisk` mode (1.07 *vs* 1.36 seconds), while no
difference were observed for accessing individual or all spectra. Thus,
for this particular setup, performance was about the same for SSD and
HDD. This might however not apply to setting in which data import is
performed in parallel from multiple files.

Data access does not prohibit interactive usage, such as plotting, for
example, as it is about 1/2 seconds, which is an operation that is
relatively rare, compared to subsetting and filtering, which are faster
for on-disk data:

\
`i`` ``<-`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`length`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)`(``inmem``)``, ``100``)`\
[`system.time`](https://rdrr.io/r/base/system.time.html)`(``inmem``[``i``]``)`

    ##    user  system elapsed 
    ##   0.133   0.001   0.134

\
[`system.time`](https://rdrr.io/r/base/system.time.html)`(``ondisk``[``i``]``)`

    ##    user  system elapsed 
    ##   0.011   0.000   0.010

Operations on the spectra data, such as peak picking, smoothing,
cleaning, … are cleverly cached and only applied when the data is
accessed, to minimise file access overhead. Finally, specific operations
such as for example quantitation (see next section) are optimised for
speed.

### MS2 quantitation

Below, we perform TMT 6-plex reporter ions quantitation on the first 100
spectra and verify that the results are identical (ignoring feature
names).

\
[`system.time`](https://rdrr.io/r/base/system.time.html)`(``eim`` ``<-`` `[`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)`(``inmem``[``1``:``100``]``, reporters ``=`` ``TMT6``,`\
`                            method ``=`` ``"max"``)``)`

    ##    user  system elapsed 
    ##   2.114   1.223   1.428

\
[`system.time`](https://rdrr.io/r/base/system.time.html)`(``eod`` ``<-`` `[`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)`(``ondisk``[``1``:``100``]``, reporters ``=`` ``TMT6``,`\
`                            method ``=`` ``"max"``)``)`

    ##    user  system elapsed 
    ##   1.677   0.267   1.812

\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``eim``, ``eod``, check.attributes ``=`` ``FALSE``)`

    ## [1] TRUE

## Notable differences *on-disk* and *in-memory* implementations

The `MSnExp` and `OnDiskMSnExp` documentation files and the *MSnbase
developement* vignette provide more information about implementation
details.

### MS levels

*On-disk* support multiple MS levels in one object, while *in-memory*
only supports a single level. While support for multiple MS levels could
be added to the in-memory back-end, memory constrains make this
pretty-much useless and will most likely never happen.

### Serialisation

*In-memory* objects can be
[`save()`](https://rdrr.io/r/base/save.html)ed and
[`load()`](https://rdrr.io/r/base/load.html)ed, while *on-disk* can’t.
As a workaround, the latter can be coerced to *in-memory* instances with
`as(, "MSnExp")`. We would need `mzML` write support in
*[mzR](https://bioconductor.org/packages/3.23/mzR)* to be able to
implement serialisation for *on-disk* data.

### Data processing

Whenever possible, accessing and processing *on-disk* data is delayed
(*lazy* processing). These operations are stored in a *processing queue*
until the spectra are effectively instantiated.

### Validity

The *on-disk* `validObject` method doesn’t verify the validity on the
spectra (as there aren’t any to check). The `validateOnDiskMSnExp`
function, on the other hand, instantiates all spectra and checks their
validity (in addition to calling `validObject`).

## Conclusions

This document focuses on speed and size improvements of the new on-disk
`MSnExp` representation. The extend of these improvements will
substantially increase for larger data.

For general functionality about the on-disk `MSnExp` data class and
*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* in general,
see other vignettes available with

\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``package ``=`` ``"MSnbase"``)`
