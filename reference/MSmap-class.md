# Class `MSmap`

A class to store mass spectrometry data maps, i.e intensities collected
along the M/Z and retention time space during a mass spectrometry
acquisition.

## Objects from the Class

Objects can be created with the `MSmap` constructor. The constructor has
the following arguments:

- object:

  An object created by
  [`mzR::openMSfile`](https://rdrr.io/pkg/mzR/man/openMSfile.html) or an
  instance of class
  [`OnDiskMSnExp`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md).
  If the latter contains data from multiple files, a warning will be
  issued and the first one will be used.

- lowMz:

  A `numeric` of length 1 defining the lower bound of the M/Z range of
  the MS map.

- highMz:

  A `numeric` of length 1 defining the upper bound of the M/Z range of
  the MS map.

- resMz:

  The resolution along the M/Z range.

- hd:

  An optional `data.frame` as produced by `mzR::header(object)`. If
  missing, will be computer within the function. Ignored when `object`
  is an `OnDiskMSnExp`.

- zeroIsNA:

  Set 0 intensities to `NA`. This can be used to clarify the 3
  dimensional plot produce by `plot3D`.

## Slots

- `call`::

  Object of class `"call"` - the call used to generate the instance.

- `map`::

  Object of class `"matrix"` containing the actual MS map.

- `mz`::

  Object of class `"numeric"` with the M/Z sampling bins.

- `res`::

  Object of class `"numeric"` storing the the M/Z resolution used to
  create the map.

- `rt`::

  Object of class `"numeric"` with the retention times of the map
  spectra.

- `ms`::

  Object of class `"numeric"` with the MS levels of the spectra.

- `t`::

  Object of class `"logical"` indicating if the instance has been
  transposed.

- `filename`::

  Object of class `"character"` specifying the filename of the original
  raw MS data.

## Methods

- coerce:

  `signature(from = "MSmap", to = "data.frame")`: convert the `MSmap`
  instance in a `data.frame`. Useful for plotting with `lattice` or
  `ggplot2`.

- fileName:

  `signature(object = "MSmap")`: returns the raw data filename.

- msLevel:

  `signature(object = "MSmap")`: returns the MS level of the map
  spectra.

- msMap:

  `signature(object = "MSmap")`: returns the actual map `matrix`.

- mz:

  `signature(object = "MSmap", ...)`: returns the M/Z values of the map.
  Additional arguments are currently ignored.

- rtime:

  `signature(object = "MSmap", ...)`: returns retention time values of
  the map. Additional arguments are currently ignored.

- mzRes:

  `signature(object = "MSmap")`: returns the resolution with which the
  sample along the M/Z range was done.

- dim:

  `signature(x = "MSmap")`: returns the dimensions of the map. `ncol`
  and `nrow` return the number of columns and rows respectively.

- t:

  `signature(x = "MSmap")`: transposes the map.

- show:

  `signature(object = "MSmap")`: prints a summary of the map.

- plot:

  `signature(x = "MSmap", allTicks = "logical")`: produces an image of
  the map using
  [`lattice::levelplot`](https://rdrr.io/pkg/lattice/man/levelplot.html).
  By default, `allTicks` is `TRUE` and all M/Z and retention times ticks
  of drawn. If set to `FALSE`, only 10 ticks in each dimension are
  plotted.

- plot3D:

  `signature(object = "MSmap", rgl = "logical")`: produces an three
  dimensional view of the map using `lattice::cloude(..., type = "h")`.
  If `rgl` is `TRUE`, the map is visualised on a `rgl` device and can be
  rotated with the mouse.

## Author

Laurent Gatto

## Examples

``` r

if (FALSE) { # \dontrun{
    ## downloads the data
    library("rpx")
    px1 <- PXDataset("PXD000001")
    (i <- grep("TMT.+mzML", pxfiles(px1), value = TRUE))
    mzf <- pxget(px1, i)

    ## Using an mzRpwiz object
    ## reads the data
    ms <- openMSfile(mzf)
    hd <- header(ms)

    ## a set of spectra of interest: MS1 spectra eluted
    ## between 30 and 35 minutes retention time
    ms1 <- which(hd$msLevel == 1)
    rtsel <- hd$retentionTime[ms1] / 60 > 30 &
        hd$retentionTime[ms1] / 60 < 35

    ## the map
    M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd)

    plot(M, aspect = 1, allTicks = FALSE)
    plot3D(M)
    if (require("rgl") & interactive())
        plot3D(M, rgl = TRUE)

    ## With some MS2 spectra
    i <- ms1[which(rtsel)][1]
    j <- ms1[which(rtsel)][2]
    M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
    plot3D(M2)

    ## Using an OnDiskMSnExp object and accessors
    msn <- readMSData(mzf, mode = "onDisk")

    ## a set of spectra of interest: MS1 spectra eluted
    ## between 30 and 35 minutes retention time
    ms1 <- which(msLevel(msn) == 1)
    rtsel <- rtime(msn)[ms1] / 60 > 30 &
        rtime(msn)[ms1] / 60 < 35

    ## the map
    M3 <- MSmap(msn, ms1[rtsel], 521, 523, .005)
    plot(M3, aspect = 1, allTicks = FALSE)

    ## With some MS2 spectra
    i <- ms1[which(rtsel)][1]
    j <- ms1[which(rtsel)][2]
    M4 <- MSmap(msn, i:j, 100, 1000, 1)
    plot3D(M4)
} # }
```
