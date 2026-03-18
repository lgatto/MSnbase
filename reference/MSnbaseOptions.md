# MSnbase options

MSnbase defined a few options globally using the standard R options
mechanism. The current values of these options can be queried with
`MSnbaseOptions`. The options are:

- `verbose`: defines a session-wide verbosity flag, that is used if the
  `verbose` argument in individual functions is not set.

- `PARALLEL_THRESH`: defines the minimum number of spectra per file
  necessary before using parallel processing.

- `fastLoad`: `logical(1)`. If `TRUE` performs faster data loading for
  all methods of
  [OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
  that load data from the original files (such as
  [`spectrapply()`](https://lgatto.github.io/MSnbase/reference/pSet-class.md)).
  Users experiencing data I/O errors (observed mostly on macOS systems)
  should set this option to `FALSE`.

## Usage

``` r
MSnbaseOptions()

isMSnbaseVerbose()

setMSnbaseVerbose(opt)

setMSnbaseParallelThresh(opt = 1000)

setMSnbaseFastLoad(opt = TRUE)

isMSnbaseFastLoad()
```

## Arguments

- opt:

  The value of the new option

## Value

A `list` of MSnbase options and the single option values for the
individual accessors.

## Details

`isMSnbaseVerbose` is one wrapper for the verbosity flag, also available
through `options("MSnbase")$verbose`.

There are also setters to set options individually. When run without
argument, the verbosity setter inverts the current value of the option.
