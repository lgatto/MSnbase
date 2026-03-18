# Export an MzTab object as mzTab file.

`writeMzTabData` exports an
[MzTab](https://lgatto.github.io/MSnbase/reference/MzTab-class.md)
object as mzTab file. Note that the comment section "COM" are not
written out.

## Usage

``` r
writeMzTabData(
  object,
  file,
  what = c("MT", "PEP", "PRT", "PSM", "SML", "SMF", "SME")
)
```

## Arguments

- object:

  [MzTab](https://lgatto.github.io/MSnbase/reference/MzTab-class.md)
  object, either read in by
  [`MzTab()`](https://lgatto.github.io/MSnbase/reference/MzTab-class.md)
  or assembled.

- file:

  `character(1)` with the file name.

- what:

  `character` with names of the sections to be written out. Expected
  sections are `"MT"`, `"PEP"`, `"PRT"`, `"PSM"`, `"SML"`, `"SMF"`, or
  `"SME"`.

## Author

Steffen Neumann
