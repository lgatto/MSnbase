# Read an 'mzTab' file

This function can be used to create a
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
by reading and parsing an `mzTab` file. The metadata section is always
used to populate the `MSnSet`'s `experimentData` slot.

## Usage

``` r
readMzTabData_v0.9(file, what = c("PRT", "PEP"), verbose = isMSnbaseVerbose())
```

## Arguments

- file:

  A `character` with the `mzTab` file to be read in.

- what:

  One of `"PRT"` or `"PEP"`, defining which of protein of peptide
  section should be parse. The metadata section, when available, is
  always used to populate the `experimentData` slot.

- verbose:

  Produce verbose output.

## Value

An instance of class `MSnSet`.

## See also

[`writeMzTabData`](https://lgatto.github.io/MSnbase/reference/writeMzTabData.md)
to save an
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
as an `mzTab` file.

## Author

Laurent Gatto
