# Returns the matching column names of indices.

Given a text spread sheet `f` and a `pattern` to be matched to its
header (first line in the file), the function returns the matching
columns names or indices of the corresponding `data.frame`.

The function starts by reading the first line of the file (or
connection) `f` with
[`readLines`](https://rdrr.io/r/base/readLines.html), then splits it
according to the optional `...` arguments (it is important to correctly
specify [`strsplit`](https://rdrr.io/r/base/strsplit.html)'s `split`
character vector here) and then matches `pattern` to the individual
column names using
[`grep`](https://rdrr.io/pkg/BiocGenerics/man/grep.html).

Similarly, `getEcols` can be used to explore the column names and decide
for the appropriate `pattern` value.

These functions are useful to check the parameters to be provided to
[`readMSnSet2`](https://lgatto.github.io/MSnbase/reference/readMSnSet.md).

## Usage

``` r
grepEcols(f, pattern, ..., n = 1)

getEcols(f, ..., n = 1)
```

## Arguments

- f:

  A connection object or a `character` string to be read in with
  `readLines(f, n = 1)`.

- pattern:

  A `character` string containing a regular expression to be matched to
  the file's header.

- ...:

  Additional parameters passed to
  [`strsplit`](https://rdrr.io/r/base/strsplit.html) to split the file
  header into individual column names.

- n:

  An `integer` specifying which line in file `f` to grep (get). Default
  is 1. Note that this argument must be named.

## Value

Depending on `value`, the matching column names of indices. In case of
`getEcols`, a `character` of column names.

## See also

[`readMSnSet2`](https://lgatto.github.io/MSnbase/reference/readMSnSet.md)

## Author

Laurent Gatto
