# Overview of missing value

Visualise missing values as a heatmap and barplots along the samples and
features.

## Usage

``` r
naplot(
  object,
  verbose = isMSnbaseVerbose(),
  reorderRows = TRUE,
  reorderColumns = TRUE,
  ...
)
```

## Arguments

- object:

  An object of class `MSnSet`.

- verbose:

  If verbose (default is
  [`isMSnbaseVerbose()`](https://lgatto.github.io/MSnbase/reference/MSnbaseOptions.md)),
  print a table of missing values.

- reorderRows:

  If reorderRows (default is `TRUE`) rows are ordered by number of NA.

- reorderColumns:

  If reorderColumns (default is `TRUE`) columns are ordered by number of
  NA.

- ...:

  Additional parameters passed to `image2`.

## Value

Used for its side effect. Invisibly returns `NULL`

## Author

Laurent Gatto

## Examples

``` r
data(naset)
naplot(naset)
```
