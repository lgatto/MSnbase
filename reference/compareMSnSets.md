# Compare two MSnSets

Compares two
[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
instances. The `qual` and `processingData` slots are generally omitted.

## Usage

``` r
compareMSnSets(x, y, qual = FALSE, proc = FALSE)
```

## Arguments

- x:

  First MSnSet

- y:

  Second MSnSet

- qual:

  Should the `qual` slots be compared? Default is `FALSE`.

- proc:

  Should the `processingData` slots be compared? Default is `FALSE`.

## Value

A `logical`

## Author

Laurent Gatto
