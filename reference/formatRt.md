# Format Retention Time

This function is used to convert retention times. Conversion is seconds
to/from the more human friendly format "mm:sec". The implementation is
from
[`MsCoreUtils::formatRt()`](https://rdrr.io/pkg/MsCoreUtils/man/formatRt.html).

## Usage

``` r
formatRt(rt)
```

## Arguments

- rt:

  retention time in seconds (`numeric`) or "mm:sec" (`character`).

## Value

A vector of same length as `rt`.

## Author

Laurent Gatto and Sebastian Gibb

## Examples

``` r

formatRt(1524)
#> [1] "25:24"
formatRt("25:24")
#> [1] 1524
```
