# Convert to camel case by replacing dots by captial letters

Convert a `vector` of characters to camel case by replacing dots by
captial letters.

## Usage

``` r
makeCamelCase(x, prefix)
```

## Arguments

- x:

  A `vector` to be transformed to camel case.

- prefix:

  An optional `character` of length one. Any additional elements are
  ignores.

## Value

A `character` of same length as `x`.

## Author

Laurent Gatto

## Examples

``` r
nms <- c("aa.foo", "ab.bar")
makeCamelCase(nms)
#> [1] "aaFoo" "abBar"
makeCamelCase(nms, prefix = "x")
#> [1] "xAaFoo" "xAbBar"
```
