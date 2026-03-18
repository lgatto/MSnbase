# Tests equality of list elements class

Compares equality of all members of a list.

## Usage

``` r
listOf(x, class, valid = TRUE)
```

## Arguments

- x:

  A `list`.

- class:

  A `character` defining the expected class.

- valid:

  A `logical` defining if all elements should be tested for validity.
  Default is `TRUE`.

## Value

`TRUE` is all elements of `x` inherit from `class`.

## Author

Laurent Gatto

## Examples

``` r
listOf(list(), "foo")
#> [1] TRUE
listOf(list("a", "b"), "character")
#> [1] TRUE
listOf(list("a", 1), "character")
#> [1] FALSE
```
