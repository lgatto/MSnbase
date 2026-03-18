# Return a variable name

Return the name of variable `varname` in call `match_call`.

## Usage

``` r
getVariableName(match_call, varname)
```

## Arguments

- match_call:

  An object of class `call`, as returned by `match.call`.

- varname:

  An `character` of length 1 which is looked up in `match_call`.

## Value

A `character` with the name of the variable passed as parameter
`varname` in parent close of `match_call`.

## Author

Laurent Gatto

## Examples

``` r
a <- 1
f <- function(x, y)
 MSnbase:::getVariableName(match.call(), "x")
f(x = a)
#> [1] "a"
f(y = a)
#> character(0)
```
