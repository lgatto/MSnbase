# Non-parametric coefficient of variation

Calculates a non-parametric version of the coefficient of variation
where the standard deviation is replaced by the median absolute
deviations (see `mad` for details) and divided by the absolute value of
the mean.

Note that the `mad` of a single value is 0 (as opposed to `NA` for the
standard deviation, see example below).

## Usage

``` r
npcv(x, na.rm = TRUE)
```

## Arguments

- x:

  A `numeric`.

- na.rm:

  A `logical` (default is `TRUE` indicating whether `NA` values should
  be stripped before the computation of the median absolute deviation
  and mean.

## Value

A `numeric`.

## Author

Laurent Gatto

## Examples

``` r
set.seed(1)
npcv(rnorm(10))
#> [1] 5.852407
replicate(10, npcv(rnorm(10)))
#>  [1]  3.112625  6.088948  4.928480  5.534139  6.124995  2.354319  2.996464
#>  [8]  6.049891 11.716089  3.265193
npcv(1)
#> [1] 0
mad(1)
#> [1] 0
sd(1)
#> [1] NA
```
