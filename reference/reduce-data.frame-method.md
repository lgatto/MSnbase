# Reduce a data.frame

Reduce a data.frame so that the (primary) key column contains only
unique entries and other columns pertaining to that entry are combined
into semicolon-separated values into a single row/observation.

An important side-effect of reducing a `data.frame` is that all columns
other than the key are converted to characters when they are collapsed
to a semi-column separated value (even if only one value is present) as
soon as one observation of transformed.

## Usage

``` r
# S4 method for class 'data.frame'
reduce(x, key, sep = ";")
```

## Arguments

- x:

  A `data.frame`.

- key:

  The column name (currenly only one is supported) to be used as primary
  key.

- sep:

  The separator. Default is `;`.

## Value

A reduced `data.frame`.

## Author

Laurent Gatto

## Examples

``` r
dfr <- data.frame(A = c(1, 1, 2),
                  B = c("x", "x", "z"),
                  C = LETTERS[1:3])
dfr
#>   A B C
#> 1 1 x A
#> 2 1 x B
#> 3 2 z C
dfr2 <- reduce(dfr, key = "A")
dfr2
#>   A   B   C
#> 1 1 x;x A;B
#> 2 2   z   C
## column A used as key is still num
str(dfr2)
#> 'data.frame':    2 obs. of  3 variables:
#>  $ A: num  1 2
#>  $ B: chr  "x;x" "z"
#>  $ C: chr  "A;B" "C"
dfr3 <- reduce(dfr, key = "B")
dfr3
#>     A B   C
#> 1 1;1 x A;B
#> 2   2 z   C
## A is converted to chr; B remains factor
str(dfr3)
#> 'data.frame':    2 obs. of  3 variables:
#>  $ A: chr  "1;1" "2"
#>  $ B: chr  "x" "z"
#>  $ C: chr  "A;B" "C"
dfr4 <- data.frame(A = 1:3,
                   B = LETTERS[1:3],
                   C = c(TRUE, FALSE, NA))
## No effect of reducing, column classes are maintained
str(reduce(dfr4, key = "B"))
#> 'data.frame':    3 obs. of  3 variables:
#>  $ A: int  1 2 3
#>  $ B: chr  "A" "B" "C"
#>  $ C: logi  TRUE FALSE NA
```
