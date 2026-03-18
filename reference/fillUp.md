# Fills up a vector

This function replaces all the empty characters `""` and/or `NA`s with
the value of the closest preceding the preceding non-`NA`/`""` element.
The function is used to populate dataframe or matrice columns where only
the cells of the first row in a set of partially identical rows are
explicitly populated and the following are empty.

## Usage

``` r
fillUp(x)
```

## Arguments

- x:

  a vector.

## Value

A vector as `x` with all empty characters `""` and `NA` values replaced
by the preceding non-`NA`/`""` value.

## Author

Laurent Gatto

## Examples

``` r
d <- data.frame(protein=c("Prot1","","","Prot2","",""),
                peptide=c("pep11","","pep12","pep21","pep22",""),
                score=c(1:2,NA,1:3))
d
#>   protein peptide score
#> 1   Prot1   pep11     1
#> 2                     2
#> 3           pep12    NA
#> 4   Prot2   pep21     1
#> 5           pep22     2
#> 6                     3
e <- apply(d,2,fillUp)
e
#>      protein peptide score
#> [1,] "Prot1" "pep11" " 1" 
#> [2,] "Prot1" "pep11" " 2" 
#> [3,] "Prot1" "pep12" " 2" 
#> [4,] "Prot2" "pep21" " 1" 
#> [5,] "Prot2" "pep22" " 2" 
#> [6,] "Prot2" "pep22" " 3" 
data.frame(e)
#>   protein peptide score
#> 1   Prot1   pep11     1
#> 2   Prot1   pep11     2
#> 3   Prot1   pep12     2
#> 4   Prot2   pep21     1
#> 5   Prot2   pep22     2
#> 6   Prot2   pep22     3
fillUp(d[,1])
#> [1] "Prot1" "Prot1" "Prot1" "Prot2" "Prot2" "Prot2"
```
