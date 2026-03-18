# Create a data with missing values

These functions take an instance of class
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
and sets randomly selected values to `NA`.

## Usage

``` r
makeNaData(object, nNA, pNA, exclude)

makeNaData2(object, nRows, nNAs, exclude)

whichNA(x)
```

## Arguments

- object:

  An instance of class `MSnSet`.

- nNA:

  The absolute number of missing values to be assigned.

- pNA:

  The proportion of missing values to be assignmed.

- exclude:

  A `vector` to be used to subset `object`, defining rows that should
  not be used to set `NA`s.

- nRows:

  The number of rows for each set.

- nNAs:

  The number of missing values for each set.

- x:

  A `matrix` or an instance of class `MSnSet`.

## Value

An instance of class `MSnSet`, as `object`, but with the appropriate
number/proportion of missing values. The returned object has an
additional feature meta-data columns, `nNA`

## Details

`makeNaData` randomly selects a number `nNA` (or a proportion `pNA`) of
cells in the expression matrix to be set to `NA`.

`makeNaData2` will select `length(nRows)` sets of rows from `object`,
each with `nRows[i]` rows respectively. The first set will be assigned
`nNAs[1]` missing values, the second `nNAs[2]`, ... As opposed to
`makeNaData`, this permits to control the number of `NAs` per rows.

The `whichNA` can be used to extract the indices of the missing values,
as illustrated in the example.

## Author

Laurent Gatto

## Examples

``` r
## Example 1
library(pRolocdata)
data(dunkley2006)
sum(is.na(dunkley2006))
#> [1] 0
dunkleyNA <- makeNaData(dunkley2006, nNA = 150)
processingData(dunkleyNA)
#> - - - Processing information - - -
#> Loaded on Thu Jul 16 22:53:08 2015. 
#> Normalised to sum of intensities. 
#> Added markers from  'mrk' marker vector. Thu Jul 16 22:53:08 2015 
#> Set 150 values to NA Wed Mar 18 17:49:18 2026 
#>  MSnbase version: 1.17.12 
sum(is.na(dunkleyNA))
#> [1] 150
table(fData(dunkleyNA)$nNA)
#> 
#>   0   1   2   3 
#> 558 113  17   1 
naIdx <- whichNA(dunkleyNA)
head(naIdx)
#>      [,1] [,2]
#> [1,]    4    1
#> [2,]   57    1
#> [3,]  110    1
#> [4,]  119    1
#> [5,]  147    1
#> [6,]  273    1
## Example 2
dunkleyNA <- makeNaData(dunkley2006, nNA = 150, exclude = 1:10)
processingData(dunkleyNA)
#> - - - Processing information - - -
#> Set 150 values to NA Wed Mar 18 17:49:18 2026
#>   (excluding 10 features) 
#>  MSnbase version: 1.17.12 
table(fData(dunkleyNA)$nNA[1:10])
#> 
#>  0 
#> 10 
table(fData(dunkleyNA)$nNA)
#> 
#>   0   1   2 
#> 549 130  10 
## Example 3
nr <- rep(10, 5)
na <- 1:5
x <- makeNaData2(dunkley2006[1:100, 1:5],
                 nRows = nr,
                 nNAs = na)
processingData(x)
#> - - - Processing information - - -
#> Loaded on Thu Jul 16 22:53:08 2015. 
#> Normalised to sum of intensities. 
#> Added markers from  'mrk' marker vector. Thu Jul 16 22:53:08 2015 
#> Subset [689,16][100,5] Wed Mar 18 17:49:18 2026 
#> Set (1,2,3,4,5) NAs in (10,10,10,10,10) rows,
#>   respectively Wed Mar 18 17:49:18 2026 
#>  MSnbase version: 1.17.12 
(res <- table(fData(x)$nNA))
#> 
#>  0  1  2  3  4  5 
#> 50 10 10 10 10 10 
stopifnot(as.numeric(names(res)[-1]) ==  na)
stopifnot(res[-1] ==  nr)
## Example 3
nr2 <- c(5, 12, 11, 8)
na2 <- c(3, 8, 1, 4)
x2 <- makeNaData2(dunkley2006[1:100, 1:10],
                  nRows = nr2,
                  nNAs = na2)
processingData(x2)
#> - - - Processing information - - -
#> Loaded on Thu Jul 16 22:53:08 2015. 
#> Normalised to sum of intensities. 
#> Added markers from  'mrk' marker vector. Thu Jul 16 22:53:08 2015 
#> Subset [689,16][100,10] Wed Mar 18 17:49:18 2026 
#> Set (3,8,1,4) NAs in (5,12,11,8) rows,
#>   respectively Wed Mar 18 17:49:18 2026 
#>  MSnbase version: 1.17.12 
(res2 <- table(fData(x2)$nNA))
#> 
#>  0  1  3  4  8 
#> 64 11  5  8 12 
stopifnot(as.numeric(names(res2)[-1]) ==  sort(na2))
stopifnot(res2[-1] ==  nr2[order(na2)])
## Example 5
nr3 <- c(5, 12, 11, 8)
na3 <- c(3, 8, 1, 3)
x3 <- makeNaData2(dunkley2006[1:100, 1:10],
                  nRows = nr3,
                  nNAs = na3)
processingData(x3)
#> - - - Processing information - - -
#> Loaded on Thu Jul 16 22:53:08 2015. 
#> Normalised to sum of intensities. 
#> Added markers from  'mrk' marker vector. Thu Jul 16 22:53:08 2015 
#> Subset [689,16][100,10] Wed Mar 18 17:49:18 2026 
#> Set (3,8,1,3) NAs in (5,12,11,8) rows,
#>   respectively Wed Mar 18 17:49:18 2026 
#>  MSnbase version: 1.17.12 
(res3 <- table(fData(x3)$nNA))
#> 
#>  0  1  3  8 
#> 64 11 13 12 
```
