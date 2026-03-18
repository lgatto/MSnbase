# Number of precursor selection events

`precSelection` computes the number of selection events each precursor
ions has undergone in an tandem MS experiment. This will be a function
of amount of peptide loaded, chromatography efficiency, exclusion
time,... and is useful when optimising and experimental setup. This
function returns a named integer vector or length equal to the number of
unique precursor MZ values in the original experiment. See `n` parameter
to set the number of MZ significant decimals.

`precSelectionTable` is a wrapper around `precSelection` and returns a
table with the number of single, 2-fold, ... selection events.

## Usage

``` r
precSelection(object,n)
```

## Arguments

- object:

  An instane of class
  `"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`.

- n:

  The number of decimal places to round the precursor MZ to. Is passed
  to the [round](https://rdrr.io/r/base/Round.html) function.

## Value

A named integer in case of `precSelection` and a `table` for
`precSelectionTable`.

## Author

Laurent Gatto

## Examples

``` r
precSelection(itraqdata)
#> 520.78326416 573.95385742 401.73919678 567.83392334 488.32687378   782.871521 
#>            1            1            1            1            1            1 
#> 773.35400391 671.33288574 650.87664795 667.37115479 648.34741211 459.75946045 
#>            1            1            1            1            1            1 
#> 495.81802368 641.32049561 596.34204102 521.60577393 529.61083984 509.32354736 
#>            1            1            1            1            1            1 
#> 662.84246826 551.85644531 1097.5622559 974.49163818 645.37414551 662.69012451 
#>            1            1            1            1            1            1 
#> 776.90661621 649.99621582  502.6104126 518.27526855 506.80905151 1078.0333252 
#>            1            1            1            1            1            1 
#> 686.73516846 541.33233643  824.4263916 546.95861816 744.04071045 599.31079102 
#>            1            1            1            1            1            1 
#> 803.94683838 646.67144775  1236.133667 1115.5579834 682.05841064 1022.5847168 
#>            1            1            1            1            1            1 
#> 585.33795166  580.8616333 819.93371582 877.50390625   529.347229 651.91625977 
#>            1            1            1            1            1            1 
#>  434.9473877 768.43304443 472.28570557 716.34051514 437.80401611 525.29980469 
#>            1            1            1            1            1            1 
#>  538.8326416 
#>            1 
precSelection(itraqdata,n=2)
#>  520.78  573.95  401.74  567.83  488.33  782.87  773.35  671.33  650.88  667.37 
#>       1       1       1       1       1       1       1       1       1       1 
#>  648.35  459.76  495.82  641.32  596.34  521.61  529.61  509.32  662.84  551.86 
#>       1       1       1       1       1       1       1       1       1       1 
#> 1097.56  974.49  645.37  662.69  776.91     650  502.61  518.28  506.81 1078.03 
#>       1       1       1       1       1       1       1       1       1       1 
#>  686.74  541.33  824.43  546.96  744.04  599.31  803.95  646.67 1236.13 1115.56 
#>       1       1       1       1       1       1       1       1       1       1 
#>  682.06 1022.58  585.34  580.86  819.93   877.5  529.35  651.92  434.95  768.43 
#>       1       1       1       1       1       1       1       1       1       1 
#>  472.29  716.34   437.8   525.3  538.83 
#>       1       1       1       1       1 
precSelectionTable(itraqdata)
#> x
#>  1 
#> 55 
## only single selection event in this reduced exeriment
```
