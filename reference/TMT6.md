# TMT 6/10-plex sets

This instance of class
`"`[`ReporterIons`](https://lgatto.github.io/MSnbase/reference/ReporterIons-class.md)`"`
corresponds to the TMT 6-plex set, i.e the 126, 127, 128, 129, 130 and
131 isobaric tags. In the `TMT7` data set, an unfragmented tag, i.e
reporter and attached isobaric tag, is also included at MZ 229. A second
`TMT6b` has slightly different values.

The `TMT10` instance corresponds to the 10-plex version. There are
spectific HCD (`TMT10HCD`, same as `TMT10`) and ETD (`TMT10ETD`) sets.

These objects are used to plot the reporter ions of interest in an MSMS
spectra (see
`"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`)
as well as for quantification (see
[`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)).

## Usage

``` r
TMT6
TMT6b
TMT7
TMT7b
TMT10
TMT10ETD
TMT10HCD
TMT11
TMT11HCD
```

## References

Thompson A, Schäfer J, Kuhn K, Kienle S, Schwarz J, Schmidt G, Neumann
T, Johnstone R, Mohammed AK, Hamon C. "Tandem mass tags: a novel
quantification strategy for comparative analysis of complex protein
mixtures by MS/MS." *Anal Chem.* 2003 Apr 15;75(8):1895-904. *Erratum*
in: *Anal Chem.* 2006 Jun 15;78(12):4235. Mohammed, A Karim A \[added\]
and *Anal Chem.* 2003 Sep 15;75(18):4942. Johnstone, R \[added\]. PubMed
PMID: 12713048.

## See also

[`iTRAQ4`](https://lgatto.github.io/MSnbase/reference/iTRAQ4.md).

## Examples

``` r
TMT6
#> Object of class "ReporterIons"
#> TMT6: '6-plex TMT tags' with 6 reporter ions
#>  - [TMT6.126] 126.1277 +/- 0.05 (red)
#>  - [TMT6.127] 127.1311 +/- 0.05 (purple)
#>  - [TMT6.128] 128.1344 +/- 0.05 (blue)
#>  - [TMT6.129] 129.1378 +/- 0.05 (steelblue)
#>  - [TMT6.130] 130.1411 +/- 0.05 (green)
#>  - [TMT6.131] 131.1382 +/- 0.05 (yellow)
TMT6[1:2]
#> Object of class "ReporterIons"
#> TMT6[1:2]: 'subset of 6-plex TMT tags' with 2 reporter ions
#>  - [TMT6.126] 126.1277 +/- 0.05 (red)
#>  - [TMT6.127] 127.1311 +/- 0.05 (purple)

TMT10
#> Object of class "ReporterIons"
#> TMT10HCD: '10-plex TMT HCD' with 10 reporter ions
#>  - [126] 126.1277 +/- 0.002 (#8DD3C7)
#>  - [127N] 127.1248 +/- 0.002 (#FFFFB3)
#>  - [127C] 127.1311 +/- 0.002 (#BEBADA)
#>  - [128N] 128.1281 +/- 0.002 (#FB8072)
#>  - [128C] 128.1344 +/- 0.002 (#80B1D3)
#>  - [129N] 129.1315 +/- 0.002 (#FDB462)
#>  - [129C] 129.1378 +/- 0.002 (#B3DE69)
#>  - [130N] 130.1348 +/- 0.002 (#FCCDE5)
#>  - [130C] 130.1411 +/- 0.002 (#D9D9D9)
#>  - [131] 131.1382 +/- 0.002 (#BC80BD)

newReporter <- new("ReporterIons",
                   description="an example",
                   name="my reporter ions",
                   reporterNames=c("myrep1","myrep2"),
                   mz=c(121,122),
                   col=c("red","blue"),
                   width=0.05)
newReporter
#> Object of class "ReporterIons"
#> my reporter ions: 'an example' with 2 reporter ions
#>  - [myrep1] 121 +/- 0.05 (red)
#>  - [myrep2] 122 +/- 0.05 (blue)
```
