# iTRAQ 4-plex set

This instance of class
`"`[`ReporterIons`](https://lgatto.github.io/MSnbase/reference/ReporterIons-class.md)`"`
corresponds to the iTRAQ 4-plex set, i.e the 114, 115, 116 and 117
isobaric tags. In the iTRAQ5 data set, an unfragmented tag, i.e reporter
and attached isobaric tag, is also included at MZ 145. These objects are
used to plot the reporter ions of interest in an MSMS spectra (see
`"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`)
as well as for quantification (see
[`quantify`](https://lgatto.github.io/MSnbase/reference/quantify-methods.md)).

## Usage

``` r
iTRAQ4
iTRAQ5
iTRAQ8
iTRAQ9
```

## References

Ross PL, Huang YN, Marchese JN, Williamson B, Parker K, Hattan S,
Khainovski N, Pillai S, Dey S, Daniels S, Purkayastha S, Juhasz P,
Martin S, Bartlet-Jones M, He F, Jacobson A, Pappin DJ. "Multiplexed
protein quantitation in Saccharomyces cerevisiae using amine-reactive
isobaric tagging reagents." *Mol Cell Proteomics*, 2004
Dec;3(12):1154-69. Epub 2004 Sep 22. PubMed PMID: 15385600.

## See also

[`TMT6`](https://lgatto.github.io/MSnbase/reference/TMT6.md).

## Examples

``` r
iTRAQ4
#> Object of class "ReporterIons"
#> iTRAQ4: '4-plex iTRAQ' with 4 reporter ions
#>  - [iTRAQ4.114] 114.1112 +/- 0.05 (red)
#>  - [iTRAQ4.115] 115.1083 +/- 0.05 (green)
#>  - [iTRAQ4.116] 116.1116 +/- 0.05 (blue)
#>  - [iTRAQ4.117] 117.115 +/- 0.05 (yellow)
iTRAQ4[1:2]
#> Object of class "ReporterIons"
#> iTRAQ4[1:2]: 'subset of 4-plex iTRAQ' with 2 reporter ions
#>  - [iTRAQ4.114] 114.1112 +/- 0.05 (red)
#>  - [iTRAQ4.115] 115.1083 +/- 0.05 (green)

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
