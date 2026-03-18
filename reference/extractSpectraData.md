# Conversion between objects from the Spectra and MSnbase packages

The [Spectra](https://bioconductor.org/packages/Spectra) package
provides a more robust and efficient infrastructure for mass
spectrometry data handling and analysis. So, wherever possible, the
newer *Spectra* package should be used instead of the *MSnbase*. The
functions listed here allow to convert between objects from the
*MSnbase* and *Spectra* packages.

`extractSpectraData` extracts the spectra data (m/z and intensity values
including metadata) from
[MSnExp](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md),
[OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md),
[Spectrum1](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md),
[Spectrum2](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)
objects (or `list` of such objects) and returns these as a `DataFrame`
that can be used to create a
[Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object.This function enables thus to convert data from the *old*
`MSnbase` package to the newer `Spectra` package.

To convert a `Spectra` object to a `MSpectra` object use
`as(sps, "MSpectra")` where `sps` is a `Spectra` object.

## Usage

``` r
extractSpectraData(x)
```

## Arguments

- x:

  a `list` of
  [Spectrum](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
  objects or an object extending
  [MSnExp](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)
  or a
  [MSpectra](https://lgatto.github.io/MSnbase/reference/MSpectra.md)
  object.

## Value

- `extracSpectraData()` returns a
  [`DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with the full spectrum data that can be passed to the
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  function to create a `Spectra` object.

- `as(x, "MSpectra")` returns a `MSpectra` object with the content of
  the `Spectra` object `x`.

## Note

Coercion from `Spectra` to a `MSpectra` will only assign values to the
contained `Spectrum1` and `Spectrum2` objects, but will not add all
eventually spectra variables present in `Spectra`.

## Author

Johannes Rainer

## Examples

``` r

## Read an mzML file with MSnbase
fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
    package = "msdata")
data <- filterRt(readMSData(fl, mode = "onDisk"), rt = c(1, 6))

## Extract the data as a DataFrame
res <- extractSpectraData(data)
res
#> DataFrame with 50 rows and 38 columns
#>           fromFile     spIdx  smoothed scanIndex acquisitionNum   msLevel
#>          <integer> <integer> <logical> <integer>      <integer> <integer>
#> F1.S0010         1        10        NA        10             10         2
#> F1.S0011         1        11        NA        11             11         2
#> F1.S0012         1        12        NA        12             12         2
#> F1.S0013         1        13        NA        13             13         2
#> F1.S0014         1        14        NA        14             14         2
#> ...            ...       ...       ...       ...            ...       ...
#> F1.S0055         1        55        NA        55             55         2
#> F1.S0056         1        56        NA        56             56         2
#> F1.S0057         1        57        NA        57             57         2
#> F1.S0058         1        58        NA        58             58         2
#> F1.S0059         1        59        NA        59             59         2
#>           polarity originalPeaksCount totIonCurrent     rtime basePeakMZ
#>          <integer>          <integer>     <numeric> <numeric>  <numeric>
#> F1.S0010         1                340         83853     1.101    141.946
#> F1.S0011         1                491         61072     1.198    217.106
#> F1.S0012         1                296         28639     1.295    217.108
#> F1.S0013         1                301         26285     1.392    149.024
#> F1.S0014         1                292         27660     1.489    281.049
#> ...            ...                ...           ...       ...        ...
#> F1.S0055         1                373         89856     5.601    141.946
#> F1.S0056         1                532         66566     5.698    217.108
#> F1.S0057         1                301         30835     5.795    217.106
#> F1.S0058         1                256         23034     5.892    279.098
#> F1.S0059         1                310         27343     5.989    281.051
#>          basePeakIntensity collisionEnergy electronBeamEnergy ionisationEnergy
#>                  <numeric>       <numeric>          <numeric>        <numeric>
#> F1.S0010              6580               0                 NA                0
#> F1.S0011              1429               0                 NA                0
#> F1.S0012               501               0                 NA                0
#> F1.S0013               413               0                 NA                0
#> F1.S0014               714               0                 NA                0
#> ...                    ...             ...                ...              ...
#> F1.S0055              7477               0                 NA                0
#> F1.S0056              1400               0                 NA                0
#> F1.S0057               432               0                 NA                0
#> F1.S0058               186               0                 NA                0
#> F1.S0059               496               0                 NA                0
#>              lowMZ    highMZ precursorScanNum precursorMZ precursorCharge
#>          <numeric> <numeric>        <integer>   <numeric>       <integer>
#> F1.S0010         0         0               NA      163.75               0
#> F1.S0011         0         0               NA      208.95               0
#> F1.S0012         0         0               NA      244.05               0
#> F1.S0013         0         0               NA      270.85               0
#> F1.S0014         0         0               NA      299.10               0
#> ...            ...       ...              ...         ...             ...
#> F1.S0055         0         0               NA      163.75               0
#> F1.S0056         0         0               NA      208.95               0
#> F1.S0057         0         0               NA      244.05               0
#> F1.S0058         0         0               NA      270.85               0
#> F1.S0059         0         0               NA      299.10               0
#>          precursorIntensity mergedScan mergedResultScanNum
#>                   <numeric>  <integer>           <integer>
#> F1.S0010                  0         NA                  NA
#> F1.S0011                  0         NA                  NA
#> F1.S0012                  0         NA                  NA
#> F1.S0013                  0         NA                  NA
#> F1.S0014                  0         NA                  NA
#> ...                     ...        ...                 ...
#> F1.S0055                  0         NA                  NA
#> F1.S0056                  0         NA                  NA
#> F1.S0057                  0         NA                  NA
#> F1.S0058                  0         NA                  NA
#> F1.S0059                  0         NA                  NA
#>          mergedResultStartScanNum mergedResultEndScanNum injectionTime
#>                         <integer>              <integer>     <numeric>
#> F1.S0010                       NA                     NA             0
#> F1.S0011                       NA                     NA             0
#> F1.S0012                       NA                     NA             0
#> F1.S0013                       NA                     NA             0
#> F1.S0014                       NA                     NA             0
#> ...                           ...                    ...           ...
#> F1.S0055                       NA                     NA             0
#> F1.S0056                       NA                     NA             0
#> F1.S0057                       NA                     NA             0
#> F1.S0058                       NA                     NA             0
#> F1.S0059                       NA                     NA             0
#>          filterString    spectrumId centroided ionMobilityDriftTime
#>           <character>   <character>  <logical>            <numeric>
#> F1.S0010           NA sample=1 p...       TRUE                   NA
#> F1.S0011           NA sample=1 p...       TRUE                   NA
#> F1.S0012           NA sample=1 p...       TRUE                   NA
#> F1.S0013           NA sample=1 p...       TRUE                   NA
#> F1.S0014           NA sample=1 p...       TRUE                   NA
#> ...               ...           ...        ...                  ...
#> F1.S0055           NA sample=1 p...       TRUE                   NA
#> F1.S0056           NA sample=1 p...       TRUE                   NA
#> F1.S0057           NA sample=1 p...       TRUE                   NA
#> F1.S0058           NA sample=1 p...       TRUE                   NA
#> F1.S0059           NA sample=1 p...       TRUE                   NA
#>          isolationWindowTargetMZ isolationWindowLowerOffset
#>                        <numeric>                  <numeric>
#> F1.S0010                  163.75                      24.25
#> F1.S0011                  208.95                      21.95
#> F1.S0012                  244.05                      14.15
#> F1.S0013                  270.85                      13.65
#> F1.S0014                  299.10                      15.60
#> ...                          ...                        ...
#> F1.S0055                  163.75                      24.25
#> F1.S0056                  208.95                      21.95
#> F1.S0057                  244.05                      14.15
#> F1.S0058                  270.85                      13.65
#> F1.S0059                  299.10                      15.60
#>          isolationWindowUpperOffset scanWindowLowerLimit scanWindowUpperLimit
#>                           <numeric>            <numeric>            <numeric>
#> F1.S0010                      24.25                  140                  900
#> F1.S0011                      21.95                  140                  900
#> F1.S0012                      14.15                  140                  900
#> F1.S0013                      13.65                  140                  900
#> F1.S0014                      15.60                  140                  900
#> ...                             ...                  ...                  ...
#> F1.S0055                      24.25                  140                  900
#> F1.S0056                      21.95                  140                  900
#> F1.S0057                      14.15                  140                  900
#> F1.S0058                      13.65                  140                  900
#> F1.S0059                      15.60                  140                  900
#>           spectrum                          mz
#>          <integer>               <NumericList>
#> F1.S0010        10 140.888,140.949,140.974,...
#> F1.S0011        11 140.901,140.913,140.958,...
#> F1.S0012        12 140.992,141.090,141.104,...
#> F1.S0013        13 140.872,141.957,143.089,...
#> F1.S0014        14 140.872,140.960,141.970,...
#> ...            ...                         ...
#> F1.S0055        55 140.949,140.967,140.983,...
#> F1.S0056        56 140.045,140.912,140.955,...
#> F1.S0057        57 141.094,141.963,141.981,...
#> F1.S0058        58 140.876,141.955,141.972,...
#> F1.S0059        59 141.098,141.957,141.972,...
#>                                  intensity
#>                              <NumericList>
#> F1.S0010 0.0999316,2.6204529,0.5139834,...
#> F1.S0011 0.0999386,0.1832319,0.5464582,...
#> F1.S0012 0.0599841,0.2358523,0.0800102,...
#> F1.S0013 0.0199862,0.0200630,0.0604285,...
#> F1.S0014 0.0799445,0.0199924,0.0200639,...
#> ...                                    ...
#> F1.S0055    1.384358,0.183261,0.337385,...
#> F1.S0056 0.0597824,0.1823982,0.2357384,...
#> F1.S0057    0.288360,0.120379,0.100323,...
#> F1.S0058 0.0999327,0.1237221,0.0200640,...
#> F1.S0059    0.180017,0.020063,0.140449,...

library(Spectra)
#> 
#> Attaching package: ‘Spectra’
#> The following objects are masked from ‘package:MSnbase’:
#> 
#>     combineSpectra, pickPeaks, plotMzDelta
## This can be used as an input for the Spectra constructor of the
## Spectra package:
sps <- Spectra::Spectra(res)
sps
#> MSn data (Spectra) with 50 spectra in a MsBackendMemory backend:
#>            msLevel     rtime scanIndex
#>          <integer> <numeric> <integer>
#> F1.S0010         2     1.101        10
#> F1.S0011         2     1.198        11
#> F1.S0012         2     1.295        12
#> F1.S0013         2     1.392        13
#> F1.S0014         2     1.489        14
#> ...            ...       ...       ...
#> F1.S0055         2     5.601        55
#> F1.S0056         2     5.698        56
#> F1.S0057         2     5.795        57
#> F1.S0058         2     5.892        58
#> F1.S0059         2     5.989        59
#>  ... 42 more variables/columns.

## A Spectra object can be coerced to a MSnbase MSpectra object using
msps <- as(sps, "MSpectra")
```
