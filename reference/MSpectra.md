# List of Spectrum objects along with annotations

`MSpectra` (Mass Spectra) objects allow to collect one or more
[Spectrum](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
object(s)
([Spectrum1](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)
or
[Spectrum2](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md))
in a `list`-like structure with the possibility to add arbitrary
annotations to each individual `Spectrum` object. These can be
accessed/set with the
[`S4Vectors::mcols()`](https://rdrr.io/pkg/S4Vectors/man/Vector-class.html)
method.

`MSpectra` objects can be created with the `MSpectra` function.

Functions to access the individual spectra's attributes are available
(listed below).

`writeMgfData` exports a `MSpectra` object to a file in MGF format. All
metadata columns present in `mcols` are exported as additional fields
with the capitalized column names used as field names (see examples
below).

## Usage

``` r
MSpectra(..., elementMetadata = NULL)

# S4 method for class 'MSpectra'
mz(object)

# S4 method for class 'MSpectra'
intensity(object)

# S4 method for class 'MSpectra'
rtime(object)

# S4 method for class 'MSpectra'
precursorMz(object)

# S4 method for class 'MSpectra'
precursorCharge(object)

# S4 method for class 'MSpectra'
precScanNum(object)

# S4 method for class 'MSpectra'
precursorIntensity(object)

# S4 method for class 'MSpectra'
acquisitionNum(object)

# S4 method for class 'MSpectra'
scanIndex(object)

# S4 method for class 'MSpectra,ANY'
peaksCount(object)

# S4 method for class 'MSpectra'
msLevel(object)

# S4 method for class 'MSpectra'
tic(object)

# S4 method for class 'MSpectra'
ionCount(object)

# S4 method for class 'MSpectra'
collisionEnergy(object)

# S4 method for class 'MSpectra'
fromFile(object)

# S4 method for class 'MSpectra'
polarity(object)

# S4 method for class 'MSpectra'
smoothed(object)

# S4 method for class 'MSpectra'
isEmpty(x)

# S4 method for class 'MSpectra'
centroided(object)

# S4 method for class 'MSpectra'
isCentroided(object)

# S4 method for class 'MSpectra'
writeMgfData(object, con = "spectra.mgf", COM = NULL, TITLE = NULL)

# S4 method for class 'MSpectra'
clean(object, all = FALSE, msLevel. = msLevel., ...)

# S4 method for class 'MSpectra'
removePeaks(object, t, msLevel., ...)

# S4 method for class 'MSpectra'
filterMz(object, mz, msLevel., ...)

# S4 method for class 'MSpectra'
pickPeaks(
  object,
  halfWindowSize = 3L,
  method = c("MAD", "SuperSmoother"),
  SNR = 0L,
  refineMz = c("none", "kNeighbors", "kNeighbours", "descendPeak"),
  msLevel. = unique(msLevel(object)),
  ...
)

# S4 method for class 'MSpectra'
smooth(
  x,
  method = c("SavitzkyGolay", "MovingAverage"),
  halfWindowSize = 2L,
  ...
)

# S4 method for class 'MSpectra'
filterMsLevel(object, msLevel.)
```

## Arguments

- ...:

  For `MSpectra`:
  [Spectrum](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
  object(s) or a `list` of
  [Spectrum](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
  objects. For all other methods optional arguments passed along.

- elementMetadata:

  For `MSpectra`:
  [S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with optional information that should be added as metadata information
  (`mcols`) to the object. The number of rows has to match the number of
  [Spectrum](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
  objects, each row is expected to represent additional metadata
  information for one spectrum.

- object:

  For all functions: a `MSpectra` object.

- x:

  For all functions: a `MSpectra` object.

- con:

  For `writeMgfData`: `character(1)` defining the file name of the MGF
  file.

- COM:

  For `writeMgfData`: optional `character(1)` providing a comment to be
  added to the file.

- TITLE:

  For `writeMgfData`: optional `character(1)` defining the title for the
  MGF file.

- all:

  For `clean`: if `FALSE` original 0-intensity values are retained
  around peaks.

- msLevel.:

  For `clean`, `removePeaks`, `filterMz`, `pickPeaks`: optionally
  specify the MS level(s) of the spectra on which the operation should
  be performed. For `filterMsLevels`: MS level(s) to which the
  `MSpectra` should be reduced.

- t:

  For `removePeaks`: `numeric(1)` specifying the threshold below which
  intensities are set to 0.

- mz:

  For `filterMz`: `numeric(2)` defining the lower and upper m/z for the
  filter. See
  [`filterMz()`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
  for details.

- halfWindowSize:

  For `pickPeaks` and `smooth`: see
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md)
  and
  [`smooth()`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md)
  for details.

- method:

  For `pickPeaks` and `smooth`: see
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md)
  and
  [`smooth()`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md)
  for details.

- SNR:

  For `pickPeaks`: see
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md)
  for details.

- refineMz:

  For `pickPeaks`: see
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md)
  for details.

## Details

`MSpectra` inherits all methods from the
[S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
class of the `S4Vectors` package. This includes `lapply` and other data
manipulation and subsetting operations.

## Note

Note that the [Spectra](https://bioconductor.org/packages/Spectra)
package provides a more robust and efficient infrastructure for mass
spectrometry data handling and analysis. So, wherever possible, the
newer *Spectra* package should be used instead of the *MSnbase*.

For backward compatibility, it is however possible to convert between
the `MSpectra` and the newer `Spectra` objects:

- A `Spectra` object can be coerced to a `MSpectra` using
  `as(sps, "MSpectra")` where `sps` is a `Spectra` object.

- The
  [`extractSpectraData()`](https://lgatto.github.io/MSnbase/reference/extractSpectraData.md)
  function can be used to extract the data from a `MSpectra` as a
  `DataFrame`, which can then be used to create a `Spectra` object.

## Constructor

New MSpectra can be created with the `MSpectra(...)` function where
`...` can either be a single
[Spectrum](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)
object or a `list` of `Spectrum` objects
([Spectrum1](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)
and/or
[Spectrum2](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)).

## Accessing spectrum attributes

These methods allow to access the attributes and values of the
individual `Spectrum`
([Spectrum1](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)
or
[Spectrum2](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md))
objects within the list.

- `mz` return the m/z values of each spectrum as a `list` of `numeric`
  vectors.

- `intensity` return the intensity values of each spectrum as a `list`
  of `numeric` vectors.

- `rtime` return the retention time of each spectrum as a `numeric`
  vector with length equal to the length of `object`.

- `precursorMz`, `precursorCharge`, `precursorIntensity`, `precScanNum`
  return precursor m/z values, charge, intensity and scan number for
  each spectrum as a `numeric` (or `integer`) vector with length equal
  to the length of `object`. Note that for
  [Spectrum1](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)
  objects `NA` will be returned.

- `acquisitionNum` and `scanIndex` return the acquisition number of each
  spectrum and its scan index as an `integer` vector with the same
  length than `object`.

- `ionCount` and `tic` return the ion count and total ion current of
  each spectrum.

- `peaksCount` returns the number of peaks for each spectrum as a
  `integer` vector.

- `msLevel` returns the MS level of each spectrum.

- `collisionEnergy` returns the collision energy for each spectrum or
  `NA` for
  [Spectrum1](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)
  objects.

- `polarity` returns the spectra's polarity.

- `fromFile` returns the index from the (e.g. mzML) file the spectra
  where from. This applies only for spectra read using the
  [`readMSData()`](https://lgatto.github.io/MSnbase/reference/readMSData.md)
  function.

- `smoothed` whether spectra have been smoothed (i.e. processed with the
  [`smooth()`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md)
  method. Returns a `logical` of length equal to the number of spectra.

- `isEmpty` returns `TRUE` for spectra without peak data.

- `centroided`, `isCentroided` returns for each spectrum whether it
  contains *centroided* data. While `centroided` returns the internal
  attribute of each spectrum, `isCentroided` tries to guess whether
  spectra are centroided from the actual peak data.

## Data manipulation methods

- `clean` *cleans* each spectrum. See
  [`clean()`](https://lgatto.github.io/MSnbase/reference/clean-methods.md)
  for more details.

- `pickPeaks` performs peak picking to generate centroided spectra. See
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.md)
  for more details.

- `removePeaks` removes peaks lower than a threshold `t`. See
  [`removePeaks()`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.md)
  for more details.

- `smooth` *smooths* spectra. See
  [`smooth()`](https://lgatto.github.io/MSnbase/reference/smooth-methods.md)
  for more details.

## Filtering and subsetting

- `[` can be used to subset the `MSpectra` object.

- `filterMsLevel` filters `MSpectra` to retain only spectra from certain
  MS level(s).

- `filterMz` filters the spectra by the specified `mz` range. See
  [`filterMz()`](https://lgatto.github.io/MSnbase/reference/trimMz-methods.md)
  for details.

## Author

Johannes Rainer

## Examples

``` r

## Create from Spectrum objects
sp1 <- new("Spectrum1", mz = c(1, 2, 4), intensity = c(4, 5, 2))
sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
    precursorMz = 2)

spl <- MSpectra(sp1, sp2)
spl
#> MSpectra with 2 spectra and 0 metadata column(s):
#>       msLevel     rtime peaksCount
#>     <integer> <numeric>  <integer>
#>   1         1        NA          3
#>   2         2        NA          4
spl[[1]]
#> Object of class "Spectrum1"
#>  MSn level: 1 
#>  Total ion count: 3 
#>  Polarity:  

## Add also metadata columns
mcols(spl)$id <- c("a", "b")
mcols(spl)
#> DataFrame with 2 rows and 1 column
#>            id
#>   <character>
#> 1           a
#> 2           b

## Create a MSpectra with metadata
spl <- MSpectra(sp1, sp2, elementMetadata = DataFrame(id = c("a", "b")))

mcols(spl)
#> DataFrame with 2 rows and 1 column
#>            id
#>   <character>
#> 1           a
#> 2           b
mcols(spl)$id
#> [1] "a" "b"

## Extract the mz values for the individual spectra
mz(spl)
#> $`1`
#> [1] 1 2 4
#> 
#> $`2`
#> [1] 1 2 3 4
#> 

## Extract the intensity values for the individual spectra
intensity(spl)
#> $`1`
#> [1] 4 5 2
#> 
#> $`2`
#> [1] 5 3 2 5
#> 

## Extract the retention time values for the individual spectra
rtime(spl)
#>  1  2 
#> NA NA 

## Extract the precursor m/z of each spectrum.
precursorMz(spl)
#>  1  2 
#> NA  2 

## Extract the precursor charge of each spectrum.
precursorCharge(spl)
#>  1  2 
#> NA NA 

## Extract the precursor scan number for each spectrum.
precScanNum(spl)
#>  1  2 
#> NA NA 

## Extract the precursor intensity of each spectrum.
precursorIntensity(spl)
#>  1  2 
#> NA NA 

## Extract the acquisition number of each spectrum.
acquisitionNum(spl)
#>  1  2 
#> NA NA 

## Extract the scan index of each spectrum.
scanIndex(spl)
#>  1  2 
#> NA NA 

## Get the number of peaks per spectrum.
peaksCount(spl)
#> 1 2 
#> 3 4 

## Get the MS level of each spectrum.
msLevel(spl)
#> 1 2 
#> 1 2 

## Get the total ion current for each spectrum.
tic(spl)
#>  1  2 
#> 11 15 

## Get the total ion current for each spectrum.
ionCount(spl)
#>  1  2 
#> 11 15 

## Extract the collision energy for each spectrum.
collisionEnergy(spl)
#>  1  2 
#> NA NA 

## Extract the file index for each spectrum.
fromFile(spl)
#>  1  2 
#> NA NA 

## Get the polarity for each spectrum.
polarity(spl)
#>  1  2 
#> NA NA 

## Whether spectra are smoothed (i.e. processed with the `smooth`
## function).
smoothed(spl)
#>  1  2 
#> NA NA 

## Are spectra empty (i.e. contain no peak data)?
isEmpty(spl)
#>     1     2 
#> FALSE FALSE 

## Do the spectra contain centroided data?
centroided(spl)
#>  1  2 
#> NA NA 

## Do the spectra contain centroided data? Whether spectra are centroided
## is estimated from the peak data.
isCentroided(spl)
#>  1  2 
#> NA NA 

## Export the spectrum list to a MGF file. Values in metadata columns are
## exported as additional field for each spectrum.
tmpf <- tempfile()
writeMgfData(spl, tmpf)

## Evaluate the written output. The ID of each spectrum (defined in the
## "id" metadata column) is exported as field "ID".
readLines(tmpf)
#>  [1] "COM=Experimentexported by MSnbase on Fri Apr 10 15:48:24 2026"   
#>  [2] "BEGIN IONS"                                                      
#>  [3] "SCANS=NA"                                                        
#>  [4] "TITLE=msLevel 1; retentionTime ; scanNum NA"                     
#>  [5] "RTINSECONDS="                                                    
#>  [6] "ID=a"                                                            
#>  [7] "1 4"                                                             
#>  [8] "2 5"                                                             
#>  [9] "4 2"                                                             
#> [10] "END IONS"                                                        
#> [11] ""                                                                
#> [12] "BEGIN IONS"                                                      
#> [13] "SCANS="                                                          
#> [14] "TITLE=msLevel 2; retentionTime ; scanNum ; precMz 2; precCharge "
#> [15] "RTINSECONDS="                                                    
#> [16] "PEPMASS=2"                                                       
#> [17] "ID=b"                                                            
#> [18] "1 5"                                                             
#> [19] "2 3"                                                             
#> [20] "3 2"                                                             
#> [21] "4 5"                                                             
#> [22] "END IONS"                                                        

## Set mcols to NULL to avoid export of additional data fields.
mcols(spl) <- NULL
file.remove(tmpf)
#> [1] TRUE

writeMgfData(spl, tmpf)
readLines(tmpf)
#>  [1] "COM=Experimentexported by MSnbase on Fri Apr 10 15:48:24 2026"   
#>  [2] "BEGIN IONS"                                                      
#>  [3] "SCANS=NA"                                                        
#>  [4] "TITLE=msLevel 1; retentionTime ; scanNum NA"                     
#>  [5] "RTINSECONDS="                                                    
#>  [6] "1 4"                                                             
#>  [7] "2 5"                                                             
#>  [8] "4 2"                                                             
#>  [9] "END IONS"                                                        
#> [10] ""                                                                
#> [11] "BEGIN IONS"                                                      
#> [12] "SCANS="                                                          
#> [13] "TITLE=msLevel 2; retentionTime ; scanNum ; precMz 2; precCharge "
#> [14] "RTINSECONDS="                                                    
#> [15] "PEPMASS=2"                                                       
#> [16] "1 5"                                                             
#> [17] "2 3"                                                             
#> [18] "3 2"                                                             
#> [19] "4 5"                                                             
#> [20] "END IONS"                                                        

## Filter the object by MS level
filterMsLevel(spl, msLevel. = 1)
#> MSpectra with 1 spectra and 0 metadata column(s):
#>       msLevel     rtime peaksCount
#>     <integer> <numeric>  <integer>
#>   1         1        NA          3
```
