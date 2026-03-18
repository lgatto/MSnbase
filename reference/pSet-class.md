# Class to Contain Raw Mass-Spectrometry Assays and Experimental Metadata

Container for high-throughput mass-spectrometry assays and experimental
metadata. This class is based on Biobase's
`"`[`eSet`](https://rdrr.io/pkg/Biobase/man/class.eSet.html)`"` virtual
class, with the notable exception that 'assayData' slot is an
environment contain objects of class
`"`[`Spectrum`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.md)`"`.

## Objects from the Class

A virtual Class: No objects may be created from it. See
`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
for instantiatable sub-classes.

## Slots

- `assayData`::

  Object of class `"environment"` containing the MS spectra (see
  `"`[`Spectrum1`](https://lgatto.github.io/MSnbase/reference/Spectrum1-class.md)`"`
  and
  `"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`).

- `phenoData`::

  Object of class
  `"`[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`"`
  containing experimenter-supplied variables describing sample (i.e the
  individual tags for an labelled MS experiment) See
  [`phenoData`](https://rdrr.io/pkg/Biobase/man/phenoData.html) for more
  details.

- `featureData`::

  Object of class
  `"`[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`"`
  containing variables describing features (spectra in our case), e.g.
  identificaiton data, peptide sequence, identification score,...
  (inherited from
  `"`[`eSet`](https://rdrr.io/pkg/Biobase/man/class.eSet.html)`"`). See
  [`featureData`](https://rdrr.io/pkg/Biobase/man/featureData.html) for
  more details.

- `experimentData`::

  Object of class
  `"`[`MIAPE`](https://lgatto.github.io/MSnbase/reference/MIAPE-class.md)`"`,
  containing details of experimental methods. See
  [`experimentData`](https://rdrr.io/pkg/Biobase/man/abstract.html) for
  more details.

- `protocolData`::

  Object of class
  `"`[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`"`
  containing equipment-generated variables (inherited from
  `"`[`eSet`](https://rdrr.io/pkg/Biobase/man/class.eSet.html)`"`). See
  [`protocolData`](https://rdrr.io/pkg/Biobase/man/protocolData.html)
  for more details.

- `processingData`::

  Object of class
  `"`[`MSnProcess`](https://lgatto.github.io/MSnbase/reference/MSnProcess-class.md)`"`
  that records all processing.

- `.cache`::

  Object of class `environment` used to cache data. Under development.

- `.__classVersion__`::

  Object of class
  `"`[`Versions`](https://rdrr.io/pkg/Biobase/man/class.Versions.html)`"`
  describing the versions of the class.

## Extends

Class
`"`[`VersionedBiobase`](https://rdrr.io/pkg/Biobase/man/class.VersionedBiobase.html)`"`,
directly. Class
`"`[`Versioned`](https://rdrr.io/pkg/Biobase/man/class.Versioned.html)`"`,
by class "VersionedBiobase", distance 2.

## Methods

Methods defined in derived classes may override the methods described
here.

- \[:

  `signature(x = "pSet")`: Subset current object and return object of
  same class.

- \[\[:

  `signature(x = "pSet")`: Direct access to individual spectra.

- \$:

  `signature(x = "pSet")`: directly access a specific sample annotation
  column from the `pData`.

- \$\<-:

  `signature(x = "pSet")`: replace or add a sample annotation column in
  the `pData`.

- abstract:

  Access abstract in `experimentData`.

- assayData:

  `signature(object = "pSet")`: Access the `assayData` slot. Returns an
  `environment`.

- desciption:

  `signature(x = "pSet")`: Synonymous with experimentData.

- dim:

  `signature(x = "pSet")`: Returns the dimensions of the `phenoData`
  slot.

- experimentData:

  `signature(x = "pSet")`: Access details of experimental methods.

- featureData:

  `signature(x = "pSet")`: Access the `featureData` slot.

- fData:

  `signature(x = "pSet")`: Access feature data information.

- featureNames:

  `signature(x = "pSet")`: Coordinate access of feature names (e.g
  spectra, peptides or proteins) in `assayData` slot.

- fileNames:

  `signature(object = "pSet")`: Access file names in the
  `processingData` slot.

- fromFile:

  `signature(object = "pSet")`: Access raw data file indexes (to be
  found in the `processingData` slot) from which the individual object's
  spectra where read from.

- centroided:

  `signature(object = "pSet")`: Indicates whether individual spectra are
  centroided ('TRUE') of uncentroided ('FALSE'). Use
  `centroided(object) <- value` to update a whole experiment, ensuring
  that `object` and `value` have the same length.

- smoothed:

  `signature(object = "pSet")`: Indicates whether individual spectra are
  smoothed ('TRUE') of unsmoothed ('FALSE'). Use
  `smoothed(object) <- value` to update a whole experiment, ensuring
  that `object` and `value` have the same length.

- fvarMetadata:

  `signature(x = "pSet")`: Access metadata describing features reported
  in `fData`.

- fvarLabels:

  `signature(x = "pSet")`: Access variable labels in `featureData`.

- length:

  `signature(x = "pSet")`: Returns the number of features in the
  `assayData` slot.

- notes:

  `signature(x = "pSet")`: Retrieve and unstructured notes associated
  with `pSet` in the `experimentData` slot.

- pData:

  `signature(x = "pSet")`: Access sample data information.

- pData\<-:

  `signature(x = "pSet", value)`: Replace sample data information with
  `value`, value being a `data.frame`.

- phenoData:

  `signature(x = "pSet")`: Access the `phenoData` slot.

- phenoData\<-:

  `signature(x = "pSet", value)`: Replace sample data information with
  `value`. `value` can be a `data.frame` or an `AnnotatedDataFrame`.

- processingData:

  `signature(object = "pSet")`: Access the `processingData` slot.

- protocolData:

  `signature(x = "pSet")`: Access the `protocolData` slot.

- pubMedIds:

  `signature(x = "pSet")`: Access PMIDs in `experimentData`.

- sampleNames:

  `signature(x = "pSet")`: Access sample names in `phenoData`. A
  replacement method is also available.

- spectra:

  `signature(x = "pSet", ...)`: Access the `assayData` slot, returning
  the features as a `list`. Additional arguments are currently ignored.

- varMetadata:

  `signature(x = "pSet")`: Access metadata describing variables reported
  in `pData`.

- varLabels:

  `signature(x = "pSet")`: Access variable labels in `phenoData`.

- acquisitionNum:

  `signature(object = "pSet")`: Accessor for spectra acquisition
  numbers.

- scanIndex:

  `signature(object = "pSet")`: Accessor for spectra scan indices.

- collisionEnergy:

  `signature(object = "pSet")`: Accessor for MS2 spectra collision
  energies.

- intensity:

  `signature(object = "pSet", ...)`: Accessor for spectra instenities,
  returned as named list. Additional arguments are currently ignored.

- msInfo:

  `signature(object = "pSet")`: Prints the MIAPE-MS meta-data stored in
  the `experimentData` slot.

- msLevel:

  `signature(object = "pSet")`: Accessor for spectra MS levels.

- mz:

  `signature(object = "pSet", ...)`: Accessor for spectra M/Z values,
  returned as a named list. Additional arguments are currently ignored.

- peaksCount:

  `signature(object = "pSet")`: Accessor for spectra preak counts.

- peaksCount:

  `signature(object = "pSet", scans = "numeric")`: Accessor to `scans`
  spectra preak counts.

- polarity:

  `signature(object = "pSet")`: Accessor for MS1 spectra polarities.

- precursorCharge:

  `signature(object = "pSet")`: Accessor for MS2 precursor charges.

- precursorIntensity:

  `signature(object = "pSet")`: Accessor for MS2 precursor intensity.

- precursorMz:

  `signature(object = "pSet")`: Accessor for MS2 precursor M/Z values.

- precAcquisitionNum:

  `signature(object = "pSet")`: Accessor for MS2 precursor scan numbers.

- precScanNum:

  see `precAcquisitionNum`.

- rtime:

  `signature(object = "pSet", ...)`: Accessor for spectra retention
  times. Additional arguments are currently ignored.

- tic:

  `signature(object = "pSet", ...)`: Accessor for spectra total ion
  counts. Additional arguments are currently ignored.

- ionCount:

  `signature(object = "pSet")`: Accessor for spectra total ion current.

- header:

  `signature(object = "pSet")`: Returns a data frame containing all
  available spectra parameters (MSn only).

- header:

  `signature(object = "pSet", scans = "numeric")`: Returns a data frame
  containing `scans` spectra parameters (MSn only).

- spectrapply:

  `spectrapply(signature(object = "pSet"), FUN = NULL, BPPARAM = bpparam(), ...)`:
  applies the function `FUN` to each spectrum passing additional
  parameters in `...` to that function and return its results. For
  `FUN = NULL` it returns the list of spectra (same as a call to
  `spectra`). Parameter `BPPARAM` allows to specify how and if parallel
  processing should be enabled.

  Returns a list with the result for each of spectrum.

- isolationWindowLowerMz:

  `isolationWindowLowerMz(object = "pSet")`: return the lower m/z
  boundary for the isolation window. Note that this method is at present
  only available for
  [`OnDiskMSnExp`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
  objects.

- isolationWindowUpperMz:

  `isolationWindowUpperMz(object = "pSet")`: return the upper m/z
  boundary for the isolation window. Note that this method is at present
  only available for
  [`OnDiskMSnExp`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
  objects.

Additional accessors for the experimental metadata (`experimentData`
slot) are defined. See
`"`[`MIAPE`](https://lgatto.github.io/MSnbase/reference/MIAPE-class.md)`"`
for details.

## References

The `"`[`eSet`](https://rdrr.io/pkg/Biobase/man/class.eSet.html)`"`
class, on which `pSet` is based.

## Author

Laurent Gatto

## See also

`"`[`MSnExp`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md)`"`
for an instantiatable application of `pSet`.

## Examples

``` r
showClass("pSet")
#> Virtual Class "pSet" [package "MSnbase"]
#> 
#> Slots:
#>                                                                
#> Name:           assayData          phenoData        featureData
#> Class:        environment AnnotatedDataFrame AnnotatedDataFrame
#>                                                                
#> Name:      experimentData       protocolData     processingData
#> Class:              MIAxE AnnotatedDataFrame         MSnProcess
#>                                             
#> Name:              .cache  .__classVersion__
#> Class:        environment           Versions
#> 
#> Extends: "Versioned"
#> 
#> Known Subclasses: 
#> Class "MSnExp", directly
#> Class "OnDiskMSnExp", by class "MSnExp", distance 2, with explicit coerce
```
