# Write MS data to mzML or mzXML files

The `writeMSData,MSnExp` and `writeMSData,OnDiskMSnExp` saves the
content of a
[MSnExp](https://lgatto.github.io/MSnbase/reference/MSnExp-class.md) or
[OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
object to MS file(s) in either *mzML* or *mzXML* format.

## Usage

``` r
# S4 method for class 'MSnExp,character'
writeMSData(
  object,
  file,
  outformat = c("mzml", "mzxml"),
  merge = FALSE,
  verbose = isMSnbaseVerbose(),
  copy = FALSE,
  software_processing = NULL
)
```

## Arguments

- object:

  `OnDiskMSnExp` or `MSnExp` object.

- file:

  `character` with the file name(s). Its length has to match the number
  of samples/files of `x`.

- outformat:

  `character(1)` defining the format of the output files. Default output
  format is `"mzml"`.

- merge:

  `logical(1)` whether the data should be saved into a single *mzML*
  file. Default is `merge = FALSE`, i.e. each sample is saved to a
  separate file. **Note**: `merge = TRUE` is not yet implemented.

- verbose:

  `logical(1)` if progress messages should be displayed.

- copy:

  `logical(1)` if metadata (data processings, original file names etc)
  should be copied from the original files. See details for more
  information.

- software_processing:

  optionally provide specific data processing steps. See documentation
  of the `software_processing` parameter of
  [`mzR::writeMSData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html).

## Details

The `writeMSData` method uses the *proteowizard* libraries through the
`mzR` package to save the MS data. The data can be written to *mzML* or
*mzXML* files with or without copying additional metadata information
from the original files from which the data was read by the
[`readMSData()`](https://lgatto.github.io/MSnbase/reference/readMSData.md)
function. This can be set using the `copy` parameter. Note that
`copy = TRUE` requires the original files to be available and is not
supported for input files in other than mzML or mzXML format. All
metadata related to the run is copied, such as instrument information,
data processings etc. If `copy = FALSE` only processing information
performed in R (using `MSnbase`) are saved to the mzML file.

Currently only spectrum data is supported, i.e. if the original mzML
file contains also chromatogram data it is not copied/saved to the new
mzML file.

## Note

General spectrum data such as total ion current, peak count, base peak
m/z or base peak intensity are calculated from the actual spectrum data
before writing the data to the files.

For MSn data, if the `OnDiskMSnExp` or `MSnExp` does not contain also
the precursor scan of a MS level \> 1 spectrum (e.g. due to filtering on
the MS level) `precursorScanNum` is set to 0 in the output file to avoid
potentially linking to a wrong spectrum.

The exported `mzML` file *should* be valid according to the mzML 1.1.2
standard. For exported `mzXML` files it can not be guaranteed that they
are valid and can be opened with other software than `mzR`/`MSnbase`.

## Author

Johannes Rainer
