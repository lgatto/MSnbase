# Read 'MSnSet'

This function reads data files to generate an
[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
instance. It is a wrapper around `Biobase`'s
[`readExpressionSet`](https://rdrr.io/pkg/Biobase/man/readExpressionSet.html)
function with an additional `featureDataFile` parameter to include
feature data. See also
[`readExpressionSet`](https://rdrr.io/pkg/Biobase/man/readExpressionSet.html)
for more details. `readMSnSet2` is a simple version that takes a single
text spreadsheet as input and extracts the expression data and feature
meta-data to create and `MSnSet`.

Note that when using `readMSnSet2`, one should not set `rownames` as
additional argument to defined feature names. It is ignored and used to
set `fnames` if not provided otherwise.

## Usage

``` r
readMSnSet(exprsFile,
           phenoDataFile,
           featureDataFile,
           experimentDataFile,
           notesFile,
           path, annotation,
           exprsArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, ...),
           phenoDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, stringsAsFactors = stringsAsFactors, ...),
           featureDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, stringsAsFactors = stringsAsFactors, ...),
           experimentDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, stringsAsFactors = stringsAsFactors, ...),
           sep = "\t",
           header = TRUE,
           quote = "",
           stringsAsFactors = FALSE,
           row.names = 1L,
           widget = getOption("BioC")$Base$use.widgets, ...)

readMSnSet2(file, ecol, fnames, ...)
```

## Arguments

Arguments direclty passed to `readExpressionSet`. The description is
from the `readExpressionSet` documentation page.

- exprsFile:

  (character) File or connection from which to read expression values.
  The file should contain a matrix with rows as features and columns as
  samples. [`read.table`](https://rdrr.io/r/utils/read.table.html) is
  called with this as its `file` argument and further arguments given by
  `exprsArgs`.

- phenoDataFile:

  (character) File or connection from which to read phenotypic data.
  [`read.AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/read.AnnotatedDataFrame.html)
  is called with this as its `file` argument and further arguments given
  by `phenoDataArgs`.

- experimentDataFile:

  (character) File or connection from which to read experiment data.
  [`read.MIAME`](https://rdrr.io/pkg/Biobase/man/read.MIAME.html) is
  called with this as its `file` argument and further arguments given by
  `experimentDataArgs`.

- notesFile:

  (character) File or connection from which to read notes;
  [`readLines`](https://rdrr.io/r/base/readLines.html) is used to input
  the file.

- path:

  (optional) directory in which to find all the above files.

- annotation:

  (character) A single character string indicating the annotation
  associated with this ExpressionSet.

- exprsArgs:

  A list of arguments to be used with
  [`read.table`](https://rdrr.io/r/utils/read.table.html) when reading
  in the expression matrix.

- phenoDataArgs:

  A list of arguments to be used (with
  [`read.AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/read.AnnotatedDataFrame.html))
  when reading the phenotypic data.

- experimentDataArgs:

  A list of arguments to be used (with
  [`read.MIAME`](https://rdrr.io/pkg/Biobase/man/read.MIAME.html)) when
  reading the experiment data.

- sep, header, quote, stringsAsFactors, row.names:

  arguments used by the
  [`read.table`](https://rdrr.io/r/utils/read.table.html)-like
  functions.

- widget:

  A boolean value indicating whether widgets can be used. Widgets are
  NOT yet implemented for `read.AnnotatedDataFrame`.

- ...:

  Further arguments that can be passed on to the
  [`read.table`](https://rdrr.io/r/utils/read.table.html)-like
  functions.

Additional argument, specific to `readMSnSet`:

- featureDataFile:

  (character) File or connection from which to read feature data.
  [`read.AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/read.AnnotatedDataFrame.html)
  is called with this as its `file` argument and further arguments given
  by `phenoDataArgs`.

- featureDataArgs:

  A list of arguments to be used (with
  [`read.AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/read.AnnotatedDataFrame.html))
  when reading the phenotypic data.

Arguments for `readMSnSet2`:

- file:

  A `character` indicating the spreadsheet file or a `data.frame` (new
  in version 1.19.8). Default, when `file` is a `character`, is to read
  the file as a comma-separated values (csv). If different, use the
  additional arguments, passed to
  [`read.csv`](https://rdrr.io/r/utils/read.table.html), to parametrise
  file import.

  Passing a `data.frame` can be particularly useful if the spreadsheet
  is in Excel format. The appropriate sheet can first be read into R as
  a `data.frame` using, for example `readxl::read_excel`, and then pass
  it to `readMSnSet2`.

- ecol:

  A `numeric` indicating the indices of the columns to be used as
  expression values. Can also be a `character` indicating the names of
  the columns. Caution must be taken if the column names are composed of
  special characters like `(` or `-` that will be converted to a `.`. If
  `ecol` does not match, the error message will dislpay the column names
  are see by `R`.

- fnames:

  An optional `character` or `numeric` of length 1 indicating the column
  to be used as feature names.

## Value

An instance of the
[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)
class.

## Author

Laurent Gatto \<laurent.gatto@uclouvain.be\>

## See also

The
[`grepEcols`](https://lgatto.github.io/MSnbase/reference/grepEcols.md)
and
[`getEcols`](https://lgatto.github.io/MSnbase/reference/grepEcols.md)
helper functions to identify the `ecol` values. The `MSnbase-io`
vignette illustrates these functions in detail. It can be accessed with
`vignette("MSnbase-io")`.

## Examples

``` r
if (FALSE) { # \dontrun{
exprsFile <- "path_to_intensity_file.csv"
fdatafile <- "path_to_featuredata_file.csv"
pdatafile <- "path_to_sampledata_file.csv"
## Read ExpressionSet with appropriate parameters
res <- readMSnSet(exprsFile, pdataFile, fdataFile, sep = "\t", header=TRUE)
} # }

library("pRolocdata")
f0 <- dir(system.file("extdata", package = "pRolocdata"),
          full.names = TRUE,
          pattern = "Dunkley2006")
basename(f0)
#> [1] "Dunkley2006.csv.gz"
res <- readMSnSet2(f0, ecol = 5:20)
res
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 689 features, 16 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData: none
#> featureData
#>   featureNames: 1 2 ... 689 (689 total)
#>   fvarLabels: Protein.ID Loc.Predicted ... pd.markers (6 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:  
#> - - - Processing information - - -
#>  MSnbase version: 2.37.2 
head(exprs(res)) ## columns 5 to 20
#>      M1F1A    M1F4A    M1F7A   M1F11A    M1F2B    M1F5B    M1F8B   M1F11B
#> 1 0.323250 0.275500 0.216000 0.185250 0.465667 0.199667 0.174333 0.160333
#> 2 0.332000 0.279667 0.222000 0.166000 0.451500 0.200375 0.177250 0.171125
#> 3 0.397250 0.246500 0.168250 0.188250 0.459750 0.214500 0.183250 0.142250
#> 4 0.336733 0.303267 0.201133 0.158800 0.487167 0.201833 0.165333 0.145333
#> 5 0.328800 0.302900 0.192100 0.176400 0.542000 0.180750 0.151250 0.126250
#> 6 0.343714 0.295286 0.195000 0.165714 0.500111 0.207000 0.160333 0.132444
#>      M2F1A    M2F4A    M2F7A   M2F11A    M2F2B    M2F5B    M2F8B M2F11B
#> 1 0.370667 0.317444 0.154333 0.157444 0.379500 0.333000 0.161000 0.1270
#> 2 0.371923 0.290923 0.168000 0.169154 0.428800 0.285600 0.153000 0.1328
#> 3 0.390200 0.298400 0.176800 0.134200 0.413500 0.255000 0.172500 0.1590
#> 4 0.387833 0.326833 0.152667 0.133000 0.416333 0.305667 0.147333 0.1310
#> 5 0.356714 0.306857 0.172143 0.164286 0.450333 0.260667 0.158667 0.1300
#> 6 0.344087 0.326739 0.177609 0.151739 0.373750 0.295375 0.176625 0.1540
head(fData(res)) ## other columns
#>   Protein.ID Loc.Predicted Loc.Confirmed Loc.Assigned pd.2013 pd.markers
#> 1  At1g09210  predicted ER  predicted ER           ER      ER   ER lumen
#> 2  At1g21750  predicted ER  predicted ER           ER      ER   ER lumen
#> 3  At1g51760       unknown       unknown           ER      ER   ER lumen
#> 4  At1g56340  predicted ER  predicted ER           ER      ER   ER lumen
#> 5  At2g32920  predicted ER  predicted ER           ER      ER   ER lumen
#> 6  At2g47470  predicted ER  predicted ER           ER      ER   ER lumen
```
