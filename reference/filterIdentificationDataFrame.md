# Filter out unreliable PSMs.

A function to filter out PSMs matching to the decoy database, of rank
greater than one and matching non-proteotypic peptides.

## Usage

``` r
filterIdentificationDataFrame(
  x,
  decoy = "isDecoy",
  rank = "rank",
  accession = "DatabaseAccess",
  spectrumID = "spectrumID",
  verbose = isMSnbaseVerbose()
)
```

## Arguments

- x:

  A `data.frame` containing PSMs.

- decoy:

  The column name defining whether entries match the decoy database.
  Default is `"isDecoy"`. The column should be a `logical` and only PSMs
  holding a `FALSE` are retained. Ignored is set to `NULL`.

- rank:

  The column name holding the rank of the PSM. Default is `"rank"`. This
  column should be a `numeric` and only PSMs having rank equal to 1 are
  retained. Ignored is set to `NULL`.

- accession:

  The column name holding the protein (groups) accession. Default is
  `"DatabaseAccess"`. Ignored is set to `NULL`.

- spectrumID:

  The name of the spectrum identifier column. Default is `spectrumID`.

- verbose:

  A `logical` verbosity flag. Default is to take
  [`isMSnbaseVerbose()`](https://lgatto.github.io/MSnbase/reference/MSnbaseOptions.md).

## Value

A new `data.frame` with filtered out peptides and with the same columns
as the input `x`.

## Details

The PSMs should be stored in a `data.frame` such as those produced by
[`readMzIdData()`](https://lgatto.github.io/MSnbase/reference/readMzIdData.md).
Note that this function should be called before calling the
[reduce](https://lgatto.github.io/MSnbase/reference/reduce-data.frame-method.md)
method on a PSM `data.frame`.

## Author

Laurent Gatto
