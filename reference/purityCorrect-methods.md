# Performs reporter ions purity correction

Manufacturers sometimes provide purity correction values indicating the
percentages of each reporter ion that have masses differing by +/- n Da
from the nominal reporter ion mass due to isotopic variants. This
correction is generally applied after reporter peaks quantitation.

Purity correction here is applied using `solve` from the `base` package
using the purity correction values as coefficient of the linear system
and the reporter quantities as the right-hand side of the linear system.
'NA' values are ignored and negative intensities after correction are
also set to 'NA'.

A more elaborated purity correction method is described in Shadforth *et
al.*, i-Tracker: for quantitative proteomics using iTRAQ. BMC Genomics.
2005 Oct 20;6:145. (PMID 16242023).

Function `makeImpuritiesMatrix(x, filename, edit = TRUE)` helps the user
to create such a matrix. The function can be used in two ways. If given
an integer `x`, it is used as the dimension of the square matrix (i.e
the number of reporter ions). For TMT6-plex and iTRAQ4-plex, default
values taken from manufacturer's certification sheets are used as
templates, but batch specific values should be used whenever possible.
Alternatively, the `filename` of a `csv` spreadsheet can be provided.
The sheet should define the correction factors as illustrated below
(including reporter names in the first column and header row) and the
corresponding correction matrix is calculated. Examples of such `csv`
files are available in the package's `extdata` directory. Use
`dir(system.file("extdata", package = "MSnbase"), pattern = "PurityCorrection", full.names = TRUE)`
to locate them. If `edit = TRUE`, the the matrix can be edited before it
is returned.

## Methods

- `signature(object = "MSnSet", impurities = "matrix")`:

## Arguments

- object:

  An object of class
  `"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`.

- impurities:

  A square 'matrix' of dim equal to ncol(object) defining the correction
  coefficients to be applied. The reporter ions should be ordered along
  the columns and the relative percentages along the rows.

  As an example, below is the correction factors as provided in an ABI
  iTRAQ 4-plex certificate of analysis:

  |          |         |         |         |         |
  |----------|---------|---------|---------|---------|
  | reporter | % of -2 | % of -1 | % of +1 | % of +2 |
  | 114      | 0.0     | 1.0     | 5.9     | 0.2     |
  | 115      | 0.0     | 2.0     | 5.6     | 0.1     |
  | 116      | 0.0     | 3.0     | 4.5     | 0.1     |
  | 117      | 0.1     | 4.0     | 3.5     | 0.1     |

  The impurity table will be

  |       |       |       |       |
  |-------|-------|-------|-------|
  | 0.929 | 0.059 | 0.002 | 0.000 |
  | 0.020 | 0.923 | 0.056 | 0.001 |
  | 0.000 | 0.030 | 0.924 | 0.045 |
  | 0.000 | 0.001 | 0.040 | 0.923 |

  where, the diagonal is computed as 100 - sum of rows of the original
  table and subsequent cells are directly filled in.

  Similarly, for TMT 6-plex tags, we observe

  |          |         |         |         |           |         |         |
  |----------|---------|---------|---------|-----------|---------|---------|
  | reporter | % of -3 | % of -2 | % of -1 | % of +1 % | % of +2 | % of +3 |
  | 126      | 0       | 0       | 0       | 6.1       | 0       | 0       |
  | 127      | 0       | 0       | 0.5     | 6.7       | 0       | 0       |
  | 128      | 0       | 0       | 1.1     | 4.2       | 0       | 0       |
  | 129      | 0       | 0       | 1.7     | 4.1       | 0       | 0       |
  | 130      | 0       | 0       | 1.6     | 2.1       | 0       | 0       |
  | 131      | 0       | 0.2     | 3.2     | 2.8       | 0       | 0       |

  and obtain the following impurity correction matrix

  |       |       |       |       |       |       |
  |-------|-------|-------|-------|-------|-------|
  | 0.939 | 0.061 | 0.000 | 0.000 | 0.000 | 0.000 |
  | 0.005 | 0.928 | 0.067 | 0.000 | 0.000 | 0.000 |
  | 0.000 | 0.011 | 0.947 | 0.042 | 0.000 | 0.000 |
  | 0.000 | 0.000 | 0.017 | 0.942 | 0.041 | 0.000 |
  | 0.000 | 0.000 | 0.000 | 0.016 | 0.963 | 0.021 |
  | 0.000 | 0.000 | 0.000 | 0.002 | 0.032 | 0.938 |

  For iTRAQ 8-plex, given the following correction factors (to make such
  a matrix square, if suffices to add -4, -3, +3 and +4 columns filled
  with zeros):

  |     |     |     |     |     |
  |-----|-----|-----|-----|-----|
  | TAG | -2  | -1  | +1  | +2  |
  | 113 | 0   | 2.5 | 3   | 0.1 |
  | 114 | 0   | 1   | 5.9 | 0.2 |
  | 115 | 0   | 2   | 5.6 | 0.1 |
  | 116 | 0   | 3   | 4.5 | 0.1 |
  | 117 | 0.1 | 4   | 3.5 | 0.1 |
  | 118 | 0.1 | 2   | 3   | 0.1 |
  | 119 | 0.1 | 2   | 4   | 0.1 |
  | 121 | 0.1 | 2   | 3   | 0.1 |

  we calculate the impurity correction matrix shown below

  |                |       |       |       |       |       |       |       |       |
  |----------------|-------|-------|-------|-------|-------|-------|-------|-------|
  |                | 113   | 114   | 115   | 116   | 117   | 118   | 119   | 121   |
  | % reporter 113 | 0.944 | 0.030 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
  | % reporter 114 | 0.010 | 0.929 | 0.059 | 0.002 | 0.000 | 0.000 | 0.000 | 0.000 |
  | % reporter 115 | 0.000 | 0.020 | 0.923 | 0.056 | 0.001 | 0.000 | 0.000 | 0.000 |
  | % reporter 116 | 0.000 | 0.000 | 0.030 | 0.924 | 0.045 | 0.001 | 0.000 | 0.000 |
  | % reporter 117 | 0.000 | 0.000 | 0.001 | 0.040 | 0.923 | 0.035 | 0.001 | 0.000 |
  | % reporter 118 | 0.000 | 0.000 | 0.000 | 0.001 | 0.020 | 0.948 | 0.030 | 0.001 |
  | % reporter 119 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.020 | 0.938 | 0.040 |
  | % reporter 121 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.020 | 0.948 |

  Finally, for a TMT 10-plex impurity matrix (for example lot
  [RH239932](https://www.thermofisher.com/document-connect/document-connect.html?url=https://assets.thermofisher.com/TFS-Assets/LSG/certificate/Certificates%20of%20Analysis/RH239932%5C_90309.pdf))

  |  |  |  |  |  |  |  |  |  |  |
  |----|----|----|----|----|----|----|----|----|----|
  | . | -2 | -1 | 1 | 2 | 126 | 0.0 | 0.0 | 5.0 (127C) | 0.0 (128C) |
  | 127N | 0.0 | 0.2 | 5.8 (128N) | 0.0 (129N) | 127C | 0.0 | 0.3 (126) | 4.8 (128C) | 0.0 (129C) |
  | 128N | 0.0 | 0.4 (127N) | 4.1 (129N) | 0.0 (130N) | 128C | 0.0 (126) | 0.6 (127C) | 3.0 (129C) | 0.0 (130C) |
  | 129N | 0.0 (127N) | 0.8 (128N) | 3.5 (130N) | 0.0 (131) | 129C | 0.0 (127C) | 1.4 (128C) | 2.4 (130C) | 0.0 |
  | 130N | 0.1 (128N) | 1.5 (129N) | 2.4 (131) | 3.2 | 130C | 0.0 (128C) | 1.7 (129C) | 1.8 | 0.0 |

  (Note that a previous example, taken from lot
  [PB199188A](https://www.thermofisher.com/document-connect/document-connect.html?url=https://assets.thermofisher.com/TFS-Assets/LSG/certificate/Certificates%20of%20Analysis/PB199188A%5C_90110.pdf),
  contained a typo.)

  the impurity correction matrix is

  |                 |       |       |       |       |       |       |       |       |       |       |
  |-----------------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
  | .               | 126   | 127N  | 127C  | 128N  | 128C  | 129N  | 129C  | 130N  | 130C  | 131   |
  | % reporter 126  | 0.950 | 0.000 | 0.050 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
  | % reporter 127N | 0.000 | 0.940 | 0.000 | 0.058 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
  | % reporter 127C | 0.003 | 0.000 | 0.949 | 0.000 | 0.048 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
  | % reporter 128N | 0.000 | 0.004 | 0.000 | 0.955 | 0.000 | 0.041 | 0.000 | 0.000 | 0.000 | 0.000 |
  | % reporter 128C | 0.000 | 0.000 | 0.006 | 0.000 | 0.964 | 0.000 | 0.030 | 0.000 | 0.000 | 0.000 |
  | % reporter 129N | 0.000 | 0.000 | 0.000 | 0.008 | 0.000 | 0.957 | 0.000 | 0.035 | 0.000 | 0.000 |
  | % reporter 129C | 0.000 | 0.000 | 0.000 | 0.000 | 0.014 | 0.000 | 0.962 | 0.000 | 0.024 | 0.000 |
  | % reporter 130N | 0.000 | 0.000 | 0.000 | 0.001 | 0.000 | 0.015 | 0.000 | 0.928 | 0.000 | 0.024 |
  | % reporter 130C | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.017 | 0.000 | 0.965 | 0.000 |
  | % reporter 131  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.002 | 0.000 | 0.020 | 0.000 | 0.956 |

  These examples are provided as defaults impurity correction matrices
  in `makeImpuritiesMatrix`.

## Examples

``` r
## quantifying full experiment
data(msnset)
impurities <- matrix(c(0.929,0.059,0.002,0.000,
           0.020,0.923,0.056,0.001,
           0.000,0.030,0.924,0.045,
           0.000,0.001,0.040,0.923),
         nrow=4, byrow = TRUE)
## or, using makeImpuritiesMatrix()
if (FALSE) impurities <- makeImpuritiesMatrix(4) # \dontrun{}
msnset.crct <- purityCorrect(msnset, impurities)
head(exprs(msnset))
#>     iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> X1   1347.6158  2247.3097  3927.6931  7661.1463
#> X10   739.9861   799.3501   712.5983   940.6793
#> X11 27638.3582 33394.0252 32104.2879 26628.7278
#> X12 31892.8928 33634.6980 37674.7272 37227.7119
#> X13 26143.7542 29677.4781 29089.0593 27902.5608
#> X14  6448.0829  6234.1957  6902.8903  6437.2303
head(exprs(msnset.crct))
#>     iTRAQ4.114 iTRAQ4.115 iTRAQ4.116 iTRAQ4.117
#> X1   1402.9442  2214.0346  3762.2549  8114.4429
#> X10   779.4666   793.0792   678.8083   985.2003
#> X11 29034.3781 33271.0470 31484.7131 27279.1383
#> X12 33618.9092 33046.3075 37031.6133 38492.1376
#> X13 27508.0038 29440.9296 28390.4561 28814.2463
#> X14  6809.7600  6090.7894  6799.5030  6636.1450
processingData(msnset.crct)
#> - - - Processing information - - -
#> Data loaded: Wed May 11 18:54:39 2011 
#> iTRAQ4 quantification by trapezoidation: Wed Apr  1 21:41:53 2015 
#> Purity corrected: Fri Apr 10 14:44:58 2026 
#>  MSnbase version: 1.1.22 

## default impurity matrix for iTRAQ 8-plex
makeImpuritiesMatrix(8, edit = FALSE)
#>                  113   114   115   116   117   118   119   121
#> % reporter 113 0.944 0.030 0.001 0.000 0.000 0.000 0.000 0.000
#> % reporter 114 0.010 0.929 0.059 0.002 0.000 0.000 0.000 0.000
#> % reporter 115 0.000 0.020 0.923 0.056 0.001 0.000 0.000 0.000
#> % reporter 116 0.000 0.000 0.030 0.924 0.045 0.001 0.000 0.000
#> % reporter 117 0.000 0.000 0.001 0.040 0.923 0.035 0.001 0.000
#> % reporter 118 0.000 0.000 0.000 0.001 0.020 0.948 0.030 0.001
#> % reporter 119 0.000 0.000 0.000 0.000 0.001 0.020 0.938 0.040
#> % reporter 121 0.000 0.000 0.000 0.000 0.000 0.001 0.020 0.948

## default impurity matrix for TMT 10-plex
makeImpuritiesMatrix(10, edit = FALSE)
#>                   126  127N  127C  128N  128C  129N  129C  130N  130C   131
#> % reporter 126  0.950 0.000 0.050 0.000 0.000 0.000 0.000 0.000 0.000 0.000
#> % reporter 127N 0.000 0.940 0.000 0.058 0.000 0.000 0.000 0.000 0.000 0.000
#> % reporter 127C 0.003 0.000 0.949 0.000 0.048 0.000 0.000 0.000 0.000 0.000
#> % reporter 128N 0.000 0.004 0.000 0.955 0.000 0.041 0.000 0.000 0.000 0.000
#> % reporter 128C 0.000 0.000 0.006 0.000 0.964 0.000 0.030 0.000 0.000 0.000
#> % reporter 129N 0.000 0.000 0.000 0.008 0.000 0.957 0.000 0.035 0.000 0.000
#> % reporter 129C 0.000 0.000 0.000 0.000 0.014 0.000 0.962 0.000 0.024 0.000
#> % reporter 130N 0.000 0.000 0.000 0.001 0.000 0.015 0.000 0.928 0.000 0.024
#> % reporter 130C 0.000 0.000 0.000 0.000 0.000 0.000 0.017 0.000 0.965 0.000
#> % reporter 131  0.000 0.000 0.000 0.000 0.000 0.002 0.000 0.020 0.000 0.956
```
