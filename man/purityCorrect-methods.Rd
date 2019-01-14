\name{purityCorrect-methods}
\docType{methods}

\alias{purityCorrect-methods}
\alias{purityCorrect}
\alias{makeImpuritiesMatrix}

\title{ Performs reporter ions purity correction }

\description{
  Manufacturers sometimes provide purity correction values indicating
  the percentages of each reporter ion that have masses differing by
  +/- n Da from the nominal reporter ion mass due to isotopic
  variants. This correction is generally applied after reporter peaks
  quantitation.

  Purity correction here is applied using \code{solve} from the
  \code{base} package using the purity correction values as coefficient of
  the linear system and the reporter quantities as the right-hand side
  of the linear system. 'NA' values are ignored and negative
  intensities after correction are also set to 'NA'.

  A more elaborated purity correction method is described in
  Shadforth \emph{et al.}, i-Tracker: for quantitative proteomics using
  iTRAQ. BMC Genomics. 2005 Oct 20;6:145. (PMID 16242023).

  Function \code{makeImpuritiesMatrix(x, filename, edit = TRUE)} helps
  the user to create such a matrix. The function can be used in two ways.
  If given an integer \code{x}, it is used as the dimension of the
  square matrix (i.e the number of reporter ions). For TMT6-plex and
  iTRAQ4-plex, default values taken from manufacturer's certification
  sheets are used as templates, but batch specific values should be used
  whenever possible. Alternatively, the \code{filename} of a \code{csv}
  spreadsheet can be provided. The sheet should define the correction
  factors as illustrated below (including reporter names in the first
  column and header row) and the corresponding correction matrix is
  calculated. Examples of such \code{csv} files are available in the
  package's \code{extdata} directory. Use
  \code{dir(system.file("extdata", package = "MSnbase"), pattern =
    "PurityCorrection", full.names = TRUE)} to locate them.
  If \code{edit = TRUE}, the the matrix can be edited before
  it is returned.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnSet", impurities = "matrix")}}{
    }
  }
}

\arguments{
  \item{object}{ An object of class \code{"\linkS4class{MSnSet}"}. }
  \item{impurities}{ A square 'matrix' of dim equal to ncol(object)
    defining the correction coefficients to be applied.
    The reporter ions should be ordered along the columns
    and the relative percentages along the rows.

    As an example, below is the correction factors as provided in an ABI
    iTRAQ 4-plex certificate of analysis:
    \tabular{lrrrr}{
      reporter \tab \% of -2 \tab \% of -1 \tab \% of +1 \tab \% of +2 \cr
      114 \tab 0.0 \tab 1.0 \tab 5.9 \tab 0.2 \cr
      115 \tab 0.0 \tab 2.0 \tab 5.6 \tab 0.1 \cr
      116 \tab 0.0 \tab 3.0 \tab 4.5 \tab 0.1 \cr
      117 \tab 0.1 \tab 4.0 \tab 3.5 \tab 0.1 \cr
    }

    The impurity table will be
    \tabular{rrrr}{
      0.929 \tab 0.059 \tab 0.002 \tab 0.000 \cr
      0.020 \tab 0.923 \tab 0.056 \tab 0.001 \cr
      0.000 \tab 0.030 \tab 0.924 \tab 0.045 \cr
      0.000 \tab 0.001 \tab 0.040 \tab 0.923 \cr
    }
    where, the diagonal is computed as 100 - sum of rows of the
    original table and subsequent cells are directly filled in.

    Similarly, for TMT 6-plex tags, we observe
    \tabular{lrrrrrr}{
      reporter \tab \% of -3 \tab \% of -2 \tab \% of -1 \tab \% of +1 \% \tab \% of +2 \tab \% of +3 \cr
      126 \tab 0 \tab 0 \tab 0 \tab 6.1 \tab 0 \tab 0 \cr
      127 \tab 0 \tab 0 \tab 0.5 \tab 6.7 \tab 0 \tab 0 \cr
      128 \tab 0 \tab 0 \tab 1.1 \tab 4.2 \tab 0 \tab 0 \cr
      129 \tab 0 \tab 0 \tab 1.7 \tab 4.1 \tab 0 \tab 0 \cr
      130 \tab 0 \tab 0 \tab 1.6 \tab 2.1 \tab 0 \tab 0 \cr
      131 \tab 0 \tab 0.2 \tab 3.2 \tab 2.8 \tab 0 \tab 0 \cr
    }
    and obtain the following impurity correction matrix
    \tabular{rrrrrr}{
      0.939 \tab 0.061 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      0.005 \tab 0.928 \tab 0.067 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      0.000 \tab 0.011 \tab 0.947 \tab 0.042 \tab 0.000 \tab 0.000 \cr
      0.000 \tab 0.000 \tab 0.017 \tab 0.942 \tab 0.041 \tab 0.000 \cr
      0.000 \tab 0.000 \tab 0.000 \tab 0.016 \tab 0.963 \tab 0.021 \cr
      0.000 \tab 0.000 \tab 0.000 \tab 0.002 \tab 0.032 \tab 0.938 \cr
    }

    For iTRAQ 8-plex, given the following correction factors (to make
    such a matrix square, if suffices to add -4, -3, +3 and +4 columns
    filled with zeros):

    \tabular{lllll}{
      TAG \tab -2 \tab -1 \tab +1 \tab +2 \cr
      113 \tab	0 \tab	2.5 \tab 3 \tab	0.1 \cr
      114 \tab	0 \tab	1 \tab 5.9 \tab	0.2 \cr
      115 \tab	0 \tab	2 \tab 5.6 \tab	0.1 \cr
      116 \tab	0 \tab	3 \tab 4.5 \tab	0.1 \cr
      117 \tab	0.1 \tab 4 \tab	 3.5 \tab 0.1 \cr
      118 \tab	0.1 \tab 2 \tab	3 \tab 0.1 \cr
      119 \tab	0.1 \tab 2 \tab	4 \tab 0.1 \cr
      121 \tab	0.1 \tab 2 \tab	3 \tab 0.1 \cr
    }

    we calculate the impurity correction matrix shown below

    \tabular{lrrrrrrrr}{
		      \tab 113   \tab 114   \tab 115   \tab 116   \tab 117   \tab 118   \tab 119   \tab 121 \cr
      \% reporter 113 \tab 0.944 \tab 0.030 \tab 0.001 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 114 \tab 0.010 \tab 0.929 \tab 0.059 \tab 0.002 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 115 \tab 0.000 \tab 0.020 \tab 0.923 \tab 0.056 \tab 0.001 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 116 \tab 0.000 \tab 0.000 \tab 0.030 \tab 0.924 \tab 0.045 \tab 0.001 \tab 0.000 \tab 0.000 \cr
      \% reporter 117 \tab 0.000 \tab 0.000 \tab 0.001 \tab 0.040 \tab 0.923 \tab 0.035 \tab 0.001 \tab 0.000 \cr
      \% reporter 118 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.001 \tab 0.020 \tab 0.948 \tab 0.030 \tab 0.001 \cr
      \% reporter 119 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.001 \tab 0.020 \tab 0.938 \tab 0.040 \cr
      \% reporter 121 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.001 \tab 0.020 \tab 0.948 \cr
    }


    Finally, for a TMT 10-plex impurity matrix (for example lot \href{https://www.thermofisher.com/document-connect/document-connect.html?url=https://assets.thermofisher.com/TFS-Assets/LSG/certificate/Certificates\%20of\%20Analysis/RH239932\_90309.pdf}{RH239932})

    \tabular{llllllllll}{
      .    \tab   -2       \tab -1         \tab   1        \tab   2        \cr
      126  \tab 0.0        \tab 0.0        \tab 5.0 (127C) \tab 0.0 (128C) \cr
      127N \tab 0.0        \tab 0.2        \tab 5.8 (128N) \tab 0.0 (129N) \cr
      127C \tab 0.0        \tab 0.3 (126)  \tab 4.8 (128C) \tab 0.0 (129C) \cr
      128N \tab 0.0        \tab 0.4 (127N) \tab 4.1 (129N) \tab 0.0 (130N) \cr
      128C \tab 0.0 (126)  \tab 0.6 (127C) \tab 3.0 (129C) \tab 0.0 (130C) \cr
      129N \tab 0.0 (127N) \tab 0.8 (128N) \tab 3.5 (130N) \tab 0.0 (131) \cr
      129C \tab 0.0 (127C) \tab 1.4 (128C) \tab 2.4 (130C) \tab 0.0 \cr
      130N \tab 0.1 (128N) \tab 1.5 (129N) \tab 2.4 (131)  \tab 3.2 \cr
      130C \tab 0.0 (128C) \tab 1.7 (129C) \tab 1.8        \tab 0.0 \cr
      131  \tab 0.2 (129N) \tab 2.0 (130N) \tab 2.2        \tab 0.0 \cr
    }

    (Note that a previous example, taken from lot
    \href{https://www.thermofisher.com/document-connect/document-connect.html?url=https://assets.thermofisher.com/TFS-Assets/LSG/certificate/Certificates\%20of\%20Analysis/PB199188A\_90110.pdf}{PB199188A},
    contained a typo.)

    the impurity correction matrix is

    \tabular{lllllllllll}{

      .      \tab 126  \tab 127N  \tab 127C  \tab 128N  \tab 128C  \tab 129N \tab 129C \tab 130N \tab 130C \tab 131 \cr
      \% reporter 126  \tab 0.950 \tab 0.000 \tab 0.050 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 127N \tab 0.000 \tab 0.940 \tab 0.000 \tab 0.058 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 127C \tab 0.003 \tab 0.000 \tab 0.949 \tab 0.000 \tab 0.048 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 128N \tab 0.000 \tab 0.004 \tab 0.000 \tab 0.955 \tab 0.000 \tab 0.041 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 128C \tab 0.000 \tab 0.000 \tab 0.006 \tab 0.000 \tab 0.964 \tab 0.000 \tab 0.030 \tab 0.000 \tab 0.000 \tab 0.000 \cr
      \% reporter 129N \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.008 \tab 0.000 \tab 0.957 \tab 0.000 \tab 0.035 \tab 0.000 \tab 0.000 \cr
      \% reporter 129C \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.014 \tab 0.000 \tab 0.962 \tab 0.000 \tab 0.024 \tab 0.000 \cr
      \% reporter 130N \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.001 \tab 0.000 \tab 0.015 \tab 0.000 \tab 0.928 \tab 0.000 \tab 0.024 \cr
      \% reporter 130C \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.017 \tab 0.000 \tab 0.965 \tab 0.000 \cr
      \% reporter 131  \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.000 \tab 0.002 \tab 0.000 \tab 0.020 \tab 0.000 \tab 0.956 \cr
    }


    These examples are provided as defaults impurity correction matrices
    in \code{makeImpuritiesMatrix}.

  }
}

\examples{
## quantifying full experiment
data(msnset)
impurities <- matrix(c(0.929,0.059,0.002,0.000,
		       0.020,0.923,0.056,0.001,
		       0.000,0.030,0.924,0.045,
		       0.000,0.001,0.040,0.923),
		     nrow=4, byrow = TRUE)
## or, using makeImpuritiesMatrix()
\dontrun{impurities <- makeImpuritiesMatrix(4)}
msnset.crct <- purityCorrect(msnset, impurities)
head(exprs(msnset))
head(exprs(msnset.crct))
processingData(msnset.crct)

## default impurity matrix for iTRAQ 8-plex
makeImpuritiesMatrix(8, edit = FALSE)

## default impurity matrix for TMT 10-plex
makeImpuritiesMatrix(10, edit = FALSE)
}

\keyword{methods}
