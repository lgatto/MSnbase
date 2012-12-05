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

  Function \code{makeImpuritiesMatrix(x, edit = TRUE)} opens a simple graphical
  interface to guide the user in the creation of such a matrix. The
  function takes the dimension of the square matrix (i.e the number of
  reporter ions) as argument. When available, default values taken from
  manufacturer's certification sheets are provided, but batch specific
  values should be used whenever possible. \code{makeImpuritiesMatrix}
  returns the (possibly interactively updated) matrix to be used with
  \code{purityCorrect}. If \code{edit = FALSE}, the edition of the
  matrix is bypassed. 
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
    
    These two examples above are provided as defaults impurity correction
    matrices in \code{makeImpuritiesMatrix}.
  }
}

\examples{
## quantifying full experiment
data(itraqdata)
msnset <- quantify(itraqdata, method = "trap", reporters = iTRAQ4)
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
}

\keyword{methods}
