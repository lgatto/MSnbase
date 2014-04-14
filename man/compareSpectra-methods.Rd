\name{compareSpectra-methods}
\docType{methods}
\alias{compareSpectra-methods}
% \alias{compareSpectra,MSnExp-method}
% \alias{compareSpectra,Spectrum-method}
\alias{compareSpectra}


\title{ Compare Spectra of an 'MSnExp' or 'Spectrum' instances }

\description{
  This method compares spectra (\code{Spectrum} instances) pairwise
  or all spectra of an experiment (\code{MSnExp} instances). Currently
  the comparison is based on the number of common peaks \code{fun = "common"},
  the Pearson correlation \code{fun = "cor"} or the dot product
  \code{fun = "dotproduct"}. The user could use his own function.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object1 = "MSnExp", object2 = "missing",
      fun = "character", \dots)}}{
      Compares all spectra in an \code{MSnExp} object.
      Use \code{fun} to compare spectra by the number of common peaks
      \code{fun = "common"}, by the pearson correlation \code{fun = "cor"} and
      by the dot product \code{fun = "dotproduct"}. The \code{\dots} arguments
      are passed to the internal functions. For \code{fun = "common"} you could
      change the \code{tolerance} (default \code{25e-6}) and decide whether this
      tolerance should be applied relative (default \code{relative = TRUE}) or
      absolute (\code{relative = FALSE}).
      Instead of the three predefined functions for \code{fun} you could add
      your own comparison function. Therefore you have to define a function that
      takes two \code{\linkS4class{Spectrum}} objects as the first two arguments
      and \code{\dots} as the third. The function must return a single
      \code{numeric} value. Please see the example section.
      Before \code{fun = "cor"} and \code{fun = "dotproduct"} could be used the
      spectra need to be binned. Use the \code{binSize} argument (in Dalton) to
      control the binning precision. Please see \code{\link{bin}} for details.
      Returns a \code{n*n} matrix.
    }
    \item{\code{signature(object1 = "Spectrum", object2 = "Spectrum",
        fun = "character", \dots)}}{
      Compares two \code{Spectrum} objects. Please see the above
      explanation for \code{fun} and \code{\dots}.
      Returns a single \code{numeric} value.
    }
  }
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  \code{\link{bin}}, \code{\link{clean}}, \code{\link{pickPeaks}},
  \code{\link{smooth}}, \code{\link{removePeaks}} and \code{\link{trimMz}}
  for other spectra processing methods.
}

\references{
Stein, S. E., & Scott, D. R. (1994).
Optimization and testing of mass spectral library search algorithms for
compound identification.
Journal of the American Society for Mass Spectrometry, 5(9), 859-866.
doi: http://dx.doi.org/10.1016/1044-0305(94)87009-8

Lam, H., Deutsch, E. W., Eddes, J. S., Eng, J. K., King, N., Stein, S. E.
and Aebersold, R. (2007)
Development and validation of a spectral library searching method for peptide
identification from MS/MS.
Proteomics, 7: 655-667.
doi: http://dx.doi.org/10.1002/pmic.200600625
}

\examples{
s1 <- new("Spectrum2", mz=1:10, intensity=1:10)
s2 <- new("Spectrum2", mz=1:10, intensity=10:1)
compareSpectra(s1, s2)
compareSpectra(s1, s2, fun="cor", binSize=2)
compareSpectra(s1, s2, fun="dotproduct")

## define our own (useless) comparison function (it is just a basic example)
equalLength <- function(x, y, ...) {
  return(peaksCount(x)/(peaksCount(y)+.Machine$double.eps))
}
compareSpectra(s1, s2, fun=equalLength)
compareSpectra(s1, new("Spectrum2", mz=1:5, intensity=1:5), fun=equalLength)
compareSpectra(s1, new("Spectrum2"), fun=equalLength)

data(itraqdata)
compareSpectra(itraqdata, fun="cor")
}


\keyword{methods}
