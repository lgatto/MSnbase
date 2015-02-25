\name{calculateFragments-methods}
\alias{calculateFragments}
\alias{calculateFragments,character,missing-method}
\docType{methods}

\title{ Calculate ions produced by fragmentation. }

\description{
  These method calculates a-, b-, c-, x-, y- and z-ions produced by
  fragmentation.
}

\arguments{
  \item{sequence}{ \code{character}, peptide sequence. }
  \item{object}{ Object of class \code{"\linkS4class{Spectrum2}"} or
    \code{"missing"} . }
  \item{tolerance}{ \code{numeric} tolerance between the theoretical and
    measured MZ values (only available if \code{object} is not \code{missing}).}
  \item{method}{ \code{method} used for for duplicated matches. Choose
    \code{"highest"} or \code{"closest"} to select the peak with the highest
    intensity respectively the closest MZ in the tolerance range. If
    \code{"all"} is given all possible matches in the tolerance range are
    reported (only available if \code{object} is not \code{missing}).}
  \item{type}{ \code{character} vector of target ions; possible values:
    \code{c("a", "b", "c", "x", "y", "z")}; default: \code{type=c("b", "y")}. }
  \item{z}{ \code{numeric} desired charge state; default \code{z=1}. }
  \item{modifications}{ named \code{numeric} vector of used modifications. The
    name must correspond to the one-letter-code of the modified amino acid and
    the \code{numeric} value must represent the replacement mass, default:
    Carbamidomethyl \code{modifications=c(C=160.030649)}. }
  \item{neutralLoss}{ \code{logical} if \code{TRUE} (default) neutral loss is
    also calculated. Currently water and ammonia loss (shown with an \code{_}
    respectively a \code{*} in the results). }
  \item{verbose}{ \code{logical} if \code{TRUE} (default) the used
    modifications are printed. }
}

\section{Methods}{
  \describe{
    \item{\code{signature(sequence = "character", object = "missing", \dots)}}{
      Calculates the theoretical fragments for a peptide \code{sequence}.
      Returns a \code{data.frame} with the columns \code{c("mz", "ion", "type",
        "pos", "z", "seq")}.
    }
    \item{\code{signature(sequence = "character", object = "Spectrum2", \dots)}}{
      Calculates and matches the theoretical fragments for a peptide
      \code{sequence} and a \code{"\linkS4class{Spectrum2}"} \code{object}.
      The \code{\dots} arguments are passed to the internal functions.
      Currently \code{tolerance}, \code{method} and \code{relative} are
      supported.
      You could change the \code{tolerance} (default \code{25e-6}) and
      decide whether this tolerance should be applied relative
      (default \code{relative = TRUE}) or absolute (\code{relative = FALSE})
      to match the theoretical fragment MZ with the MZ of the spectrum.
      In cases of multiple matches use \code{method} to select the peak with
      the highest intensity (\code{method = "highest"}, default) respectively
      closest MZ (\code{method = "closes"}). If \code{method = "all"} is set
      all possible matches in the current tolerance range are reported.
      Returns the same \code{data.frame} as above but the \code{mz} column
      represents the matched MZ values of the spectrum. Additionally there
      is a column \code{error} that contains the difference between the observed
      MZ (from the spectrum) to the theoretical fragment MZ.
    }
  }
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\examples{
## find path to a mzXML file
file <- dir(system.file(package = "MSnbase", dir = "extdata"),
            full.name = TRUE, pattern = "mzXML$")

## create basic MSnExp
msexp <- readMSData(file)

## centroid them
msexp <- pickPeaks(msexp)

## calculate fragments for ACE with default modification
calculateFragments("ACE", modifications=c(C=160.030649))

## calculate fragments for ACE without any modifications
calculateFragments("ACE", modifications=NULL)

calculateFragments("VESITARHGEVLQLRPK",
                   type=c("a", "b", "c", "x", "y", "z"),
                   z=1:2)

calculateFragments("VESITARHGEVLQLRPK", msexp[[1]])
}

\keyword{methods}
