\name{addIdentificationData-methods}
\docType{methods}
\alias{addIdentificationData-methods}
\alias{addIdentificationData}

\title{ Adds Identification Data }

\description{
  This methods add identification data to an experiment
  \code{"\linkS4class{MSnExp}"} or to a \code{"\linkS4class{MSnSet}"}.
}

\section{Methods}{
  \describe{

    \item{\code{signature(object = "MSnExp", id = "character",
        fcol = "character", icol = "character", verbose = "logical")}}{
      Adds the identification data stored in mzIdentML files to a
      \code{"\linkS4class{MSnExp}"} instance. The method handles
      one or multiple mzIdentML files provided via \code{id}. \code{id} has to
      be a \code{character} vector of valid filenames.
      \code{fcol} and \code{icol} specify the columns in the \code{featureData}
      slot and the identification \code{data.frame} that are used to merge both
      data together as \code{character} vectors
      (defaults are \code{fcol = c("spectrum.file", "acquisition.number")} and
      \code{icol = c("spectrumFile", "acquisitionnum")}).
      The \code{verbose} argument (default is
      \code{TRUE}) defines whether status messages should be showed.
    }

    \item{\code{signature(object = "MSnExp", id = "mzID",
        fcol = "character", icol = "character", verbose = "logical")}}{
      Same as above but \code{id} could be an \code{mzID} object.
    }

    \item{\code{signature(object = "MSnExp", id = "mzIDCollection",
        fcol = "character", icol = "character", verbose = "logical")}}{
      Same as above but \code{id} could be an \code{mzIDCollection} object.
    }

    \item{\code{signature(object = "MSnExp", id = "data.frame",
        fcol = "character", icol = "character", verbose = "logical")}}{
      Same as above but \code{id} could be a \code{data.frame}.
    }

    \item{\code{signature(object = "MSnSet", id = "character",
	    verbose = "logical")}}{
      Adds the identification data stored in mzIdentML files to an
      \code{"\linkS4class{MSnSet}"} instance.
      The method handles one or multiple mzIdentML files provided via
      \code{id}. \code{id} has to be a \code{character} vector of valid
      filenames.

      \code{fcol} and \code{icol} specify the columns in the \code{featureData}
      slot and the identification \code{data.frame} that are used to merge both
      data together as \code{character} vectors
      (defaults are \code{fcol = c("spectrum.file", "acquisition.number")} and
      \code{icol = c("spectrumFile", "acquisitionnum")}).

      The \code{verbose} argument (default is
      \code{TRUE}) defines whether status messages should be showed.
    }

    \item{\code{signature(object = "MSnSet", id = "mzID",
        fcol = "character", icol = "character", verbose = "logical")}}{
      Same as above but \code{id} could be an \code{mzID} object.
    }

    \item{\code{signature(object = "MSnSet", id = "mzIDCollection",
        fcol = "character", icol = "character", verbose = "logical")}}{
      Same as above but \code{id} could be an \code{mzIDCollection} object.
    }

    \item{\code{signature(object = "MSnSet", id = "data.frame",
        fcol = "character", icol = "character", verbose = "logical")}}{
      Same as above but \code{id} could be a \code{data.frame}.
    }
  }
}

\details{
  The \code{featureData} slots in an \code{"\linkS4class{MSnExp}"} or an
  \code{"\linkS4class{MSnSet}"} instance provides only one row per MS2
  spectrum but the identification is not always bijective. If multiple
  possible matches are present only the highest ranked identification is
  added.

  The column \code{nprot} contains the number of members in the protein
  group; the columns \code{accession} and \code{description} contain a
  semicolon separated list of all matches sorted by their \code{rank}
  values. The columns \code{npsm.prot} and \code{npep.prot} represent
  the number of PSMs and peptides that were matched to a particular
  protein group. The column \code{npsm.pep} indicates how many PSMs were
  attributed to a peptide (as defined by its sequence \code{pepseq}).
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  \code{\linkS4class{MSnExp}} and \code{\link{MSnSet}}.
}

\examples{
## find path to a mzXML file
quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
## find path to a mzIdentML file
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")

## create basic MSnExp
msexp <- readMSData(quantFile)

## add identification information
msexp <- addIdentificationData(msexp, identFile)

## access featureData; please note the multiple identification data
## for spectrum 1 (row 1)
fData(msexp)

idSummary(msexp)
}

\keyword{methods}
