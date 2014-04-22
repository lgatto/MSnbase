\name{MSnExp-class}
\Rdversion{1.1}
\docType{class}
\concept{MSnExp}

\alias{MSnExp-class}
\alias{class:MSnExp}
\alias{MSnExp}

\alias{bin,MSnExp-method}
\alias{clean,MSnExp-method}
\alias{compareSpectra,MSnExp,missing-method}
\alias{extractPrecSpectra,MSnExp-method}
\alias{extractPrecSpectra,MSnExp,numeric-method}
%% \alias{extractSpectra,MSnExp-method}
%% \alias{extractSpectra,MSnExp,logical-method}
\alias{pickPeaks,MSnExp-method}
\alias{plot,MSnExp}
\alias{plot,MSnExp,missing-method}
\alias{quantify,MSnExp-method}
\alias{quantify,MSnExp,character-method}
\alias{removePeaks,MSnExp-method}
\alias{removeReporters,MSnExp-method}
\alias{smooth,MSnExp-method}
\alias{addIdentificationData,MSnExp-method}
\alias{removeNoId,MSnExp-method}
\alias{removeMultipleAssignment,MSnExp-method}
\alias{idSummary,MSnExp-method}
\alias{show,MSnExp-method}
\alias{spectra,MSnExp-method}
\alias{trimMz,MSnExp-method}
\alias{trimMz,MSnExp,numeric-method}

\title{The 'MSnExp' Class for MS Data And Meta-Data}

\description{
  The \code{MSnExp} class encapsulates data and meta-data for mass
  spectrometry experiments, as described in the \code{slots}
  section. Several data files (currently in \code{mzXML}) can be loaded
  together with the function \code{\link{readMSData}}.

  This class extends the virtual \code{"\linkS4class{pSet}"} class.
}


\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("MSnExp",...)}. However, it is preferred to use the
  \code{\link{readMSData}} function that will read raw mass
  spectrometry data to generate a valid \code{"MSnExp"} instance.
}

\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"environment"}
      containing the MS spectra (see \code{"\linkS4class{Spectrum1}"}
      and \code{"\linkS4class{Spectrum2}"}).
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{phenoData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      experimenter-supplied variables describing sample (i.e the
      individual tags for an labelled MS experiment)
      See \code{\link{phenoData}} for more details.
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{featureData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing variables
      describing features (spectra in our case), e.g. identificaiton data,
      peptide sequence, identification score,... (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{featureData}} for
      more details.
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAPE}"}, containing details of experimental
      methods. See \code{\link{experimentData}} for more details.
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{protocolData}} for
      more details.
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{processingData}:}{Object of class
      \code{"\linkS4class{MSnProcess}"} that records all processing.
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of R,
      the Biobase package, \code{"\linkS4class{pSet}"} and
      \code{MSnExp} of the current instance.
      Slot is inherited from \code{"\linkS4class{pSet}"}.
      Intended for developer use and debugging (inherited from
      \code{"\linkS4class{eSet}"}). }
  }
}

\section{Extends}{
Class \code{"\linkS4class{pSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "pSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "pSet", distance 3.
}

\section{Methods}{
  See the \code{"\linkS4class{pSet}"} class for documentation on
  accessors inherited from \code{pSet}, subsetting and general attribute
  accession.

  \describe{
    \item{bin}{\code{signature(object = "MSnExp")}: Bins spectra.
      See \code{\link{bin}} documentation for more details and examples. }
    \item{clean}{\code{signature(object = "MSnExp")}: Removes unused 0
      intensity data points. See \code{\link{clean}} documentation
      for more details and examples. }
    \item{compareSpectra}{\code{signature(object1 = "Spectrum",
      object2 = "missing")}: Compares spectra.  See
      \code{\link{compareSpectra}} documentation for more details and
      examples. }
    \item{extractPrecSpectra}{\code{signature(object = "MSnExp", prec =
	"numeric")}: extracts spectra with precursor MZ value equal to
      \code{prec} and returns an object of class 'MSnExp'.
      See \code{\link{extractPrecSpectra}} documentation for more
      details and examples. }
    % \item{extractSpectra}{\code{signature(object = "MSnExp", selected =
    % 	"logical")}: extracts spectra corresponding to 'TRUE' in
    %   \code{selected}.
    % See \code{\link{extractSpectra}} documentation for more details
    % and examples. }
    \item{pickPeaks}{\code{signature(object = "MSnExp")}: Performs the peak
      picking to generate centroided spectra.
      See \code{\link{pickPeaks}} documentation for more
      details and examples. }
    \item{plot}{\code{signature(x = "MSnExp", y = "missing")}: Plots all
      the spectra of the \code{MSnExp} instance. See
      \code{\link{plot.MSnExp}} documentation for more details. }
    \item{plot2d}{\code{signature(object = "MSnExp", ...)}:
      Plots retention time against precursor MZ for \code{MSnExp}
      instances. See \code{\link{plot2d}} documentation for more
      details. }
    \item{plotDensity}{\code{signature(object = "MSnExp", ...)}:
      Plots the density of parameters of interest.
      instances. See \code{\link{plotDensity}} documentation for more
      details. }
    \item{plotMzDelta}{\code{signature(object = "MSnExp", ...)}:
      Plots a histogram of the m/z difference betwee all of the highest
      peaks of all MS2 spectra of an experiment.
      See \code{\link{plotMzDelta}} documentation for more details. }
    \item{quantify}{\code{signature(object = "MSnExp")}: Performs
      quantification for all the MS2 spectra of the \code{MSnExp}
      instance. See \code{\link{quantify}} documentation for more
      details. }
    \item{removePeaks}{\code{signature(object = "MSnExp")}: Removes
      peaks lower that a threshold \code{t}. See
      \code{\link{removePeaks}} documentation for more details and
      examples. }
    \item{removeReporters}{\code{signature(object = "MSnExp", ...)}:
      Removes reporter ion peaks from all MS2 spectra of an
      experiment. See \code{\link{removeReporters}} documentation for
      more details and examples. }
    \item{smooth}{\code{signature(x = "MSnExp")}: Smooths spectra.
      See \code{\link{smooth}} documentation for more details and examples. }
    \item{addIdentificationData}{\code{signature(object = "MSnExp", ...)}:
      Adds identification data to an experiment.
      See \code{\link{addIdentificationData}} documentation for
      more details and examples. }

    \item{removeNoId}{\code{signature(object = "MSnExp", fcol =
	"pepseq", keep = NULL)}: Removes non-identified features. See
      \code{\link{removeNoId}} documentation for more details and
      examples. }
    
    \item{removeMultipleAssignment}{\code{signature(object = "MSnExp",
	fcol = "nprot")}: Removes protein groups with more than one
	member. The latter is defined by extracting a feature variable
	(default is \code{"nprot"}). }

    \item{idSummary}{\code{signature(object = "MSnExp", ...)}:
      Prints a summary that lists the percentage of identified features per
      file (called \code{coverage}).}

    \item{show}{\code{signature(object = "MSnExp")}: Displays object
      content as text. }

    \item{trimMz}{\code{signature(object = "MSnExp")}: Trims the MZ
      range of all the spectra of the \code{MSnExp} instance. See
      \code{\link{trimMz}} documentation for more details and
      examples. }
  }
}

\references{
  Information about the mzXML format as well converters from vendor
  specific formats to mzXML:
  \url{http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML}.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{"\linkS4class{pSet}"} and \code{\link{readMSData}} for loading
  \code{mzXML}, \code{mzData} or \code{mzML} files to generate an
  instance of \code{MSnExp}.
}

\examples{
mzxmlfile <- dir(system.file("extdata",package="MSnbase"),
                 pattern="mzXML",full.names=TRUE)
msnexp <- readMSData(mzxmlfile)
msnexp
}

\keyword{classes}
