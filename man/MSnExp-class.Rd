\name{MSnExp-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnExp-class}
\alias{class:MSnExp}
\alias{MSnExp}

\alias{show,MSnExp-method}
\alias{header,MSnExp-method}
\alias{header}
\alias{length,MSnExp-method}
\alias{spectra,MSnExp-method}
\alias{spectra}
\alias{spectra<-,MSnExp-method}
\alias{spectra<-}
\alias{[,MSnExp-method}
\alias{precursorMz,MSnExp-method}
\alias{precursorCharge,MSnExp-method}
\alias{acquisitionNum,MSnExp-method}
\alias{rtime,MSnExp-method}
\alias{peaksCount,MSnExp-method}
\alias{msLevel,MSnExp-method}
\alias{collisionEnergy,MSnExp-method}
\alias{intensity,MSnExp-method}
\alias{mz,MSnExp-method}

\title{The "MSnExp" Class for MS Data And Meta-Data}

\description{
  The \code{MSnExp} class encapsulates data and meta-data for mass
  spectrometry experiments, as described in the \code{slots}
  section. Several data files (currently in \code{mzXML}) can be loaded
  together with the function \code{\link{readMzXMLData}}. 

  This class contains the virtual \code{"\linkS4class{eSet}"} class and
  adds MS-based proteomics specific attributes to describe the
  experiment, while making use of \code{"\linkS4class{eSet}"}'s existing
  slots for meta-data. 
}

\section{Slots}{
  \describe{
    \item{\code{spectra}:}{Object of class \code{"list"} of
      \code{"\linkS4class{Spectrum1}"} or \code{"\linkS4class{Spectrum2}"}
      objects - the data. }
    \item{\code{process}:}{Object of class
      \code{"\linkS4class{MSnProcess}"} that records all processing. }
    \item{\code{description}:}{Object of class \code{"character"}
      describing the experiment in free form. }
    \item{\code{fromFile}:}{Object of class \code{"numeric"} of same
      length as \code{spectra} indicating the file of origin of each
      spectrum. }
    \item{\code{files}:}{Object of class \code{"character"} recording
      the name of the \code{mzXML} files that have been loaded by
      \code{\link{readMzXMLData}} to create the \code{MSnExp} instance. }
    \item{\code{featureData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing variables
      describing features (spectra in our case), e.g. identificaiton data,
      peptide sequence, identification score,... (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{featureData}} for
      more details. }
    \item{\code{phenoData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing 
      experimenter-supplied variables describing sample (i.e the
      individual tags for an labelled MS experiment) (inherited from 
      \code{"\linkS4class{eSet}"}). See \code{\link{phenoData}} for more
      details. }
    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAPE}"}, containing details of experimental
      methods, including MIAPE relevant information. See
      \code{\link{experimentData}} for more details.}
    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{protocolData}} for
      more details. }
    \item{\code{assayData}:}{ not used here. }
    \item{\code{annotation}:}{ not used here. }
    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of R,
      \code{"linkS4class{Biobase}"}, \code{"\linkS4class{eSet}"} and
      \code{MSnExp} of the current instance. Intended for developer use
      and debugging (inherited from \code{"\linkS4class{eSet}"}). }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{eSet}"}, directly.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
  Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}

\section{Methods}{
  \describe{
    \item{\code{show(object)}}{ Displays object content as text. }
    \item{\code{header(object)}}{ Returns a \code{"data.frame"} with file
      information, retention time, precursor MZ, peaks count, total ion
      count (sum of all peaks in MS2 spectrum), ms level, precursor charge
      and collision energy. This method only works for MSMS experiments. }
    \item{\code{length(object)}}{ Returns the number of spectra in the
      experiment. }
    \item{\code{spectra(object)}}{ Return a \code{"list"} of
      \code{"\linkS4class{Spectrum}"} objects. }
    \item{\code{spectra(object) <- value}}{ Sets the list of spectra for
      the current \code{MSnExp} instance using \code{value}. }
    \item{\code{object[]}}{Subsets the \code{MSnExp} instance,
      effectively subsetting the \code{spectra} and \code{featureData}
      slots. }
    \item{\code{precursorMz(object)}}{ Returns an object of class
      \code{"numeric"} with the MZ values of the precursor ions of the
      MS2 spectra of the the \code{MSnExp} instance. }
    \item{\code{precursorCharge(object)}}{ Returns an object of class
      \code{"integer"} with the charge states of the precursor ions of
      the MS2 spectra of the the \code{MSnExp} instance. }
    \item{\code{acquisitionNum(object)}}{ Returns an object of class
      \code{"integer"} with the acquisition number of all the spectra of
      the \code{MSnExp} instance. }
    \item{\code{rtime(object)}}{ Returns an object of class
      \code{"numeric"} with the retention time of all the spectra in
      the \code{MSnExp} instance. }
    \item{\code{peaksCount(object)}}{ Returns an object of class
      \code{"integer"} with the peak count of all the spectra in
      the \code{MSnExp} instance. }
    \item{\code{msLevel(object)}}{ Returns an object of class
      \code{"integer"} with the MS level of all the spectra in
      the \code{MSnExp} instance. }
    \item{\code{collisionEnergy(object)}}{ Returns an object of class
      \code{"numeric"} with the collision energy used to fragment the
      precursor ions of all the MS2 spectra of the \code{MSnExp}
      instance. }
    \item{\code{intensity(object)}}{ Returns an object of class
      \code{"list"}, where each member is of class \code{"numeric"},
      containing the intensities of the spectra of the \code{MSnExp}
      instance. } 
    \item{\code{mz(object)}}{ Returns an object of class
      \code{"list"}, where each member is of class \code{"numeric"},
      containing the MZ value of the spectra peaks of the \code{MSnExp}
      instance. } 
    %% Methods to be documentated in own Rd file
    \item{\code{curveStats(object,reporters,...)}}{ Collects curve
      statistics for the reporter ions peaks of all the MS2 spectra of
      the \code{MSnExp} instance. See
      \code{\link{curveStats,MSnExp-method}} documentation for more
      details. } 
    \item{\code{quantify(object,reporters,method,...)}}{ Performs
      quantification for all the MS2 spectra of the \code{MSnExp}
      instance. See \code{\link{quantify,MSnExp-method}} documentation
      for more details. }
    \item{\code{trimMz(object,mzlim,...)}}{ Trims the MZ range of all
      the spectra of the \code{MSnExp} instance. See
      \code{\link{trimMz,MSnExp-method}} documentation for more
      details. }
    \item{\code{bg.correct(object,bg,...)}}{ Background corrects all the
      spectra of the \code{MSnExp} instance. See
      \code{\link{bg.correct,MSnExp-method}} documentation for more
      details. } 
    \item{\code{removePeaks(object,t,...)}}{ Removed peaks lower that a
      threshold \code{t}. See \code{\link{removePeaks,MSnExp-method}}
      documentation for more details. }
    \item{\code{clean}}{ Removes unused 0 intensity data points. See
      \code{\link{clean,MSnExp-method}} documentation for more details. } 
    \item{\code{plot(object,reporters,full)}}{ Plots all the spectra of
      the \code{MSnExp} instance. See \code{\link{plot,MSnExp-method}}
      documentation for more details. } 
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
  \code{\link{readMzXMLData}} for loading \code{mzXML} data to generate
  an instance of \code{MSnExp}. See \code{"\linkS4class{eSet}"},
  \code{"\linkS4class{MSnProcess}"}, \code{"\linkS4class{MIAPE}"},
  \code{\link{featureData}}, \code{\link{phenoData}},
  \code{\link{experimentData}}, \code{"\linkS4class{MIAME}"}, 
  \code{\link{protocolData}}, for more details about the
  individual slots.
}
\examples{
mzxmlfile <- dir(system.file("extdata",package="MSnbase"),pattern="mzXML",full.names=TRUE)
msnexp <- readMzXMLData(mzxmlfile)
msnexp
}

\keyword{classes}
