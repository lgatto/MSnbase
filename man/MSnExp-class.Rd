\name{MSnExp-class}
\Rdversion{1.1}
\docType{class}
\concept{MSnExp}

\alias{MSnExp-class}
\alias{class:MSnExp}
\alias{MSnExp}


\alias{acquisitionNum,MSnExp-method}
\alias{[,MSnExp-method}
\alias{[[,MSnExp-method}
\alias{clean,MSnExp-method}
\alias{collisionEnergy,MSnExp-method}
\alias{fromFile,MSnExp-method}
\alias{header,MSnExp-method}
\alias{intensity,MSnExp-method}
\alias{length,MSnExp-method}
\alias{ms1scan,MSnExp-method}
\alias{msLevel,MSnExp-method}
\alias{mz,MSnExp-method}
\alias{peaksCount,MSnExp-method}
\alias{plot,MSnExp,missing-method}
\alias{polarity,MSnExp-method}
\alias{precursorCharge,MSnExp-method}
\alias{precursorMz,MSnExp-method}
\alias{proteomicsData<-,MSnExp,ANY-method}
\alias{proteomicsData<-,MSnExp,MIAPE-method}
\alias{proteomicsData,MSnExp-method}
\alias{quantify,MSnExp-method}
\alias{removePeaks,MSnExp-method}
\alias{rtime,MSnExp-method}
\alias{show,MSnExp-method}
\alias{spectra,MSnExp-method}
\alias{tic,MSnExp-method}
\alias{trimMz,MSnExp-method}

\title{The 'MSnExp' Class for MS Data And Meta-Data}


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


\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("MSnExp",...)}. However, it is preferred to use the
  \code{\link{readMzXMLData}} function that will read raw mass
  spectrometry data to generate a valid \code{"MSnExp"} instance.
}

\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"environment"}
      containing the MS spectra (see \code{"\linkS4class{Spectrum1}"}
      and \code{"\linkS4class{Spectrum2}"}).
      Slot is inherited from \code{"linkS4class{pSet}"}. }
    \item{\code{phenoData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing 
      experimenter-supplied variables describing sample (i.e the
      individual tags for an labelled MS experiment)
      See \code{\link{phenoData}} for more details.
      Slot is inherited from \code{"linkS4class{pSet}"}. }
    \item{\code{featureData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing variables
      describing features (spectra in our case), e.g. identificaiton data,
      peptide sequence, identification score,... (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{featureData}} for
      more details.
      Slot is inherited from \code{"linkS4class{pSet}"}. }
    \item{\code{proteomicsData}:}{Object of class
      \code{"\linkS4class{MIAPE}"} describing proteomics-specific
      experiment meta data. See also \code{\link{experimentData}}.
      Slot is inherited from \code{"linkS4class{pSet}"}. }	  
    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAPE}"}, containing details of experimental
      methods. See \code{\link{experimentData}} for more details.
      Slot is inherited from \code{"linkS4class{pSet}"}. }
    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{protocolData}} for
      more details.
      Slot is inherited from \code{"linkS4class{pSet}"}. }
    \item{\code{processindData}:}{Object of class
      \code{"\linkS4class{MSnProcess}"} that records all processing. 
      Slot is inherited from \code{"linkS4class{pSet}"}. }
    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of R,
      \code{"linkS4class{Biobase}"}, \code{"\linkS4class{pSet}"} and
      \code{MSnExp} of the current instance.
      Slot is inherited from \code{"linkS4class{pSet}"}. 
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
  \describe{
    \item{show}{\code{signature(object = "MSnExp")}: Displays object
      content as text. }
    \item{acquisitionNum}{\code{signature(object = "MSnExp")}: Returns
      an object of class \code{"integer"} with the acquisition number of
      all the spectra of the \code{MSnExp} instance. }
    \item{[}{\code{signature(x = "MSnExp")}: Subsets the \code{MSnExp}
      instance, effectively subsetting the \code{spectra} and
      \code{featureData} slots to return a new \code{MSnExp} instant. }      
    \item{[[}{\code{signature(x = "MSnExp")}: Subsets the individual
      spectrum in the \code{MSnExp} \code{AssayData} slot by name
      ('character') or index ('numeric'), returning an instance of class
      \code{"\linkS4class{Spectrum1}"} or
      \code{"\linkS4class{Spectrum2}"}, depending on the MS level of the
      object. }
    \item{header}{\code{signature(object = "MSnExp")}: Returns a
      \code{"data.frame"} with file information, retention time,
      precursor MZ, peaks count, total ion count (sum of all peaks in
      MS2 spectrum), ms level, precursor charge and collision
      energy. This method only works for MSMS experiments. }
    \item{length}{\code{signature(object = "MSnExp")}: Returns the
      number of spectra in the experiment. }
    \item{spectra}{\code{signature(object = "MSnExp")} :Return a
      \code{"list"} of \code{"\linkS4class{Spectrum}"} objects. }
    \item{clean}{\code{signature(object = "MSnExp")}: Removes unused 0
      intensity data points. See \code{\link{clean}} documentation
      for more details and examples. }  
    \item{removePeaks}{\code{signature(object = "MSnExp")}: Removed
      peaks lower that a threshold \code{t}. See
      \code{\link{removePeaks}} documentation for more details and
      examples. } 
    \item{collisionEnergy}{\code{signature(object = "MSnExp")}: Returns
      an object of class \code{"numeric"} with the collision energy used
      to fragment the  precursor ions of all the MS2 spectra of the
      \code{MSnExp} instance. }
    \item{fromFile}{\code{signature(object = "MSnExp")}: Returns an
      object of class \code{"numeric"} indexing the original file names
      (to be found in the 'code{processingData}' slot) from which the
      individual object's spectra where read from. }
    \item{precursorMz}{\code{signature(object = "MSnExp")}: Returns an
      object of class \code{"numeric"} with the MZ values of the
      precursor ions of the MS2 spectra of the the \code{MSnExp}
      instance. }
    \item{precursorCharge}{\code{signature(object = "MSnExp")}: Returns
      an object of class \code{"integer"} with the charge states of the
      precursor ions of the MS2 spectra of the the \code{MSnExp}
      instance. }
    \item{rtime}{\code{signature(object = "MSnExp")}: Returns an object
      of class \code{"numeric"} with the retention time of all the
      spectra in the \code{MSnExp} instance. }
    \item{peaksCount}{\code{signature(object = "MSnExp")}: Returns an
      object of class \code{"integer"} with the peak count of all the
      spectra in the \code{MSnExp} instance. }
    \item{msLevel}{\code{signature(object = "MSnExp")}: Returns an
      object of class \code{"integer"} with the MS level of all the
      spectra in the \code{MSnExp} instance. }
    \item{polarity}{\code{signature(object = "MSnExp")}: Returns an
      object of class \code{"integer"} with the polarity of all the MS1
      spectra in the \code{MSnExp} instance. }
    \item{tic}{\code{signature(object = "MSnExp")}: Returns an
      object of class \code{"numeric"} with the total ion count of all the
      spectra in the \code{MSnExp} instance. }
    \item{intensity}{\code{signature(object = "MSnExp")}: Returns an
      object of class \code{"list"}, where each member is of class
      \code{"numeric"}, containing the intensities of the spectra of the
      \code{MSnExp} instance. } 
    \item{mz}{\code{signature(object = "MSnExp")}: Returns an object of
      class \code{"list"}, where each member is of class
      \code{"numeric"}, containing the MZ value of the spectra peaks of
      the \code{MSnExp} instance. } 
    \item{trimMz}{\code{signature(object = "MSnExp")}: Trims the MZ
      range of all the spectra of the \code{MSnExp} instance. See
      \code{\link{trimMz}} documentation for more details and
      examples. } 


    



    
    \item{quantify}{\code{signature(object = "MSnExp")}: Performs
      quantification for all the MS2 spectra of the \code{MSnExp}
      instance. See \code{\link{quantify}} documentation for more
      details. } 

    \item{plot}{\code{signature(x = "MSnExp", y = "missing")}: Plots all
      the spectra of the \code{MSnExp} instance. See
      \code{\link{plot.MSnExp}} documentation for more details. } 

    
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
  an instance of \code{MSnExp}.
}

\examples{
mzxmlfile <- dir(system.file("extdata",package="MSnbase"),
                 pattern="mzXML",full.names=TRUE)
msnexp <- readMzXMLData(mzxmlfile)
msnexp
}

\keyword{classes}
