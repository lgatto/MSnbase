\name{pSet-class}
\Rdversion{1.1}
\docType{class}
\alias{pSet-class}
\alias{pSet}
\alias{class:pSet}

\alias{[,pSet-method}
\alias{[[,pSet-method}
\alias{abstract,pSet-method}
\alias{acquisitionNum,pSet-method}
\alias{assayData,pSet-method}
\alias{collisionEnergy,pSet-method}
\alias{description,pSet-method}
\alias{dim,pSet-method}
\alias{dim}
\alias{experimentData,pSet-method}
\alias{fData,pSet-method}
\alias{featureData,pSet-method}
\alias{featureNames,pSet-method}
\alias{fileNames,pSet-method}
\alias{fileNames}
\alias{fromFile,pSet-method}
\alias{fvarLabels,pSet-method}
\alias{fvarMetadata,pSet-method}
\alias{header,pSet-method}
\alias{header}
\alias{intensity,pSet-method}
\alias{length,pSet-method}
\alias{length}
\alias{msLevel,pSet-method}
\alias{mz,pSet-method}
\alias{notes,pSet-method}
\alias{pData,pSet-method}
\alias{peaksCount,pSet-method}
\alias{phenoData,pSet-method}
\alias{polarity,pSet-method}
\alias{precursorCharge,pSet-method}
\alias{precursorMz,pSet-method}
\alias{processingData,pSet-method}
\alias{processingData}
\alias{proteomicsData,pSet-method}
\alias{proteomicsData,MSnExp-method}
\alias{proteomicsData}
\alias{proteomicsData<-,pSet-method}
\alias{proteomicsData<-,MSnExp-method}
\alias{proteomicsData<-,MSnExp,ANY-method}
\alias{proteomicsData<-,MSnExp,MIAPE-method}
\alias{proteomicsData<-}
\alias{protocolData,pSet-method}
\alias{pubMedIds,pSet-method}
\alias{rtime,pSet-method}
\alias{sampleNames,pSet-method}
\alias{spectra,pSet-method}
\alias{spectra,MSnExp-method}
\alias{spectra}
\alias{tic,pSet-method}
\alias{varLabels,pSet-method}
\alias{varMetadata,pSet-method}

\title{
  Class to Contain Raw Mass-Spectrometry Assays and Experimental
  Metadata
}

\description{
  Container for high-throughput mass-spectrometry assays and
  experimental metadata. This class is based on Biobase's
  \code{"\linkS4class{eSet}"} virtual class, with the notable exception
  that 'assayData' slot is an environment contain objects of class
  \code{"\linkS4class{Spectrum}"}. 
}

\section{Objects from the Class}{
  A virtual Class: No objects may be created from it.  
  See \code{"\linkS4class{MSnExp}"} for instantiatable sub-classes.
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
    \item{\code{processingData}:}{Object of class
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
  Class \code{"\linkS4class{VersionedBiobase}"}, directly.
  Class \code{"\linkS4class{Versioned}"}, by class "VersionedBiobase", distance 2.
}


\section{Methods}{
  Methods defined in derived classes may override the methods described
  here. 
  
  \describe{
    \item{[}{\code{signature(x = "pSet")}: Subset current object and
      return object of same class. }
    \item{[[}{\code{signature(x = "pSet")}: Direct access to individual
      spectra. }
    \item{abstract}{Access abstract in \code{experimentData}. }
    \item{assayData}{\code{signature(object = "pSet")}: Access the
      \code{assayData} slot. Returns an \code{environment}. }    
    \item{desciption}{\code{signature(x = "pSet")}: Synonymous with
      experimentData. } 
    \item{dim}{\code{signature(x = "pSet")}: Returns the dimensions of
      the \code{phenoData} slot. }
    \item{experimentData}{\code{signature(x = "pSet")}: Access details 
      of experimental methods. }
    \item{featureData}{\code{signature(x = "pSet")}: Access the
      \code{featureData} slot. }
    \item{fData}{\code{signature(x = "pSet")}: Access feature data
      information. }
    \item{featureNames}{\code{signature(x = "pSet")}: Coordinate access
      of feature names (e.g spectra, peptides or proteins) in
      \code{assayData} slot. }
    \item{fileNames}{\code{signature(object = "pSet")}: Access file
      names in the \code{processingData} slot. }
    \item{fromFile}{\code{signature(object = "MSnExp")}: Access raw data
      file indexes (to be found in the 'code{processingData}' slot) from
      which the individual object's spectra where read from. }
    \item{fvarMetadata}{\code{signature(x = "pSet")}: Access metadata
      describing features reported in \code{fData}. }
    \item{fvarLabels}{\code{signature(x = "pSet")}: Access variable
      labels in \code{featureData}. }    
    \item{length}{\code{signature(x = "pSet")}: Returns the number of
      features in the \code{assayData} slot. }
    \item{notes}{\code{signature(x = "pSet")}:  Retrieve and
      unstructured notes associated with \code{pSet} in the
    \code{experimentData} slot. }
    \item{pData}{\code{signature(x = "pSet")}: Access sample data
      information. } 
    \item{phenoData}{\code{signature(x = "pSet")}: Access the
      \code{phenoData} slot. }
    \item{processingData}{\code{signature(object = "pSet")}: Access the
      \code{processingData} slot. }    
    \item{proteomicsData}{\code{signature(x = "pSet")}: Access the
      \code{proteomicsData} slot. }
    \item{proteomicsData<-}{\code{signature(x = "pSet", value="MIAPE")}:
      Set the \code{proteomicsData} slot. }
    \item{protocolData}{\code{signature(x = "pSet")}: Access the
      \code{protocolData} slot. }
    \item{pubMedIds}{\code{signature(x = "pSet")}: Access PMIDs in
      \code{experimentData}. } 
    \item{sampleNames}{\code{signature(x = "pSet")}: Access sample names
      in \code{phenoData}. }
    \item{spectra}{\code{signature(x = "pSet")}: Access the
      \code{assayData} slot, returning the features as a \code{list}. } 
    \item{varMetadata}{\code{signature(x = "pSet")}: Access metadata
      describing variables reported in \code{pData}. }
    \item{varLabels}{\code{signature(x = "pSet")}: Access variable
      labels in \code{phenoData}. }
    \item{acquisitionNum}{\code{signature(object = "pSet")}: Accessor
      for spectra acquisition numbers. }
    \item{collisionEnergy}{\code{signature(object = "pSet")}: Accessor
      for MS2 spectra collision energies. }
    \item{intensity}{\code{signature(object = "pSet")}: Accessor for
    spectra instenities, returned as named list. }
    \item{msLevel}{\code{signature(object = "pSet")}: Accessor for
      spectra MS levels. }
    \item{mz}{\code{signature(object = "pSet")}: Accessor for spectra
      M/Z values, returned as a named list. }
    \item{peaksCount}{\code{signature(object = "pSet")}: Accessor for
      spectra preak counts. }
    \item{polarity}{\code{signature(object = "pSet")}: Accessor for MS1
      spectra polarities. }
    \item{precursorCharge}{\code{signature(object = "pSet")}: Accessor
      for MS2 precursor charges. }
    \item{precursorMz}{\code{signature(object = "pSet")}: Accessor
      for MS2 precursor M/Z values. }
    \item{rtime}{\code{signature(object = "pSet")}: Accessor for spectra
    retention times. }
    \item{tic}{\code{signature(object = "pSet")}: Accessor for spectra
      total ion counts. }
    \item{header}{\code{signature(object = "pSet")}: Returns a data
      frame containing all available spectra parameters (MSn only). }
  }			      
}


\references{
  The \code{"\linkS4class{eSet}"} class, on which \code{pSet} is based.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{"\linkS4class{MSnExp}"} for an instantiatable application of
  \code{pSet}. 
}

\examples{
showClass("pSet")
}

\keyword{classes}
