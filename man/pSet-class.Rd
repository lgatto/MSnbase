\name{pSet-class}
\Rdversion{1.1}
\docType{class}
\alias{pSet-class}
\alias{pSet}
\alias{class:pSet}

\alias{[,pSet-method}
\alias{[,pSet,ANY,ANY-method}
\alias{[,pSet,ANY,ANY,ANY-method}
\alias{[[,pSet-method}
\alias{[[,pSet,ANY,ANY-method}
\alias{$,pSet-method}
\alias{$<-,pSet-method}
\alias{abstract,pSet-method}
\alias{acquisitionNum,pSet-method}
\alias{scanIndex,pSet-method}
\alias{assayData,pSet-method}
\alias{collisionEnergy,pSet-method}
\alias{dim,pSet-method}
\alias{dim}
\alias{experimentData,pSet-method}
\alias{fData,pSet-method}
\alias{fData<-,pSet,data.frame-method}
\alias{featureData,pSet-method}
\alias{featureNames,pSet-method}
\alias{fileNames,pSet-method}
\alias{fileNames}
\alias{fromFile,pSet-method}
\alias{centroided,pSet-method}
\alias{centroided<-,pSet,ANY-method}
\alias{centroided<-,pSet,logical-method}
\alias{smoothed,pSet-method}
\alias{smoothed<-,pSet,ANY-method}
\alias{smoothed<-,pSet,logical-method}
\alias{fvarLabels,pSet-method}
\alias{fvarMetadata,pSet-method}
\alias{header,pSet,missing-method}
\alias{header,pSet,numeric-method}
\alias{header}
\alias{intensity,pSet-method}
\alias{length,pSet-method}
\alias{length}
\alias{msInfo,pSet-method}
\alias{msLevel,pSet-method}
\alias{mz,pSet-method}
\alias{notes,pSet-method}
\alias{pData,pSet-method}
\alias{pData<-,pSet,ANY-method}
\alias{peaksCount,pSet,missing-method}
\alias{peaksCount,pSet,numeric-method}
\alias{phenoData,pSet-method}
\alias{phenoData<-,pSet,ANY-method}
\alias{polarity,pSet-method}
\alias{precursorCharge,pSet-method}
\alias{precursorIntensity,pSet-method}
\alias{precursorMz,pSet-method}
\alias{precScanNum,pSet-method}
\alias{precAcquisitionNum,pSet-method}
\alias{processingData,pSet-method}
\alias{processingData}
\alias{protocolData,pSet-method}
\alias{pubMedIds,pSet-method}
\alias{rtime,pSet-method}
\alias{sampleNames,pSet-method}
\alias{sampleNames<-,pSet,character-method}
\alias{spectra,pSet-method}
\alias{spectra}
\alias{tic,pSet-method}
\alias{ionCount,pSet-method}
\alias{varLabels,pSet-method}
\alias{varMetadata,pSet-method}
\alias{exptitle,pSet-method}
\alias{expemail,pSet-method}
\alias{ionSource,pSet-method}
\alias{ionSourceDetails,pSet-method}
\alias{analyser,pSet-method}
\alias{analyzer,pSet-method}
\alias{analyserDetails,pSet-method}
\alias{analyzerDetails,pSet-method}
\alias{instrumentModel,pSet-method}
\alias{instrumentManufacturer,pSet-method}
\alias{instrumentCustomisations,pSet-method}
\alias{detectorType,pSet-method}
\alias{description,pSet-method}
\alias{isolationWindowLowerMz}
\alias{isolationWindowUpperMz}
\alias{isolationWindowLowerMz,pSet-method}
\alias{isolationWindowUpperMz,pSet-method}

\alias{spectrapply}
\alias{spectrapply,pSet-method}

\alias{coerce,AnnotatedDataFrame,list-method}
\alias{coerce,MIAxE,list-method}

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
      and \code{"\linkS4class{Spectrum2}"}). }
    \item{\code{phenoData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      experimenter-supplied variables describing sample (i.e the
      individual tags for an labelled MS experiment)
      See \code{\link{phenoData}} for more details.}
    \item{\code{featureData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing variables
      describing features (spectra in our case), e.g. identificaiton data,
      peptide sequence, identification score,... (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{featureData}} for
      more details. }
    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAPE}"}, containing details of experimental
      methods. See \code{\link{experimentData}} for more details. }
    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{protocolData}} for
      more details. }
    \item{\code{processingData}:}{Object of class
      \code{"\linkS4class{MSnProcess}"} that records all processing. }
    \item{\code{.cache}:}{Object of class \code{environment} used to
      cache data. Under development. }
    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of the class.
    }
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
    \item{$}{\code{signature(x = "pSet")}: directly access a specific
    sample annotation column from the \code{pData}.}
    \item{$<-}{\code{signature(x = "pSet")}: replace or add a
    sample annotation column in the \code{pData}.}
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
    \item{fromFile}{\code{signature(object = "pSet")}: Access raw data
      file indexes (to be found in the \code{processingData} slot) from
      which the individual object's spectra where read from. }
    \item{centroided}{\code{signature(object = "pSet")}: Indicates
      whether individual spectra are centroided ('TRUE') of uncentroided
      ('FALSE'). Use \code{centroided(object) <- value} to update a
      whole experiment, ensuring that \code{object} and \code{value}
      have the same length. }
    \item{smoothed}{\code{signature(object = "pSet")}: Indicates
      whether individual spectra are smoothed ('TRUE') of unsmoothed
      ('FALSE'). Use \code{smoothed(object) <- value} to update a
      whole experiment, ensuring that \code{object} and \code{value}
      have the same length. }
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
    \item{pData<-}{\code{signature(x = "pSet", value)}: Replace sample data
      information with \code{value}, value being a \code{data.frame}.}
    \item{phenoData}{\code{signature(x = "pSet")}: Access the
      \code{phenoData} slot. }
    \item{phenoData<-}{\code{signature(x = "pSet", value)}: Replace
    sample data information with \code{value}. \code{value} can be a
    \code{data.frame} or an \code{AnnotatedDataFrame}.}
    \item{processingData}{\code{signature(object = "pSet")}: Access the
      \code{processingData} slot. }
    \item{protocolData}{\code{signature(x = "pSet")}: Access the
      \code{protocolData} slot. }
    \item{pubMedIds}{\code{signature(x = "pSet")}: Access PMIDs in
      \code{experimentData}. }
    \item{sampleNames}{\code{signature(x = "pSet")}: Access sample names
      in \code{phenoData}. A replacement method is also available. }
    \item{spectra}{\code{signature(x = "pSet", ...)}: Access the
      \code{assayData} slot, returning the features as a \code{list}.
      Additional arguments are currently ignored. }
    \item{varMetadata}{\code{signature(x = "pSet")}: Access metadata
      describing variables reported in \code{pData}. }
    \item{varLabels}{\code{signature(x = "pSet")}: Access variable
      labels in \code{phenoData}. }
    \item{acquisitionNum}{\code{signature(object = "pSet")}: Accessor
      for spectra acquisition numbers. }
    \item{scanIndex}{\code{signature(object = "pSet")}: Accessor
      for spectra scan indices. }
    \item{collisionEnergy}{\code{signature(object = "pSet")}: Accessor
      for MS2 spectra collision energies. }
    \item{intensity}{\code{signature(object = "pSet", ...)}: Accessor
      for spectra instenities, returned as named list. Additional
      arguments are currently ignored. }
    \item{msInfo}{\code{signature(object = "pSet")}: Prints the MIAPE-MS
      meta-data stored in the \code{experimentData} slot. }
    \item{msLevel}{\code{signature(object = "pSet")}: Accessor for
      spectra MS levels. }
    \item{mz}{\code{signature(object = "pSet", ...)}: Accessor for spectra
      M/Z values, returned as a named list.  Additional arguments are
      currently ignored. }
    \item{peaksCount}{\code{signature(object = "pSet")}: Accessor for
      spectra preak counts. }
    \item{peaksCount}{\code{signature(object = "pSet", scans =
	"numeric")}: Accessor to \code{scans} spectra preak counts. }
    \item{polarity}{\code{signature(object = "pSet")}: Accessor for MS1
      spectra polarities. }
    \item{precursorCharge}{\code{signature(object = "pSet")}: Accessor
      for MS2 precursor charges. }
    \item{precursorIntensity}{\code{signature(object = "pSet")}: Accessor
      for MS2 precursor intensity. }
    \item{precursorMz}{\code{signature(object = "pSet")}: Accessor
      for MS2 precursor M/Z values. }
    \item{precAcquisitionNum}{\code{signature(object = "pSet")}: Accessor
      for MS2 precursor scan numbers. }
    \item{precScanNum}{see \code{precAcquisitionNum}. }
    \item{rtime}{\code{signature(object = "pSet", ...)}: Accessor for spectra
    retention times. Additional arguments are currently ignored. }
    \item{tic}{\code{signature(object = "pSet", ...)}: Accessor for spectra
      total ion counts. Additional arguments are currently ignored. }
    \item{ionCount}{\code{signature(object = "pSet")}: Accessor for spectra
      total ion current. }
    \item{header}{\code{signature(object = "pSet")}: Returns a data
      frame containing all available spectra parameters (MSn only). }
    \item{header}{\code{signature(object = "pSet", scans = "numeric")}:
      Returns a data frame containing \code{scans} spectra parameters
      (MSn only). }
    \item{spectrapply}{
      \code{spectrapply(signature(object = "pSet"), FUN = NULL,
	BPPARAM = bpparam(), ...)}: applies the function \code{FUN} to each
      spectrum passing additional parameters in \code{...} to that
      function and return its results. For \code{FUN = NULL} it returns
      the list of spectra (same as a call to \code{spectra}). Parameter
      \code{BPPARAM} allows to specify how and if parallel processing
      should be enabled.

      Returns a list with the result for each of spectrum.
    }

    \item{isolationWindowLowerMz}{
      \code{isolationWindowLowerMz(object = "pSet")}: return the
      lower m/z boundary for the isolation window. Note that this method
      is at present only available for \code{\linkS4class{OnDiskMSnExp}}
      objects.
    }

    \item{isolationWindowUpperMz}{
      \code{isolationWindowUpperMz(object = "pSet")}: return the
      upper m/z boundary for the isolation window. Note that this method
      is at present only available for \code{\linkS4class{OnDiskMSnExp}}
      objects.
    }

  }

  Additional accessors for the experimental metadata
  (\code{experimentData} slot) are defined. See
  \code{"\linkS4class{MIAPE}"} for details.
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
