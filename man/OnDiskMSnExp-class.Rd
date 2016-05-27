\name{OnDiskMSnExp-class}
\Rdversion{1.1}
\docType{class}
\concept{OnDiskMSnExp}

\alias{OnDiskMSnExp-class}
\alias{class:OnDiskMSnExp}
\alias{OnDiskMSnExp}
\alias{show,OnDiskMSnExp-method}

%% Getter/Setter
\alias{[,OnDiskMSnExp,ANY,ANY,ANY-method}
\alias{[,OnDiskMSnExp,logicalOrNumeric,logicalOrNumeric,missing-method}
\alias{[,OnDiskMSnExp,logicalOrNumeric,missing,missing-method}
\alias{[,OnDiskMSnExp,missing,logicalOrNumeric,missing-method}
\alias{[[,OnDiskMSnExp-method}
\alias{[[,OnDiskMSnExp,ANY,ANY-method}
\alias{acquisitionNum,OnDiskMSnExp-method}
\alias{assayData,OnDiskMSnExp-method}
\alias{centroided,OnDiskMSnExp-method}
\alias{centroided<-,OnDiskMSnExp,logical-method}
\alias{fromFile,OnDiskMSnExp-method}
\alias{header,OnDiskMSnExp,missing-method}
\alias{header,OnDiskMSnExp,numeric-method}
\alias{ionCount,OnDiskMSnExp-method}
\alias{intensity,OnDiskMSnExp-method}
\alias{length,OnDiskMSnExp-method}
\alias{msLevel,OnDiskMSnExp-method}
\alias{mz,OnDiskMSnExp-method}
\alias{peaksCount,OnDiskMSnExp,missing-method}
\alias{peaksCount,OnDiskMSnExp,numeric-method}
\alias{polarity,OnDiskMSnExp-method}
\alias{rtime,OnDiskMSnExp-method}
\alias{scanIndex,OnDiskMSnExp-method}
\alias{spectra,OnDiskMSnExp-method}
\alias{tic,OnDiskMSnExp-method}

%% Data manipulations.
\alias{clean,OnDiskMSnExp-method}
\alias{removePeaks,OnDiskMSnExp-method}

\title{The 'OnDiskMSnExp' Class for MS Data And Meta-Data}

\description{
  Like the \code{\link{MSnExp}} class, the \code{OnDiskMSnExp} class
  encapsulates data and meta-data for mass spectrometry
  experiments, but does, in contrast to the former, not keep the
  spectrum data in memory, but fetches the M/Z and intensity values on
  demand from the raw files. This results in some instances to a
  reduced performance, has however the advantage of a much smaller
  memory footprint.

  At present only on-disk access for \code{\link{Spectrum1}} objects is
  supported.
}

\details{
  The \code{OnDiskMSnExp} object stores many spectrum related
  information into the \code{featureData}, thus, some calls, like
  \code{rtime} to retrieve the retention time of the individual scans
  does not require the raw data to be read. Only M/Z and intensity
  values are loaded on-the-fly from the original files. Extraction of
  values for individual scans is, for mzML files, very fast. Extraction
  of the full data (all spectra) are performed in a per-file parallel
  processing strategy.

  Data manipulations related to spectras' M/Z or intensity values
  (e.g. \code{\link{removePeaks}} or \code{\link{clean}}) are (for
  \code{OnDiskMSnExp} objects) not applied immediately, but are stored
  for later execution into the \code{spectraProcessingQueue}. The
  manipulations are performed \emph{on-the-fly} upon data retrieval.
  Other manipulations, like removal of individual spectra are applied
  directly, since the corresponding data is available in the object's
  \code{featureData} slot.

}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("OnDiskMSnExp",...)}. However, it is preferred to use the
  \code{\link{readMSData}} function with argument \code{backend="disk"}
  that will read raw mass spectrometry data to generate a valid
  \code{"OnDiskMSnExp"} instance.
}

\section{Slots}{
  \describe{
    \item{\code{backend}:}{
      Character string specifying the used backend.
    }
    \item{\code{spectraProcessingQueue}:}{
      \code{list} of \code{\linkS4class{ProcessingStep}} objects
      defining the functions to be applied \emph{on-the-fly} to the
      spectra data (M/Z and intensity duplets).
    }
    \item{\code{assayData}:}{Object of class \code{"environment"} that
      is however empty, as no spectrum data is stored.
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{phenoData}:}{Object of class
      \code{"\linkS4class{NAnnotatedDataFrame}"} containing
      experimenter-supplied variables describing sample (i.e the
      individual tags for an labelled MS experiment)
      See \code{\link{phenoData}} for more details.
      Slot is inherited from \code{"\linkS4class{pSet}"}. }
    \item{\code{featureData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing variables
      describing features (spectra in our case). See
      \code{\link{featureData}} for more details.
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
Class \code{"\linkS4class{MSnExp}"}, directly.
Class \code{"\linkS4class{pSet}"}, by class "MSnExp", distance 3.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "pSet", distance 4.
Class \code{"\linkS4class{Versioned}"}, by class "pSet", distance 5.
}

\section{Getter/setter methods}{
  (in alphabetical order)
  See also methods for \code{\linkS4class{MSnExp}} or
  \code{\linkS4class{pSet}} objects.

  \describe{

    \item{[}{
      \code{object[i, j]}: subset the \code{OnDiskMSnExp} by spectra
      (\code{i}) or samples (\code{j}).
      \code{i} can be a \code{numeric} or \code{logical} vector
      specifying to which spectra the data set should be reduced (with
      \code{i} being the index of the spectrum in the object's
      \code{featureData}).
      \code{j} can be a \code{numeric} or \code{logical} vector defining
      the samples/files to which the data set should be reduced.

      The method returns a \code{OnDiskMSnExp} object with the data
      sub-set.
    }

    \item{[[}{
      \code{object[[i]]}: extract s single spectrum from the
      \code{OnDiskMSnExp} object \code{object}. Argument \code{i} can be
      either numeric or character specifying the index or the name of
      the spectrum in the object (i.e. in the \code{featureData}). The
      relevant information will be extracted from the corresponding raw
      data file.

      The method returns a \code{Spectrum1} object.
    }
    \item{acquisitionNum}{
      \code{acquisitionNum(signature(object="OnDiskMSnExp"))}:
      get the acquisition number of each spectrum in each individual
      file. The relevant information is extracted from the object's
      \code{featureData} slot.

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

    \item{assayData}{
      \code{assayData(signature(object = "OnDiskMSnExp"))}:
      Extract the full data, i.e. read all spectra from the original files,
      apply all processing steps from the \code{spectraProcessingQueue}
      slot and return the data. Due to the required processing time
      accessing the full data should be avoided wherever possible.

      Returns an \code{environment}.
    }

    \item{centroided,centroided<-}{
      \code{centroided(signature(object="OnDiskMSnExp"))}:
      whether individual spectra are centroided or uncentroided. The
      relevant information is extracted from the object's
      \code{featureData} slot. Returns a logical vector with names
      corresponding to the spectrum names.
      Use \code{centroided(object) <- value} to update the information,
      with value being a logical vector of length equal to the number of
      spectra in the experiment.
    }

    \item{fromFile}{
      \code{fromFile(signature(object = "OnDiskMSnExp"))}: get the
      index of the file (in \code{fileNames(object)}) from which the
      spectra were read. The relevant information is extracted from the
      object's \code{featureData} slot.

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

    \item{header}{
      \code{header(signature(object = "OnDiskMSnExp"))}: extract
      header information for all spectra. The relevant information is
      extracted from the object's \code{featureData} slot.

      Returns a \code{data.frame} with columns \code{"file"}, \code{"spIdx"}
      (spectrum index), \code{"centroided"}, \code{"basePeakMZ"}, \code{"polarity"},
      \code{"retention.time"}, \code{"peaks.count"}, \code{"tic"},
      \code{"ms.level"} and \code{"acquisition.number"}.
    }

    \item{header}{
      \code{header(signature(object="OnDiskMSnExp", scans="numeric"))}:
      extract header information for selected scans/spectra. See
      \code{header} method above for more details.
    }

    \item{intensity}{
      \code{intensity(signature(object="OnDiskMSnExp"))}:
      return the intensities from each spectrum in the data
      set. Intensities are first read from the raw files followed by an
      optional processing (depending on the processing steps defined in
      the \code{spectraProcessingQueue}). To reduce the amount of
      required  memory, this is performed on a per-file basis.
      The \code{BPPARAM} argument allows to specify how
      and if parallel processing should be used. Information from
      individual files will be processed in parallel (one process per
      original file).

      The method returns a \code{list} of numeric intensity values. Each
      list element represents the intensities from one spectrum.
    }

    \item{ionCount}{
      \code{ionCount(signature(object="OnDiskMSnExp",
	BPPARAM=bpparam()))}:
      extract the ion count (i.e. sum of intensity values) for each
      spectrum in the data set. The relevant data has to be extracted
      from the raw files (with eventually applying processing steps).
      The \code{BPPARAM} argument can be used to define how
      and if parallel processing should be used. Information from
      individual files will be processed in parallel (one process per
      original file).

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

    \item{length}{
      \code{length(signature(object="OnDiskMSnExp"))}:
      Returns the number of spectra of the current experiment.
    }

    \item{msLevel}{
      \code{msLevel(signature(object = "OnDiskMSnExp"))}: extract the
      MS level from the spectra. The relevant information is extracted
      from the object's \code{featureData} slot.

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

    \item{mz}{
      \code{mz(signature(object="OnDiskMSnExp"))}:
      return the M/Z values from each spectrum in the data
      set. M/Z values are first read from the raw files followed by an
      optional processing (depending on the processing steps defined in
      the \code{spectraProcessingQueue}). To reduce the amount of
      required  memory, this is performed on a per-file basis.
      The \code{BPPARAM} argument allows to specify how
      and if parallel processing should be used. Information from
      individual files will be processed in parallel (one process per
      original file).

      The method returns a \code{list} of numeric M/Z values. Each
      list element represents the values from one spectrum.
    }

    \item{peaksCount}{
      \code{peaksCount(signature(object="OnDiskMSnExp",
	scans="numeric"), BPPARAM=bpparam())}:
      extrac the peaks count from each spectrum in the object. Depending
      on the eventually present \code{ProcessingStep} objects in the
      \code{spectraProcessingQueue} raw data will be loaded to calculate
      the peaks count. If no steps are present, the data is extracted
      from the \code{featureData}. Optional argument \code{scans} allows
      to specify the index of specific spectra from which the count
      should be returned. The \code{BPPARAM} argument can be used to define how
      and if parallel processing should be used. Information from
      individual files will be processed in parallel (one process per
      original file).

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

    \item{polarity}{
      \code{polarity(signature(object="OnDiskMSnExp"))}:
      returns a numeric vector with the polarity of the individual
      spectra in the data set. The relevant information is extracted
      from the \code{featureData}.
    }

    \item{rtime}{
      \code{rtime(signature(object="OnDiskMSnExp"))}:
      extrac the retention time of the individual spectra in the data
      set (from the \code{featureData}).

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

    \item{scanIndex}{
      \code{scanIndex(signature(object="OnDiskMSnExp"))}:
      get the spectra scan indices within the respective file. The
      relevant information is extracted from the object's \code{featureData} slot.
      Returns a numeric vector of indices with names corresponding to the
      spectrum names.
    }

    \item{spectra}{
      \code{spectra(signature(object="OnDiskMSnExp"), scans=NULL,
	rtlim=NULL, BPPARAM=bpparam())}:
      extract spectrum data from the individual files. This causes the
      spectrum data to be read from the original raw files. After that
      all processing steps defined in the \code{spectraProcessingQueue}
      are applied to it. The results are then returned as a \code{list}
      of \code{\linkS4class{Spectrum1}} objects.

      Argument \code{scans}: numeric, character or logical vector
      specifying the spectra that should be extracted (either index in
      the data set (not scan index within the file!), name of the
      spectrum or logical vector with length equal to the number of
      spectra in the data set).
      Argument \code{rtlim}: numeric vector of length 2 defining a
      retention time window/range from which spectra should be
      extracted.
      The \code{BPPARAM} argument can be used to define how and if
      parallel processing should be used. Information from individual
      files will be processed in parallel (one process per file).
      Note: extraction of selected spectra results in a considerable
      processing speed and should thus be preferred over whole data
      extraction.

      Returns a \code{list} of \code{\linkS4class{Spectrum1}} objects
      with names corresponding to the spectrum names.
    }

    \item{tic}{
      \code{tic(signature(object="OnDiskMSnExp"))}:
      get the total ion current of each spectrum in the data
      set. This information is extracted from the object's
      \code{featureData} and represents the tic provided in the original
      raw data files.

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

  }

}

\section{Data manipulation methods}{
  (in alphabetical order)
  See also methods for \code{\linkS4class{MSnExp}} or
  \code{\linkS4class{pSet}} objects.

  \describe{

    \item{clean}{
      \code{clean(signature(object="OnDiskMSnExp"), all=TRUE, verbose=TRUE)}:
      add an \emph{clean} processing step to the lazy processing queue
      of the \code{OnDiskMSnExp} object. The \code{clean} command will
      only be executed when spectra information (including M/Z and
      intensity values) is requested from the \code{OnDiskMSnExp}
      object. Optional arguments to the methods are \code{all=TRUE} and
      \code{verbose=TRUE}.

      The method returns an \code{OnDiskMSnExp} object.

      For more details see documentation of the \code{\link{clean}}
      method.
    }

    \item{removePeaks}{
      \code{removePeaks(signature(object="OnDiskMSnExp"), t="min", verbose=TRUE)}:
      add an \code{removePeaks} processing step to the lazy processing
      queue of the returned \code{OnDiskMSnExp} object.
    }

  }

}

\author{
  Johannes Rainer <johannes.rainer@eurac.edu>
}

\seealso{
  \code{\linkS4class{pSet}},
  \code{\linkS4class{MSnExp}},
  \code{\link{readMSData}}
}

\examples{
## Get some example mzML files
library(msdata)
mzfiles <- c(system.file("microtofq/MM14.mzML", package="msdata"),
             system.file("microtofq/MM8.mzML", package="msdata"))
## Read the data as an OnDiskMSnExp; we ignore the warning.
odmse <- readMSData(mzfiles, msLevel=1)

## Get the length of data, i.e. the total number of spectra.
length(odmse)

## Get the MS level
head(msLevel(odmse))

## Get the featureData, use fData to return as a data.frame
head(fData(odmse))

## Get to know from which file the spectra are
head(fromFile(odmse))

## And the file names:
fileNames(odmse)

## Scan index and acquisitionNum
head(scanIndex(odmse))
head(acquisitionNum(odmse))

## Extract the spectra; the data is retrieved from the raw files.
head(spectra(odmse))

## Extracting individual spectra or a subset is much faster.
spectra(odmse, scans=c(1, 50))

## Spectra can also be extracted given a retention time window
spectra(odmse, rtlim=c(2, 20))

## Alternatively, we could also subset the whole object by spectra and/or samples:
subs <- odmse[rtime(odmse) >= 2 & rtime(odmse) <= 20, 1]
fileNames(subs)
rtime(subs)

## Extract intensities and M/Z values per spectrum; the methods return a list,
## each element representing the values for one spectrum.
ints <- intensity(odmse)
mzs <- mz(odmse)


## Calling removePeaks, i.e. setting intensity values below a certain threshold to 0.
## Unlike the name suggests, this is not actually removing peaks. Such peaks with a 0
## intensity are then removed by the "clean" step.
## Also, the manipulations are not applied directly, but put into the "lazy"
## processing queue.
odmse <- removePeaks(odmse, t=10000)
odmse <- clean(odmse)

## The processing steps are only applied when actual raw data is extracted.
spectra(odmse, scans=c(1, 50))

## Get the polarity of the spectra.
head(polarity(odmse))

## Get the retention time of all spectra
head(rtime(odmse))

## Get the intensities after removePeaks and clean
intsAfter <- intensity(odmse)

head(lengths(ints))
head(lengths(intsAfter))

## The same for the M/Z values
mzsAfter <- intensity(odmse)
head(lengths(mzs))
head(lengths(mzsAfter))

}

\keyword{classes}