\name{OnDiskMSnExp-class}
\Rdversion{1.1}
\docType{class}
\concept{OnDiskMSnExp}

\alias{OnDiskMSnExp-class}
\alias{class:OnDiskMSnExp}
\alias{OnDiskMSnExp}
\alias{show,OnDiskMSnExp-method}

%% Getter/Setter
\alias{acquisitionNum,OnDiskMSnExp-method}
\alias{assayData,OnDiskMSnExp-method}
\alias{centroided,OnDiskMSnExp-method}
\alias{centroided<-,OnDiskMSnExp,logical-method}
\alias{fromFile,OnDiskMSnExp-method}
\alias{header,OnDiskMSnExp-method}
\alias{header,OnDiskMSnExp,numeric-method}
\alias{ionCount,OnDiskMSnExp-method}
\alias{length,OnDiskMSnExp-method}
\alias{msLevel,OnDiskMSnExp-method}
\alias{peaksCount,OnDiskMSnExp-method}
\alias{peaksCount,OnDiskMSnExp,numeric-method}
\alias{polarity,OnDiskMSnExp-method}
\alias{rtime,OnDiskMSnExp-method}
\alias{scanIndex,OnDiskMSnExp-method}
\alias{spectra,OnDiskMSnExp-method}
\alias{tic,OnDiskMSnExp-method}

%% Data manipulations.
\alias{clean,OnDiskMSnExp-method}
\alias{removePeaks,OnDiskMSnExp-method}

%% \alias{bin,MSnExp-method}
%% \alias{clean,MSnExp-method}
%% \alias{compareSpectra,MSnExp,missing-method}
%% \alias{extractPrecSpectra,MSnExp-method}
%% \alias{extractPrecSpectra,MSnExp,numeric-method}
%% %% \alias{extractSpectra,MSnExp-method}
%% %% \alias{extractSpectra,MSnExp,logical-method}
%% \alias{pickPeaks,MSnExp-method}
%% \alias{estimateNoise,MSnExp-method}
%% \alias{plot,MSnExp}
%% \alias{plot,MSnExp,missing-method}
%% \alias{quantify,MSnExp-method}
%% \alias{quantify,MSnExp,character-method}
%% \alias{removePeaks,MSnExp-method}
%% \alias{removeReporters,MSnExp-method}
%% \alias{smooth,MSnExp-method}
%% \alias{addIdentificationData,MSnExp,character-method}
%% \alias{addIdentificationData,MSnExp,mzIDClasses-method}
%% \alias{addIdentificationData,MSnExp,mzID-method}
%% \alias{addIdentificationData,MSnExp,mzIDCollection-method}
%% \alias{addIdentificationData,MSnExp,data.frame-method}
%% \alias{removeNoId,MSnExp-method}
%% \alias{removeMultipleAssignment,MSnExp-method}
%% \alias{idSummary,MSnExp-method}
%% \alias{show,MSnExp-method}
%% \alias{spectra,MSnExp-method}
%% \alias{trimMz,MSnExp-method}
%% \alias{trimMz,MSnExp,numeric-method}

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
	BPPARAM=bpparam())}:
      extract spectrum data from the individual files. Spectrum data is
      retrieved from the raw files, any processing steps in the
      \code{spectraProcessingQueue} are applied to the data which is
      then returned as \code{\linkS4class{Spectrum1}} objects. Argument
      \code{scans} can be used to extract only specific spectra (which
      considerably increases processing speed) by specifying their index
      in the object's \code{featureData}. The \code{BPPARAM} argument
      can be used to define how and if parallel processing should be
      used. Information from individual files will be processed in
      parallel (one process per file).

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

%% \section{Methods}{
%%   See the \code{"\linkS4class{pSet}"} class for documentation on
%%   accessors inherited from \code{pSet}, subsetting and general attribute
%%   accession.

%%   \describe{
%%     \item{bin}{\code{signature(object = "MSnExp")}: Bins spectra.
%%       See \code{\link{bin}} documentation for more details and examples. }
%%     \item{clean}{\code{signature(object = "MSnExp")}: Removes unused 0
%%       intensity data points. See \code{\link{clean}} documentation
%%       for more details and examples. }
%%     \item{compareSpectra}{\code{signature(object1 = "Spectrum",
%%       object2 = "missing")}: Compares spectra.  See
%%       \code{\link{compareSpectra}} documentation for more details and
%%       examples. }
%%     \item{extractPrecSpectra}{\code{signature(object = "MSnExp", prec =
%% 	"numeric")}: extracts spectra with precursor MZ value equal to
%%       \code{prec} and returns an object of class 'MSnExp'.
%%       See \code{\link{extractPrecSpectra}} documentation for more
%%       details and examples. }
%%     % \item{extractSpectra}{\code{signature(object = "MSnExp", selected =
%%     % 	"logical")}: extracts spectra corresponding to 'TRUE' in
%%     %   \code{selected}.
%%     % See \code{\link{extractSpectra}} documentation for more details
%%     % and examples. }
%%     \item{pickPeaks}{\code{signature(object = "MSnExp")}: Performs the peak
%%       picking to generate centroided spectra.
%%       See \code{\link{pickPeaks}} documentation for more
%%       details and examples. }
%%     \item{estimateNoise}{\code{signature(object = "MSnExp")}: Estimates
%%       the noise in all profile spectra of \code{object}.  See
%%       \code{\link{estimateNoise}} documentation for more details and
%%       examples. }
%%     \item{plot}{\code{signature(x = "MSnExp", y = "missing")}: Plots all
%%       the spectra of the \code{MSnExp} instance. See
%%       \code{\link{plot.MSnExp}} documentation for more details. }
%%     \item{plot2d}{\code{signature(object = "MSnExp", ...)}:
%%       Plots retention time against precursor MZ for \code{MSnExp}
%%       instances. See \code{\link{plot2d}} documentation for more
%%       details. }
%%     \item{plotDensity}{\code{signature(object = "MSnExp", ...)}:
%%       Plots the density of parameters of interest.
%%       instances. See \code{\link{plotDensity}} documentation for more
%%       details. }
%%     \item{plotMzDelta}{\code{signature(object = "MSnExp", ...)}:
%%       Plots a histogram of the m/z difference betwee all of the highest
%%       peaks of all MS2 spectra of an experiment.
%%       See \code{\link{plotMzDelta}} documentation for more details. }
%%     \item{quantify}{\code{signature(object = "MSnExp")}: Performs
%%       quantification for all the MS2 spectra of the \code{MSnExp}
%%       instance. See \code{\link{quantify}} documentation for more
%%       details. }
%%     \item{removePeaks}{\code{signature(object = "MSnExp")}: Removes
%%       peaks lower that a threshold \code{t}. See
%%       \code{\link{removePeaks}} documentation for more details and
%%       examples. }
%%     \item{removeReporters}{\code{signature(object = "MSnExp", ...)}:
%%       Removes reporter ion peaks from all MS2 spectra of an
%%       experiment. See \code{\link{removeReporters}} documentation for
%%       more details and examples. }
%%     \item{smooth}{\code{signature(x = "MSnExp")}: Smooths spectra.
%%       See \code{\link{smooth}} documentation for more details and examples. }
%%     \item{addIdentificationData}{\code{signature(object = "MSnExp", ...)}:
%%       Adds identification data to an experiment.
%%       See \code{\link{addIdentificationData}} documentation for
%%       more details and examples. }

%%     \item{removeNoId}{\code{signature(object = "MSnExp", fcol =
%% 	"pepseq", keep = NULL)}: Removes non-identified features. See
%%       \code{\link{removeNoId}} documentation for more details and
%%       examples. }

%%     \item{removeMultipleAssignment}{\code{signature(object = "MSnExp",
%% 	fcol = "nprot")}: Removes protein groups with more than one
%% 	member. The latter is defined by extracting a feature variable
%% 	(default is \code{"nprot"}). }

%%     \item{idSummary}{\code{signature(object = "MSnExp", ...)}:
%%       Prints a summary that lists the percentage of identified features per
%%       file (called \code{coverage}).}

%%     \item{show}{\code{signature(object = "MSnExp")}: Displays object
%%       content as text. }

%%     \item{trimMz}{\code{signature(object = "MSnExp")}: Trims the MZ
%%       range of all the spectra of the \code{MSnExp} instance. See
%%       \code{\link{trimMz}} documentation for more details and
%%       examples. }
%%   }
%% }

%% \references{
%%   Information about the mzXML format as well converters from vendor
%%   specific formats to mzXML:
%%   \url{http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML}.
%% }

\author{
  Johannes Rainer <johannes.rainer@eurac.edu>
}

%% \seealso{
%%   \code{"\linkS4class{pSet}"} and \code{\link{readMSData}} for loading
%%   \code{mzXML}, \code{mzData} or \code{mzML} files to generate an
%%   instance of \code{MSnExp}.
%% }

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

## Calling removePeaks, i.e. setting intensity values below a certain threshold to 0.
## Unlike the name suggests, this is not actually removing peaks. Such peaks with a 0
## intensity are then removed by the "clean" step.
## Also, the manipulations are not applied directly, but put into the "lazy"
## processing queue.
odmse <- removePeaks(odmse)
odmse <- clean(odmse)

## The processing steps are only applied when actual raw data is extracted.
spectra(odmse, scans=c(1, 50))

## Get the polarity of the spectra.
head(polarity(odmse))

## Get the retention time of all spectra
head(rtime(odmse))

}

\keyword{classes}
