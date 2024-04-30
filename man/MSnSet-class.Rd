\name{MSnSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnSet-class}
\alias{class:MSnSet}
\alias{MSnSet}


\alias{acquisitionNum,MSnSet-method}
\alias{exprs,MSnSet-method}
\alias{dim,MSnSet-method}
\alias{fileNames,MSnSet-method}
\alias{msInfo,MSnSet-method}
\alias{processingData,MSnSet-method}
\alias{qual,MSnSet-method}
\alias{qual}
\alias{show,MSnSet-method}
\alias{purityCorrect,MSnSet-method}
\alias{purityCorrect,MSnSet,matrix-method}
\alias{meanSdPlot,MSnSet-method}
\alias{t.MSnSet}
\alias{[,MSnSet-method}
\alias{[,MSnSet,ANY,ANY-method}
\alias{[,MSnSet,ANY,ANY,ANY-method}
\alias{as.ExpressionSet.MSnSet}
\alias{as.MSnSet.ExpressionSet}
\alias{as.data.frame.MSnSet}
\alias{ms2df}
\alias{coerce,IBSpectra,MSnSet-method}
\alias{coerce,MSnSet,ExpressionSet-method}
\alias{coerce,ExpressionSet,MSnSet-method}
\alias{coerce,MSnSet,data.frame-method}
\alias{coerce,MSnSet,SummarizedExperiment-method}
\alias{coerce,SummarizedExperiment,MSnSet-method}
\alias{addMSnSetMetadata}
\alias{write.exprs}
\alias{write.exprs,MSnSet-method}
\alias{experimentData<-,MSnSet,MIAPE-method}
\alias{pData<-,MSnSet,data.frame-method}
\alias{fData<-,MSnSet,data.frame-method}
\alias{combine,MSnSet,MSnSet-method}
\alias{topN,MSnSet,MSnSet-method}
\alias{topN,MSnSet-method}
\alias{topN,matrix-method}
\alias{topN}
\alias{filterNA,MSnSet-method}
\alias{filterNA,matrix-method}
\alias{filterNA}
\alias{filterZero,MSnSet-method}
\alias{filterZero,matrix-method}
\alias{filterZero}
\alias{filterMsLevel,MSnSet-method}
\alias{log,MSnSet-method}
\alias{image,MSnSet-method}
\alias{image2}
\alias{MAplot,MSnSet-method}
\alias{addIdentificationData,MSnSet,character-method}
\alias{addIdentificationData,MSnSet,mzIDClasses-method}
\alias{addIdentificationData,MSnSet,mzID-method}
\alias{addIdentificationData,MSnSet,mzRident-method}
\alias{addIdentificationData,MSnSet,mzIDCollection-method}
\alias{addIdentificationData,MSnSet,data.frame-method}
\alias{removeNoId,MSnSet-method}
\alias{removeMultipleAssignment-method}
\alias{removeMultipleAssignment,MSnSet-method}
\alias{removeMultipleAssignment}
\alias{idSummary,MSnSet-method}
\alias{idSummary}
\alias{fromFile,MSnSet-method}

\alias{trimws}
\alias{trimws,MSnSet-method}
\alias{trimws,data.frame-method}

\alias{exptitle,MSnSet-method}
\alias{expemail,MSnSet-method}
\alias{ionSource,MSnSet-method}
\alias{analyser,MSnSet-method}
\alias{analyzer,MSnSet-method}
\alias{detectorType,MSnSet-method}
\alias{description,MSnSet-method}

\alias{updateFvarLabels}
\alias{updateSampleNames}
\alias{updateFeatureNames}

\alias{droplevels.MSnSet}

\title{
  The "MSnSet" Class for MS Proteomics Expression Data and Meta-Data
}

\description{
  The \code{MSnSet} holds quantified expression data for MS proteomics
  data and the experimental meta-data.
  The \code{MSnSet} class is derived from the
  \code{"\linkS4class{eSet}"} class and mimics the
  \code{"\linkS4class{ExpressionSet}"} class classically used for
  microarray data.
}

\section{Objects from the Class}{

  The constructor \code{MSnSet(exprs, fData, pData)} can be used to
  create \code{MSnSet} instances. Argument \code{exprs} is a
  \code{matrix} and \code{fData} and \code{pData} must be of class
  \code{data.frame} or \code{"\linkS4class{AnnotatedDataFrame}"} and all
  must meet the dimensions and name validity constrains.

  Objects can also be created by calls of the form \code{new("MSnSet",
  exprs, ...)}.  See also \code{"\linkS4class{ExpressionSet}"} for
  helpful information.  Expression data produced from other softwares
  can thus make use of this standardized data container to benefit
  \code{R} and \code{Bioconductor} packages. Proteomics expression data
  available as spreadsheets, as produced by third-party software such as
  Proteome Discoverer, MaxQuant, ... can be imported using the
  \code{\link{readMSnSet}} and \code{\link{readMSnSet2}} functions.

  Coercion methods are also available to transform \code{MSnSet} objects
  to \code{IBSpectra}, to \code{data.frame} and to/from
  \code{ExpressionSet} and \code{SummarizedExperiment} objects. In the
  latter case, the metadata available in the \code{protocolData},
  \code{experimentData} are completely dropped, and only the logging
  information of the \code{processingData} slot is retained. All these
  metadata can be subsequently be added using the
  \code{addMSnSetMetadata} (see examples below). When converting a
  \code{SummarizedExperiment} to an \code{MSnSet}, the respective
  metadata slots will be populated if available in the
  \code{SummarizedExperiment} metadata.

  In the frame of the \code{MSnbase} package, \code{MSnSet} instances
  can be generated from \code{"\linkS4class{MSnExp}"} experiments using
  the \code{\link{quantify}} method).

}

\section{Slots}{

  \describe{

    \item{\code{qual}:}{Object of class \code{"data.frame"} that records
      peaks data for each of the reporter ions to be used as quality
      metrics. }

    \item{\code{processingData}:}{Object of class
      \code{"\linkS4class{MSnProcess}"} that records all processing. }

    \item{\code{assayData}:}{Object of class \code{"assayData"}
      containing a matrix with equal with column number equal to
      \code{nrow(phenoData)}. \code{assayData} must contain a matrix
      \code{exprs} with rows represening features (e.g., reporters ions)
      and columns representing samples. See the \code{"\linkS4class{AssayData}"}
      class, \code{\link{exprs}} and \code{\link{assayData}} accessor
      for more details. This slot in indirectly inherited from
      \code{"\linkS4class{eSet}"}. }

    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      containing experimenter-supplied variables describing sample (i.e
      the individual tags for an labelled MS experiment) (indireclty
      inherited from \code{"\linkS4class{eSet}"}). See
      \code{\link{phenoData}} and the \code{"\linkS4class{eSet}"} class
      for more details. This slot can be accessed as a \code{data.frame}
      with \code{pData} and be replaced by a new valid (i.e. of
      compatible dimensions and row names) \code{data.frame} with
      \code{pData()<-}.  }

    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"} containing variables describing
      features (spectra in our case), e.g. identificaiton data, peptide
      sequence, identification score,... (inherited indirectly from
      \code{"\linkS4class{eSet}"}). See \code{\link{featureData}} and
      the \code{"\linkS4class{eSet}"} class for more details. This slot
      can be accessed as a \code{data.frame} with \code{fData} and be
      replaced by a new valid (i.e. of compatible dimensions and row
      names) \code{data.frame} with \code{fData()<-}.
    }

    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAPE}"}, containing details of experimental
      methods (inherited from \code{"\linkS4class{eSet}"}). See
      \code{\link{experimentData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }

    \item{\code{annotation}:}{ not used here. }

    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited indirectly from
      \code{"\linkS4class{eSet}"}). See
      \code{\link{protocolData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }

    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of R,
      the Biobase package, \code{"\linkS4class{eSet}"},
      \code{"\linkS4class{pSet}"} and \code{MSnSet} of the
      current instance. Intended for developer use and debugging (inherited
      indirectly from \code{"\linkS4class{eSet}"}). }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{eSet}"}, directly.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
  Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}

\section{Methods}{
  MSnSet specific methods or over-riding it's super-class are described
  below. See also more \code{"\linkS4class{eSet}"} for
  inherited methods.

  \describe{
    \item{acquisitionNum}{
      \code{acquisitionNum(signature(object = "MSnSet"))}: Returns the
      a numeric vector with acquisition number of each spectrum. The vector
      names are the corresponding spectrum names.
      The information is extracted from the object's \code{featureData} slot.
    }

    \item{fromFile}{
      \code{fromFile(signature(object = "MSnSet"))}: get the index of
      the file (in \code{fileNames(object)}) from which the raw
      spectra from which the corresponding feature were originally
      read. The relevant information is extracted from the object's
      \code{featureData} slot.

      Returns a numeric vector with names corresponding to the spectrum
      names.
    }

    \item{dim}{\code{signature(x = "MSnSet")}: Returns the dimensions of
      object's assay data, i.e the number of samples and the number of
      features. }

    \item{fileNames}{\code{signature(object = "MSnSet")}: Access file
      names in the \code{processingData} slot. }

    \item{msInfo}{\code{signature(object = "MSnSet")}: Prints the MIAPE-MS
      meta-data stored in the \code{experimentData} slot. }

    \item{processingData}{\code{signature(object = "MSnSet")}: Access the
      \code{processingData} slot. }

    \item{show}{\code{signature(object = "MSnSet")}: Displays object
      content as text. }

    \item{qual}{\code{signature(object = "MSnSet")}: Access the reporter
      ion peaks description. }

    \item{purityCorrect}{\code{signature(object = "MSnSet", impurities =
	"matrix")}: performs reporter ions purity correction. See
      \code{\link{purityCorrect}} documentation for more details. }

    \item{normalise}{\code{signature(object = "MSnSet")}: Performs
      \code{MSnSet} normalisation. See \code{\link{normalise}} for more
      details. }

    \item{t}{\code{signature(x = "MSnSet")}: Returns a transposed
      \code{MSnSet} object where features are now aligned along columns
      and samples along rows and the \code{phenoData} and
      \code{featureData} slots have been swapped. The
      \code{protocolData} slot is always dropped. }

    \item{as(,"ExpressionSet")}{\code{signature(x = "MSnSet")}: Coerce
      object from \code{MSnSet} to
      \code{\link{ExpressionSet-class}}. The \code{experimentData} slot
      is converted to a \code{MIAME} instance. It is also possible to
      coerce an \code{ExpressionSet} to and \code{MSnSet}, in which case
      the \code{experimentData} slot is newly initialised. }

    \item{as(,"IBSpectra")}{\code{signature(x = "MSnSet")}: Coerce
      object from \code{MSnSet} to \code{IBSpectra} from the
      \code{isobar} package.
    }

    \item{as(,"data.frame")}{\code{signature(x = "MSnSet")}: Coerce
      object from \code{MSnSet} to \code{data.frame}. The \code{MSnSet}
      is transposed and the \code{PhenoData} slot is appended. }

    \item{as(,"SummarizedExperiment")}{\code{signature(x = "MSnSet")}:
      Coerce object from \code{MSnSet} to
      \code{SummarizedExperiment}. Only part of the metadata is
      retained. See \code{addMSnSetMetadata} and the example below for
      details.}

    \item{write.exprs}{\code{signature(x = "MSnSet")}: Writes expression values
      to a tab-separated file (default is \code{tmp.txt}). The
      \code{fDataCols} parameter can be used to specify which
      \code{featureData} columns (as column names, column number or
      \code{logical}) to append on the right of the expression matrix.
      The following arguments are the same as \code{write.table}.}

    \item{combine}{\code{signature(x = "MSnSet", y = "MSnSet", ...)}: Combines
      2 or more \code{MSnSet} instances according to their feature names.
      Note that the \code{qual} slot and the processing information are
      silently dropped. }

    \item{topN}{\code{signature(object = "MSnSet", groupBy, n = 3, fun,
	..., verbose = isMSnbaseVerbose())}:
      Selects the \code{n} most intense features (typically peptides or
      spectra) out of all available for each set defined by
      \code{groupBy} (typically proteins) and creates a new instance of
      class \code{MSnSet}. If less than \code{n} features are available,
      all are selected. The \code{ncol(object)} features are summerised
      using \code{fun} (default is \code{sum}) prior to be ordered in
      decreasing order. Additional parameters can be passed to
      \code{fun} through \code{...}, for instance to control the
      behaviour of \code{topN} in case of \code{NA} values.
      (Works also with \code{matrix} instances.)

      See also the \code{\link{nQuants}} function to retrieve the
      actual number of retained peptides out of \code{n}.

      A complete use case using \code{topN} and \code{nQuants} is
      detailed in the \code{synapter} package vignette.
    }

    \item{filterNA}{\code{signature(object = "MSnSet", pNA = "numeric",
      pattern = "character", droplevels = "logical")}: This method
      subsets \code{object} by removing features that have (strictly)
      more than \code{pNA} percent of NA values. Default \code{pNA} is
      0, which removes any feature that exhibits missing data.
      The method can also be used with a character pattern composed of
      \code{0} or \code{1} characters only. A \code{0} represent a
      column/sample that is allowed a missing values, while
      columns/samples with and \code{1} must not have \code{NA}s.

      This method also accepts \code{matrix}
      instances. \code{droplevels} defines whether unused levels in the
      feature meta-data ought to be lost. Default is \code{TRUE}.
      See the \code{droplevels} method below.

      See also the \code{\link{is.na.MSnSet}} and \code{\link{plotNA}}
      methods for missing data exploration.
    }

    \item{filterZero}{\code{signature(object = "MSnSet", pNA = "numeric",
      pattern = "character", droplevels = "logical")}: As
      \code{filterNA}, but for zeros.
    }

    \item{filterMsLevel}{\code{signature(object = "MSnSet", msLevel. =
      "numeric", fcol = "character")} Keeps only spectra with level
      \code{msLevel.}, as defined by the \code{fcol} feature variable
      (default is \code{"msLevel"}).
    }

    \item{log}{\code{signature(object = "MSnSet", base = "numeric")} Log
      transforms \code{exprs(object)} using
      \code{base::log}. \code{base} (defaults is \code{e='exp(1)'}) must
      be a positive or complex number, the base with respect to which
      logarithms are computed.
    }

    \item{droplevels}{\code{signature(x = "MSnSet", ...)}Drops the unused
      factor levels in the \code{featureData} slot. See
      \code{\link{droplevels}} for details.
    }

    \item{impute}{\code{signature(object = "MSnSet", ...)}
      Performs data imputation on the \code{MSnSet} object.
      See \code{\link{impute}} for more details.
    }

    \item{trimws}{\code{signature(object = "MSnSet", ...)}Trim leading and/or
      trailing white spaces in the feature data slot. Also available for
      \code{data.frame} objects. See \code{?base::\link[base]{trimws}}
      for details.
    }

  }

  Additional accessors for the experimental metadata
  (\code{experimentData} slot) are defined. See
  \code{"\linkS4class{MIAPE}"} for details.
}

\section{Plotting}{
  \describe{
    \item{meanSdPlot}{\code{signature(object = "MSnSet")} Plots row
      standard deviations versus row means. See
      \code{\link{meanSdPlot}} (\code{vsn} package) for more details.
    }

    \item{image}{\code{signature(x = "MSnSet", facetBy = "character",
	sOrderBy = "character", legend = "character", low = "character",
	high = "character", fnames = "logical", nmax =
	"numeric")} Produces an heatmap of expression values in the
	\code{x} object. Simple horizontal facetting is enabled by
	passing a single character as \code{facetBy}. Arbitrary
	facetting can be performed manually by saving the return value
	of the method (see example below). Re-ordering of the samples is
	possible by providing the name of a phenotypic variable to
	\code{sOrderBy}. The title of the legend can be set with
	\code{legend} and the colours with the \code{low} and
	\code{high} arguments. If any negative value is detected in the
	data, the values are considered as log fold-changes and a
	divergent colour scale is used. Otherwise, a gradient from low
	to high is used. To scale the quantitative data in \code{x}
	prior to plotting, please see the \code{scale} method.

	When there are more than \code{nmax} (default is 50)
	features/rows, these are not printed. This behaviour can be
	controlled by setting \code{fnames} to \code{TRUE} (always
	print) or \code{FALSE} (never print). See examples below.

	The code is based on Vlad Petyuk's
	\code{vp.misc::image_msnset}. The previous version of this
	method is still available through the \code{image2} function.
      }

    \item{plotNA}{\code{signature(object = "MSnSet", pNA =
	"numeric")}
	Plots missing data for an \code{MSnSet} instance. \code{pNA} is a
	\code{numeric} of length 1 that specifies the percentage
	of accepted missing data values per features. This value will be
	highlighted with a point on the figure, illustrating the overall
	percentage of NA values in the full data set and the number of
	proteins retained. Default is 1/2. See also
	\code{\link{plotNA}}.
      }

    \item{MAplot}{\code{signature(object = "MSnSet", log.it = "logical",
	base = "numeric", ...)}
      Produces MA plots (Ratio as a function
      of average intensity) for the samples in \code{object}. If
      \code{ncol(object) == 2}, then one MA plot is produced using the
      \code{\link{ma.plot}} function from the \code{affy} package. If
      \code{object} has more than 2 columns, then
      \code{\link{mva.pairs}}. \code{log.it} specifies is the data
      should be log-transformed (default is \code{TRUE}) using
      \code{base}. Further \code{...} arguments will be passed to the
      respective functions.
    }

    \item{addIdentificationData}{\code{signature(object = "MSnSet", ...)}:
      Adds identification data to a \code{MSnSet} instance.
      See \code{\link{addIdentificationData}} documentation for
      more details and examples.
    }

    \item{removeNoId}{\code{signature(object = "MSnSet", fcol =
	"pepseq", keep = NULL)}: Removes non-identified features. See
      \code{\link{removeNoId}} documentation for more details and
      examples.
    }

    \item{removeMultipleAssignment}{\code{signature(object = "MSnSet",
	fcol = "nprot")}: Removes protein groups (or feature belong to
	protein groups) with more than one member. The latter is defined
	by extracting a feature variable (default is
	\code{"nprot"}). Also removes non-identified features.
      }

    \item{idSummary}{\code{signature(object = "MSnSet", ...)}: Prints a
      summary that lists the percentage of identified features per file
      (called \code{coverage}).
    }

  }
}

\section{Functions}{
  \describe{

    \item{updateFvarLabels}{\code{signature(object, label, sep)} This
      function updates \code{object}'s featureData variable labels by
      appending \code{label}. By default, \code{label} is the variable
      name and the separator \code{sep} is \code{.}.
    }

    \item{updateSampleNames}{\code{signature(object, label, sep)} This
      function updates \code{object}'s sample names by appending
      \code{label}. By default, \code{label} is the variable name and
      the separator \code{sep} is \code{.}.}

    \item{updateFeatureNames}{\code{signature(object, label, sep)} This
      function updates \code{object}'s feature names by appending
      \code{label}. By default, \code{label} is the variable name and
      the separator \code{sep} is \code{.}.}

    \item{ms2df}{\code{signature(x, fcols)} Coerces the \code{MSnSet} instance
      to a \code{data.frame}. The direction of the data is retained and
      the feature variable labels that match \code{fcol} are appended to
      the expression values. See also \code{as(x, "data.frame")} above.
      }

    \item{addMSnSetMetadata}{\code{signature(x, y)} When coercing an
      \code{MSnSet} \code{y} to a \code{SummarizedExperiment} \code{x}
      with \code{x <- as(y, "SummarizedExperiment")}, most of \code{y}'s
      metadata is lost. Only the file names, the processing log and the
      MSnbase version from the \code{processingData} slots are passed
      along. The \code{addMSnSetMetadata} function can be used to add
      the complete \code{processingData}, \code{experimentData} and
      \code{protocolData} slots. The downside of this is that MSnbase is
      now required to use the \code{SummarizedExperiment} object.
    }

  }
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{

  \code{"\linkS4class{eSet}"}, \code{"\linkS4class{ExpressionSet}"} and
  \code{\link{quantify}}. \code{MSnSet} quantitation values and
  annotation can be exported to a file with
  \code{\link{write.exprs}}. See \code{\link{readMSnSet}} to
  create and \code{MSnSet} using data available in a spreadsheet or
  \code{data.frame}.

}

\examples{
data(msnset)
msnset <- msnset[10:15]

exprs(msnset)[1, c(1, 4)] <- NA
exprs(msnset)[2, c(1, 2)] <- NA
is.na(msnset)
featureNames(filterNA(msnset, pNA = 1/4))
featureNames(filterNA(msnset, pattern = "0110"))

M <- matrix(rnorm(12), 4)
pd <- data.frame(otherpdata = letters[1:3])
fd <- data.frame(otherfdata = letters[1:4])
x0 <- MSnSet(M, fd, pd)
sampleNames(x0)

M <- matrix(rnorm(12), 4)
colnames(M) <- LETTERS[1:3]
rownames(M) <- paste0("id", LETTERS[1:4])
pd <- data.frame(otherpdata = letters[1:3])
rownames(pd) <- colnames(M)
fd <- data.frame(otherfdata = letters[1:4])
rownames(fd) <- rownames(M)
x <- MSnSet(M, fd, pd)
sampleNames(x)


## Visualisation

library("pRolocdata")
data(dunkley2006)
image(dunkley2006)
## Changing colours
image(dunkley2006, high = "darkgreen")
image(dunkley2006, high = "darkgreen", low = "yellow")
## Forcing feature names
image(dunkley2006, fnames = TRUE)
## Facetting
image(dunkley2006, facetBy = "replicate")
p <- image(dunkley2006)
library("ggplot2") ## for facet_grid
p + facet_grid(replicate ~ membrane.prep, scales = 'free', space = 'free')
p + facet_grid(markers ~ replicate)
## Fold-changes
dd <- dunkley2006
exprs(dd) <- exprs(dd) - 0.25
image(dd)
image(dd, low = "green", high = "red")
## Feature names are displayed by default for smaller data
dunkley2006 <- dunkley2006[1:25, ]
image(dunkley2006)
image(dunkley2006, legend = "hello")


## Coercion
if (require("SummarizedExperiment")) {
    data(msnset)
    se <- as(msnset, "SummarizedExperiment")
    metadata(se) ## only logging
    se <- addMSnSetMetadata(se, msnset)
    metadata(se) ## all metadata
    msnset2 <- as(se, "MSnSet")
    processingData(msnset2)
}

as(msnset, "ExpressionSet")

}

\keyword{classes}
