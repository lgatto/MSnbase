\name{quantify-methods}
\docType{methods}
\alias{quantify-methods}
\alias{quantify}

\title{ Quantifies 'MSnExp' and 'Spectrum' objects }


\description{
  This method quantifies individual \code{"\linkS4class{Spectrum}"}
  objects or full \code{"\linkS4class{MSnExp}"} experiments. Current,
  MS2-level isobar tagging using iTRAQ and TMT (or any arbitrary peaks
  of interest, see \code{"\linkS4class{ReporterIons}"}) and MS2-level
  label-free quantitation (spectral counting, spectral index or spectral
  abundance factor) are available.

  Isobaric tag peaks of single spectra or complete experiments can be
  quantified using appropriate \code{methods}. Label-free quantitation
  is available only for \code{MSnExp} experiments.

  Since version 1.13.5, parallel quantitation is supported by the
  \code{BiocParallel} package and controlled by the \code{BPPARAM}
  argument.
  
}

\arguments{
  \item{object}{
    An instance of class \code{"\linkS4class{Spectrum}"} (isobaric
    tagging only) or \code{"\linkS4class{MSnExp}"}.
  }

  \item{method}{
    Peak quantitation method. For isobaric tags, one of, possibly
    abreviated \code{"trapezoidation"}, \code{"max"}, or
    \code{"sum"}. These methods return respectively the area under the
    peak(s), the maximum of the peak(s) or the sum of all intensities of
    the peak(s).

    For label-free quantitation, one of \code{"SI"} (spectral index),
    \code{"SIgi"} (global intensity spectral index), \code{"SIn"}
    (normalised spectral index), \code{"SAF"} (spectral abundance
    factor) or \code{"NSAF"} (normalised spectral abundance factor).

    Finally, the simple \code{"count"} method counts the occurrence of
    the respective spectra (at this stage all 1s) that can then be used
    as input to \code{\link{combineFeatures}} to implement spectra
    counting. 

  }

  \item{reporters}{
    An instance of class \code{"\linkS4class{ReporterIons}"} that defines
    the peak(s) to be quantified. For isobaric tagging only.
  }

  \item{strict}{
    For isobaric tagging only.  If strict is \code{FALSE} (default), the
    quantitation is performed using data points along the entire width
    of a peak. If strict is set to \code{TRUE}, once the apex(es) is/are
    identified, only data points within apex +/- width of reporter (see
    \code{"\linkS4class{ReporterIons}"}) are used for quantitation.

  }

  \item{BPPARAM}{    
    Support for parallel processing using the \code{BiocParallel}
    infrastructure. When missing (default), the default registered
    \code{BiocParallelParam} parameters are applied using
    \code{bpparam()}. Alternatively, one can pass a valid
    \code{BiocParallelParam} parameter instance: \code{SnowParam},
    \code{MulticoreParam}, \code{DoparParam}, \ldots see the
    \code{BiocParallel} package for details.
  }
  
  \item{parallel}{
    Deprecated. Please see \code{BPPARAM}.   
  }

  \item{qual}{
    Should the \code{qual} slot be populated. Default is \code{TRUE}.
  }

  \item{pepseq}{
    A \code{character} giving the peptide sequence column in the feature
    data. Default is \code{"sequence"}.
  }

  
  \item{verbose}{
    Verbose of the output (only for \code{MSnExp} objects).
  }

  \item{...}{
    Further arguments passed to the quantitation functions.
  }
}

\details{

  \code{"\linkS4class{ReporterIons}"} define specific MZ at which peaks
  are expected and a window around that MZ value.  A peak of interest is
  searched for in that window.  Since version 1.1.2, warnings are not
  thrown anymore in case no data is found in that region or if the peak
  extends outside the window.  This can be checked manually after
  quantitation, by inspecting the quantitation data (using the
  \code{exprs} accessor) for \code{NA} values or by comaring the
  \code{lowerMz} and \code{upperMz} columns in the
  \code{"\linkS4class{MSnSet}"} \code{qual} slot against the respective
  expected \code{mz(reporters)} +/- \code{width(reporters)}.

  Once the range of the curve is found, quantification is performed.  If
  no data points are found in the expected region, \code{NA} is returned
  for the reporter peak MZ.

  Note that for label-free, spectra that have not been identified (the
  corresponding fields in the feature data are populated with \code{NA}
  values) or that have been uniquely assigned to a protein (the
  \code{nprot} feature data is greater that 1) are removed prior to
  quantitation. The latter does not apply for \code{method = "count"}
  but can be applied manually with
  \code{\link{removeMultipleAssignment}}.

}


\section{Methods}{
  \describe{
    \item{
      \code{signature(object = "MSnExp", method = "character", reporters
	= "ReporterIons", verbose = "logical", ...)}}{

      For isobaric tagging, quantifies peaks defined in \code{reporters}
      using \code{method} in all spectra of the \code{MSnExp} object. If
      verbose is set to \code{TRUE}, a progress bar will be displayed.

      For label-free quantitation, the respective quantitation methods
      and normalisations are applied to the spectra. These methods
      require two additional arguments (\code{...}), namely the protein
      accession of identifiers (\code{fcol}, with detault value
      \code{"DatabaseAccess"}) and the protein lengths (\code{plength},
      with default value \code{"DBseqLength"}). These values are
      available of the identification data had been collated using
      \code{\link{addIdentificationData}}.

      An object of class \code{"\linkS4class{MSnSet}"} is returned
      containing the quantified feature expression and all meta data
      inherited from the \code{MSnExp} \code{object} argument.

    }

    \item{\code{signature(object = "Spectrum", method = "character",
	reporters = "ReporterIons")}}{

      Quantifies peaks defined in \code{reporters} using \code{method}
      in the \code{Spectrum} object (isobaric tagging only).

      A list of length 2 will be returned. The first element, named
      \code{peakQuant}, is a 'numeric' of length equal to
      \code{length(reporters)} with quantitation of the reporter peaks
      using \code{method}.

      The second element, names \code{curveStats}, is a 'data.frame' of
      dimension \code{length(reporters)} times 7 giving, for each
      reporter curve parameters: maximum intensity (\code{'maxInt'}),
      number of maxima (\code{'nMaxInt'}), number of data points defined
      the curve (\code{'baseLength'}), lower and upper MZ values for the
      curve (\code{'lowerMz'} and \code{'upperMz'}), reporter
      (\code{'reporter'}) and precursor MZ value (\code{'precursor'})
      when available.

    }
  }
}

\references{

  For details about the spectral index (SI), see Griffin NM, Yu J, Long
  F, Oh P, Shore S, Li Y, Koziol JA, Schnitzer JE.  \emph{Label-free,
  normalized quantification of complex mass spectrometry data for
  proteomic analysis}. Nat Biotechnol. 2010 Jan;28(1):83-9. doi:
  10.1038/nbt.1592. PMID: 20010810; PubMed Central PMCID: PMC2805705.

  For details about the spectra abundance factor, see Paoletti AC,
  Parmely TJ, Tomomori-Sato C, Sato S, Zhu D, Conaway RC, Conaway JW,
  Florens L, Washburn MP. \emph{Quantitative proteomic analysis of
  distinct mammalian Mediator complexes using normalized spectral
  abundance factors}. PNAS. 2006 Dec 12;103(50):18928-33. PMID:
  17138671; PubMed Central PMCID: PMC1672612.

}

\author{
  Laurent Gatto <lg390@cam.ac.uk> and
  Sebastian Gibb <mail@sebastiangibb.de>
}

\examples{

## Quantifying a full experiment using iTRAQ4-plex tagging
data(itraqdata)
msnset <- quantify(itraqdata, method = "trap", reporters = iTRAQ4)
msnset

## specifying a custom parallel framework
## bp <- MulticoreParam(2L) # on Linux/OSX
## bp <- SnowParam(2L) # on Windows
## quantify(itraqdata[1:10], method = "trap", iTRAQ4, BPPARAM = bp)

## Checking for non-quantified peaks
sum(is.na(exprs(msnset)))

## Quantifying a single spectrum
qty <- quantify(itraqdata[[1]], method = "trap", iTRAQ4[1])
qty$peakQuant
qty$curveStats


## Label-free quantitation
## Raw (mzXML) and identification (mzid) files
quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")

msexp <- readMSData(quantFile)
msexp <- addIdentificationData(msexp, identFile)
fData(msexp)$DatabaseAccess

si <- quantify(msexp, method = "SIn")
processingData(si)
exprs(si)

saf <- quantify(msexp, method = "NSAF")
processingData(saf)
exprs(saf)

}


\keyword{methods}
