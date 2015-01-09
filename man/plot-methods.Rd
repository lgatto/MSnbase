\name{plot-methods}
\docType{methods}
\alias{plot,Spectrum-method}
\alias{plot,Spectrum,missing-method}
\alias{plot,Spectrum2,character-method}
\alias{plot-methods}
\alias{plot.MSnExp}
\alias{plot.Spectrum}
\alias{plot.Spectrum.character}
\alias{plot}

\title{ Plotting 'Spectrum' object(s) }

\description{

  These method plot mass spectra MZ values against the intensities. Full
  spectra (using the \code{full} parameter) or specific peaks of
  interest can be plotted using the \code{reporters} parameter. If
  \code{reporters} are specified and \code{full} is set to 'TRUE', a
  sub-figure of the reporter ions is inlaid inside the full spectrum.

  If an \code{"\linkS4class{MSnExp}"} is provided as argument, all the
  spectra are aligned vertically. Experiments can be subset to
  extract spectra of interest using the \code{[} operator or
  \code{\link{extractPrecSpectra}} methods.

  The methods make use the \code{ggplot2} system. An object of class
  'ggplot' is returned invisibly.

  If a single \code{"\linkS4class{Spectrum2}"} and a \code{"character"}
  representing a valid peptide sequence are passed as argument, the
  expected fragement ions are calculated and matched/annotated on the
  spectum plot.
  
}

\arguments{

  \item{x}{ Objects of class \code{"\linkS4class{Spectrum}"},
    \code{"\linkS4class{Spectrum2}"} or \code{"\linkS4class{MSnExp}"} to
    be plotted. }

  \item{y}{ Missing, \code{"\linkS4class{Spectrum}"} or
  \code{"character"}. }
  
  \item{reporters}{ An object of class
    \code{"\linkS4class{ReporterIons}"} that defines the peaks to be
    plotted. If not specified, \code{full} must be set to 'TRUE'.} 
  \item{full}{Logical indicating whether full spectrum (respectively
    spectra) of  only reporter ions of interest should be
    plotted. Default is 'FALSE', in which case \code{reporters} must be
    defined. }
  \item{centroided.}{Logical indicating if spectrum or spectra are in
    centroided mode, in which case peaks are plotted as histograms,
    rather than curves.}
  \item{plot}{Logical specifying whether plot should be printed to
    current device. Default is 'TRUE'.}
  \item{w1}{Width of sticks for full centroided spectra. Default is to
    use maximum MZ value divided by 500. }
  \item{w2}{Width of histogram bars for centroided reporter ions
    plots. Default is 0.01. }

  See below for more details.
  
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "MSnExp", y = "missing", reporters =
	"ReporterIons", full = "logical", plot = "logical")}}{ Plots
      all the spectra in the \code{MSnExp} object vertically. One of
      \code{reporters} must be defined or \code{full} set to 'TRUE'. In
      case of \code{MSnExp} objects, repoter ions are not inlaid when
      \code{full} is 'TRUE'. 
    }
    \item{\code{signature(x = "Spectrum", y = "missing", reporters =
	"ReporterIons", full = "logical", centroided. = "logical", plot =
	"logical", w1, w2)}}{ Displays the MZs against intensities of
      the \code{Spectrum} object as a line plot.  
      At least one of \code{reporters} being defined or \code{full}
      set to 'TRUE' is required. 
      \code{reporters} and \code{full} are used only for
      \code{"\linkS4class{Spectrum2}"}
      objects. Full \code{"\linkS4class{Spectrum1}"} spectra are plotted
      by default.
    }

    \item{\code{signature(x = "Spectrum2", y = "character", orientation
	= "numeric", add = "logical", col = "character", pch, xlab =
	"character", ylab = "character", xlim = "numeric", ylim =
	"numeric", tolerance = "numeric", relative = "logical", type =
	"character", modifications = "numeric", x = "numeric", fragments
	= "function", fragments.cex = "numeric", ... )}}{Plots a single
	MS2 spectrum and annotates the fragment ions based on the
	matching between the peaks in \code{x} and the fragment peaks
	calculated from the peptide sequence \code{y}. The detault
	values are \code{orientation=1}, \code{add=FALSE},
	\code{col="#74ADD1"}, \code{pch=NA}, \code{xlab="m/z"},
	\code{ylab="intensity"}, \code{ylim=c(0, 1)},
	\code{tolerance=0.1}, relative=FALSE, type=c("b", "y"),
	\code{modifications=c(C=160.030649)}, \code{z=1},
	fragments=\code{MSnbase:::calculateFragments_Spectrum2} and
	\code{fragments.cex=0.75}. Additional arguments \code{...} are
	passed to \code{plot.default}.  }

  }  
}


\seealso{
  \code{\link{calculateFragments}} to calculate ions produced by
  fragmentation and \code{\link{plot.Spectrum.Spectrum}} to plot and
  compare 2 spectra and their shared peaks.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk> and Sebastian Gibb
}

\examples{
data(itraqdata)
## plotting experiments
plot(itraqdata[1:2], reporters = iTRAQ4)
plot(itraqdata[1:2], full = TRUE)
## plotting spectra
plot(itraqdata[[1]],reporters = iTRAQ4, full = TRUE)

itraqdata2 <- pickPeaks(itraqdata)
i <- 14
s <- as.character(fData(itraqdata2)[i, "PeptideSequence"])
plot(itraqdata2[[i]], s, main = s)
}

\keyword{methods}
