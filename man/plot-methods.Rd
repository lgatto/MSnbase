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

\title{ Plotting 'MSnExp' and 'Spectrum' object(s) }

\description{

  These methods provide the functionality to plot mass spectrometry data
  provided as \code{\linkS4class{MSnExp}},
  \code{\linkS4class{OnDiskMSnExp}} or \code{\linkS4class{Spectrum}}
  objects. Most functions plot mass spectra M/Z values against
  intensities.

  Full spectra (using the \code{full} parameter) or specific peaks of
  interest can be plotted using the \code{reporters} parameter. If
  \code{reporters} are specified and \code{full} is set to 'TRUE', a
  sub-figure of the reporter ions is inlaid inside the full spectrum.

  If an \code{"\linkS4class{MSnExp}"} is provided as argument, all the
  spectra are aligned vertically. Experiments can be subset to
  extract spectra of interest using the \code{[} operator or
  \code{\link{extractPrecSpectra}} methods.

  Most methods make use the \code{ggplot2} system in which case an
  object of class 'ggplot' is returned invisibly.

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
    \item{\code{plot(signature(x = "MSnExp", y = "missing"),
	type = c("spectra", "XIC"), reporters = "ReporterIons",
	full = "logical", plot = "logical", ...)}}{

      For \code{type = "spectra"}: Plots all the spectra in the
      \code{MSnExp} object vertically. One of \code{reporters} must be
      defined or \code{full} set to 'TRUE'. In case of \code{MSnExp}
      objects, repoter ions are not inlaid when \code{full} is 'TRUE'.

      For \code{type = "XIC"}: Plots a combined plot of retention time
      against m/z values and retention time against largest signal per
      spectrum for each file. Data points are colored by intensity. The
      lower part of the plot represents the location of the individual
      signals in the retention time - m/z space, the upper part the base
      peak chromatogram of the data (i.e. the largest signal for each
      spectrum). This plot type is restricted to MS level 1 data and is
      most useful for LC-MS data.
      Ideally, the \code{MSnExp} (or \code{OnDiskMSnExp})
      object should be filtered first using the \code{\link{filterRt}}
      and \code{\link{filterMz}} functions to narrow on an ion of
      interest. See examples below. This plot uses base R
      plotting. Additional arguments to the \code{plot} function can be
      passed with \code{...}.
      
      Additional arguments for \code{type = "XIC"} are:
      \describe{
	\item{\code{col}}{color for the border of the points. Defaults to
	  \code{col = "grey"}.}

	\item{\code{colramp}}{color function/ramp to be used for the
	  intensity-dependent background color of data points. Defaults
	  to \code{colramp = topo.colors}.}

	\item{\code{grid.color}}{color for the grid lines. Defaults to
	  \code{grid.color = "lightgrey"}; use \code{grid.color = NA} to
	  disable grid lines altogether.}

	\item{\code{pch}}{point character. Defaults to \code{pch = 21}}.

	\item{\code{...}}{additional parameters for the low-level
	  \code{plot} function.}
	
      }
    }
    \item{\code{plot(signature(x = "Spectrum", y = "missing"), reporters =
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

    \item{\code{plot(signature(x = "Spectrum2", y = "character"), orientation
	= "numeric", add = "logical", col = "character", pch, xlab =
	"character", ylab = "character", xlim = "numeric", ylim =
	"numeric", tolerance = "numeric", relative = "logical", type =
	"character", modifications = "numeric", x = "numeric", fragments
	= "data.frame", fragments.cex = "numeric", ... )}}{Plots a single
	MS2 spectrum and annotates the fragment ions based on the
	matching between the peaks in \code{x} and the fragment peaks
	calculated from the peptide sequence \code{y}. The default
	values are \code{orientation=1}, \code{add=FALSE},
	\code{col="#74ADD1"}, \code{pch=NA}, \code{xlab="m/z"},
	\code{ylab="intensity"}, \code{ylim=c(0, 1)},
	\code{tolerance=25e-6}, relative=TRUE, type=c("b", "y"),
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

  \code{\link{Chromatogram}} for plotting of chromatographic data.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>, Johannes Rainer and Sebastian Gibb
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

## Load profile-mode LC-MS files
library(msdata)
od <- readMSData(dir(system.file("sciex", package = "msdata"),
                     full.names = TRUE), mode = "onDisk")
## Restrict the MS data to signal for serine
serine <- filterMz(filterRt(od, rt = c(175, 190)), mz = c(106.04, 106.06))
plot(serine, type = "XIC")

## Same plot but using heat.colors, rectangles and no point border
plot(serine, type = "XIC", pch = 22, colramp = heat.colors, col = NA)

}

\keyword{methods}
