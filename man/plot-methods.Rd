\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot.MSnExp}
\alias{plot.Spectrum}
\alias{plot}

\title{ Plotting 'Spectrum' object(s) }

\description{
  These method plot mass spectra MZ values against the intensities as
  line plots. Full spectra (using the \code{full} parameter) or specific
  peaks of interest can be plotted using the \code{reporters}
  parameter. If \code{reporters} are specified and \code{full} is set to
  'TRUE', a sub-figure of the reporter ions is inlaid inside the full
  spectrum.

  If an \code{"\linkS4class{MSnExp}"} is provided as argument, all the
  spectra are aligned vertically. Experiments can be subset to
  extract spectra of interest using the \code{[} operator or
  \code{\link{extractPrecSpectra}} methods.

  The methods make use the \code{ggplot2} system. An object of class
  'ggplot' is returned invisibly. 
}

\arguments{
  \item{x}{ Objects of class \code{"\linkS4class{Spectrum}"} or
    \code{"\linkS4class{MSnExp}"} to be plotted. }
  \item{y}{ Not used in these methods. }
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
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
data(itraqdata)
## plotting experiments
plot(itraqdata[1:2], reporters = iTRAQ4)
plot(itraqdata[1:2], full = TRUE)
## plotting spectra
plot(itraqdata[[1]],reporters = iTRAQ4, full = TRUE)
}

\keyword{methods}
