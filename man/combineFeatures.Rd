\name{combineFeatures}
\alias{combineFeatures}

\title{ Combines features in an 'MSnSet' object }


\description{
  This function combines the features in an
  \code{"\linkS4class{MSnSet}"} instance applying a summarisation
  function (see \code{fun} argument) to sets of features as
  defined by a factor (see \code{groupBy} argument). 
}

\usage{
combineFeatures(object, groupBy, fun = c("mean", "median", "weighted.mean", "sum", "medpolish"), ..., verbose = TRUE)
}

\arguments{
  \item{object}{ An instance of class \code{"\linkS4class{MSnSet}"}
    whose features will be summerised. }
  \item{groupBy}{ An object of class factor defining how to summerise
    the features. }
  \item{fun}{ The summerising function. Currently, mean, median,
    weighted mean, sum and median polish are implemented, but
    user-defined functions can also be supplied. }
  \item{\dots}{ Additional arguments for the \code{fun} function. }
  \item{verbose}{ A \code{logical} indicating whether verbose output is
    to be printed out. }
}

% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }

\value{
  A new \code{"\linkS4class{MSnSet}"} instance is returned with
  \code{ncol} (i.e. number of samples) is unchanged, but \code{nrow}
  (i.e. the number od features) is now equals to the number of levels in
  \code{groupBy}. The feature metadata (\code{featureData} slot) is
  updated accordingly and only the first occurrence of a feature in
  the original feature meta-data is kept. 
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }

\examples{
file <- dir(system.file(package="MSnbase",dir="extdata"),
            full.name=TRUE,
            pattern="mzXML$")
aa <- readMzXMLData(file)
aa.quant <- quantify(aa,method="max",reporters=iTRAQ4)
dim(aa.quant)
## arbitrary grouping into two groups
grp <- as.factor(c(1,1,2,2,2))
aa.quant.comb <- combineFeatures(aa.quant,grp,"sum")
dim(aa.quant.comb)
exprs(aa.quant.comb)
}

%% \keyword{ ~kwd1 }
