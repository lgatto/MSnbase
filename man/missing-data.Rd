\name{missing-data}
\alias{missing-data}

\title{Documenting missing data visualisation}

\description{

  There is a need for adequate handling of missing value impuation in
  quantitative proteomics. Before developing a framework to handle
  missing data imputation optimally, we propose a set of visualisation
  tools. This document serves as an internal notebook for current
  progress and ideas that will eventually materialise in exported
  functionality in the \code{MSnbase} package.
}

\details{
  The explore the structure of missing values, we propose to

  1. Explore missing values in the frame of the experimental design. The
  \code{\link{imageNA2}} function offers such a simple visualisation. It
  is currently limited to 2-group designs/comparisons. In case of time
  course experiments or sub-cellular fractionation along a density
  gradient, we propose to split the time/gradient into 2 groups
  (early/late, top/bottom) as a first approximation.

  2. Explore the proportion of missing values in each group.

  3. Explore the total and group-wise feature intensity distributions.

  The existing \code{\link{plotNA}} function illustrates the
  completeness/missingness of the data.
  
}

\seealso{\code{\link{plotNA}}, \code{\link{imageNA2}}.}

\author{Laurent Gatto <lg390@cam.ac.uk>, Vasile-Cosmin Lavar, Samuel
Wieczorek and Thomas Burger}

\keyword{documentation, internal}