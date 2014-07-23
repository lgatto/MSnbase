\name{impute-methods}
\docType{methods}
\alias{impute}
\alias{impute-methods}
\alias{impute,MSnSet-method}

\title{ Quantitative proteomics data imputation }

\description{ The \code{impute} method performs data imputation on an
  \code{MSnSet} instance. Currently, nearest neighbor averaging
  (\code{knn}) and Bayesian missing value imputation (\code{bpca},
  default) are available, as implemented in the
  \code{impute::impute.knn} and \code{pcaMethods::pca} functions
  respectively. Additional methods might be added at a later stage.
  
  The imputation and the parameters are logged into the
  \code{processingData(object)} slot. 
  
  Users should proceed with care when imputing data and take precautions
  to assure that the imputation produce valid results. 
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnSet",
	method = c("bpca", "knn"), ...)}}{ This method performs data
      imputation on the \code{object} \code{MSnSet} instance using the
      \code{method} algorithm. \code{...} is used to pass parameters to the
      imputation function. }
  } 
}

\examples{
data(itraqdata)
qnt <- quantify(itraqdata, reporters = iTRAQ4)
sum(is.na(qnt))
iqnt <- impute(qnt)
sum(is.na(iqnt))
processingData(iqnt)
}

\references{

Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor
Hastie, Robert Tibshirani, David Botstein and Russ B.  Altman, Missing
value estimation methods for DNA microarrays Bioinformatics (2001) 17
(6): 520-525.

Oba et al., A Bayesian missing value estimation method for gene
expression profile data, Bioinformatics (2003) 19 (16): 2088-2096.
}

\author{
  Laurent Gatto, Samuel Wieczorek and Vasile-Cosmin Lazar
}

\keyword{methods}

