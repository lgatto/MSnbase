\name{impute-methods}
\docType{methods}
\alias{impute}
\alias{impute-methods}
\alias{impute,MSnSet-method}

\title{ Quantitative proteomics data imputation }

\description{
  The \code{impute} method performs data imputation on an \code{MSnSet}
  instance. Currently, nearest neighbor averaging is available, as
  implemented in the \code{impute.knn} function. Additional methods
  might be added at a later stage.
  
  The imputation and the parameters are logged into the
  \code{processingData(object)} slot. 
  
  Users should proceed with care when imputing data and take precautions
  to assure that the imputation produce valid results. 
}

\section{Methods}{
  \describe{   
    \item{\code{signature(object = "MSnSet", method = "knn", ...)}}{ 
      This method performs data imputation on the \code{object}
      \code{MSnSet} instance using the \code{method}
      algorithm. Currently, only nearest neighbor averaging (\code{knn})
      from package \code{impute} (see
      \code{?impute::inpute.knn} for details) is available. 
      \code{...} is used to pass parameters to the imputation function.
    }
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

\keyword{methods}

