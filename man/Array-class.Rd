\name{Array-class}
\docType{class}

\alias{class:Array}
\alias{Array-class}
\alias{Array}

\alias{length,Array-method}
\alias{isEmpty,Array-method}
\alias{[[,Array-method}

\alias{as.matrix.Array}
\alias{as.matrix,Array-method}

\alias{t.Array}
\alias{t,Array-method}

\alias{showAsCell,Array-method}

\title{Array objects}

\description{
  Array is a virtual class with no slots intended to be extended by
  concrete subclasses with an array-like semantic.
}

\details{
  Some examples of Array derivatives:
  \itemize{
    \item \link[SparseArray]{SparseArray} objects implemented
          in the \pkg{SparseArray} package.

    \item \link[DelayedArray]{DelayedArray} objects implemented
          in the \pkg{DelayedArray} package.

    \item \link{ArrayGrid} and \link{ArrayViewport} objects implemented
          in this package (the \pkg{S4Arrays} package).
  }
}

\seealso{
  \itemize{
    \item \link[base]{array} and \link[base]{matrix} objects in base R.
  }
}

\examples{
showClass("Array")  # virtual class with no slots
}
\keyword{internal}
