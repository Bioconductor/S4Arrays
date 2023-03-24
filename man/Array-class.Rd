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
  Array is a virtual class intended to be extended by concrete subclasses
  with an array-like semantic.
}

\details{
  Some examples of Array derivatives:
  \itemize{
    \item The \link[SparseArray]{SparseArray} class defined in the
          \pkg{SparseArray} package.

    \item The \link[DelayedArray]{DelayedArray} class defined in the
          \pkg{DelayedArray} package.

    \item The \link{ArrayGrid} and \link{ArrayViewport} classes defined in
          this package (the \pkg{S4Arrays} package).
  }
}

\examples{
showClass("Array")  # virtual class with no slots
}
\keyword{internal}
