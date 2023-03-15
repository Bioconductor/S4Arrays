\name{type}

\alias{type}
\alias{type,ANY-method}
\alias{type,DataFrame-method}

\title{Get the type of the elements of an array-like object}

\description{
  The \pkg{S4Arrays} package defines a couple of \code{type()} methods to get
  the type of the \emph{elements} of a matrix-like or array-like object.
}

\usage{
\S4method{type}{ANY}(x)

\S4method{type}{DataFrame}(x)
}

\arguments{
  \item{x}{
    For the default \code{type()} method: An array-like object.
    This can be an ordinary array, a \link[SparseArray]{SparseArray}
    object from the \pkg{SparseArray} package, a \linkS4class{dgCMatrix}
    object from the \pkg{Matrix} package, a \link[DelayedArray]{DelayedArray}
    object from the \pkg{DelayedArray} package, or any object with an array
    semantic (i.e. an object for which \code{dim(x)} is not NULL).

    For the method for \link[S4Vectors]{DataFrame} objects:
    A \link[S4Vectors]{DataFrame} derivative for which
    \code{as.data.frame(x)} preserves the number of columns.
    See below for more information.
  }
}

\details{
  Note that for an ordinary matrix or array \code{x}, \code{type(x)} is
  the same as \code{typeof(x)}.

  On an array-like object \code{x} that is not an ordinary array,
  \code{type(x)} is \emph{semantically equivalent} to
  \code{typeof(as.array(x))}. However, the actual implementation is
  careful to avoid turning the full array-like object \code{x} into
  an ordinary array, as this would tend to be very inefficient in general.
  For example, doing so on a big \link[DelayedArray]{DelayedArray} object
  could easily eat all the memory available on the machine.

  On a \link[S4Vectors]{DataFrame} object, \code{type(x)} only
  works if \code{as.data.frame(x)} preserves the number of columns,
  in which case it is \emph{semantically equivalent} to
  \code{typeof(as.matrix(as.data.frame(x)))}. Here too, the actual
  implementation is careful to avoid turning the full object into a
  data frame, then into a matrix, for efficiency reasons.
}

\value{
  A single string indicating the type of the array elements in \code{x}.
}

\seealso{
  \itemize{
    \item The \code{\link[BiocGenerics]{type}} generic function defined
          in the \pkg{BiocGenerics} package.

    \item \link[SparseArray]{SparseArray} objects implemented in the
          \pkg{SparseArray} package.

    \item \link[DelayedArray]{DelayedArray} objects implemented in the
          \pkg{DelayedArray} package.

    \item \link[S4Vectors]{DataFrame} objects implemented in the
          \pkg{S4Vectors} package.
  }
}

\examples{
m <- matrix(rpois(54e6, lambda=0.4), ncol=1200)
type(m)   # integer

x1 <- as(m, "dgCMatrix")
type(x1)  # double

library(SparseArray)
x2 <- SparseArray(m)
type(x2)  # integer
}
\keyword{methods}
