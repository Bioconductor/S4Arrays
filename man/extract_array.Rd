\name{extract_array}

\alias{extract_array}
\alias{extract_array,ANY-method}
\alias{extract_array,array-method}
\alias{extract_array,data.frame-method}
\alias{extract_array,DataFrame-method}

\alias{as.array.Array}
\alias{as.array,Array-method}

\alias{as.data.frame.Array}
\alias{as.data.frame,Array-method}
\alias{as.vector.Array}
\alias{as.vector,Array-method}
\alias{as.logical.Array}
\alias{as.logical,Array-method}
\alias{as.integer.Array}
\alias{as.integer,Array-method}
\alias{as.numeric.Array}
\alias{as.numeric,Array-method}
\alias{as.complex.Array}
\alias{as.complex,Array-method}
\alias{as.character.Array}
\alias{as.character,Array-method}
\alias{as.raw.Array}
\alias{as.raw,Array-method}

\title{extract_array}

\description{
  \code{extract_array} is an internal generic function not aimed to be used
  directly by the user. It has methods defined for array, data.frame,
  \link[S4Vectors]{DataFrame} objects, and other array-like objects.

  Note that \code{extract_array} is part of the \emph{seed contract} as
  defined in the \emph{Implementing A DelayedArray Backend} vignette from
  the \pkg{DelayedArray} package.
}

\usage{
## The extract_array() S4 generic:

extract_array(x, index)

## extract_array() methods defined in the S4Arrays package:

\S4method{extract_array}{ANY}(x, index)

\S4method{extract_array}{array}(x, index)

\S4method{extract_array}{data.frame}(x, index)

\S4method{extract_array}{DataFrame}(x, index)
}

\arguments{
  \item{x}{
    An array-like object.
    This can be an ordinary array, a \link[SparseArray]{SparseArray}
    object from the \pkg{SparseArray} package, a \linkS4class{dgCMatrix}
    object from the \pkg{Matrix} package, a \link[DelayedArray]{DelayedArray}
    object from the \pkg{DelayedArray} package, or any object with an array
    semantic (i.e. an object for which \code{dim(x)} is not NULL).

    Note that data.frame and \link[S4Vectors]{DataFrame} objects are also
    supported.
  }
  \item{index}{
    An unnamed list of subscripts as positive integer vectors, one vector
    per dimension in \code{x}. Empty and missing subscripts (represented
    by \code{integer(0)} and \code{NULL} list elements, respectively) are
    allowed. The subscripts can contain duplicated indices. They cannot
    contain NAs or non-positive values.
  }
}

\details{
  \code{extract_array} methods need to support empty and missing subscripts
  e.g. \code{extract_array(x, list(NULL, integer(0)))} must return an M x 0
  matrix and \code{extract_array(x, list(integer(0), integer(0)))} a 0 x 0
  matrix.

  Also subscripts are allowed to contain duplicated indices so things like
  \code{extract_array(seed, list(c(1:3, 3:1), 2L))} need to be supported.
}

\value{
  An \emph{ordinary} array of the appropriate type (i.e. integer, double,
  etc...). For example, if \code{x} is an object representing an M x N matrix
  of complex numbers, \code{extract_array(x, list(NULL, 2L))} must return
  its 2nd column as an \emph{ordinary} M x 1 matrix of type complex.
}

\seealso{
  \itemize{
    \item \code{\link{read_block}} to read a block of data from an
          array-like object.

    \item \link[base]{array} and \link[base]{data.frame} objects in base R.

    \item \link[SparseArray]{SparseArray} objects implemented in the
          \pkg{SparseArray} package.

    \item \link[DelayedArray]{DelayedArray} objects implemented in the
          \pkg{DelayedArray} package.

    \item \link[S4Vectors]{DataFrame} objects implemented in the
          \pkg{S4Vectors} package.
  }
}

\examples{
## Coming soon...
}
\keyword{internal}
