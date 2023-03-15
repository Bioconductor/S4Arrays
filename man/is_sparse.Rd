\name{is_sparse}

\alias{is_sparse}
\alias{is_sparse<-}
\alias{is_sparse,ANY-method}
\alias{is_sparse,dgCMatrix-method}
\alias{is_sparse,lgCMatrix-method}
\alias{is_sparse,dgRMatrix-method}
\alias{is_sparse,lgRMatrix-method}

\title{Is an array-like object sparse?}

\description{
  \code{is_sparse} indicates whether an array-like object is
  \emph{structurally} sparse or not.

  Note that this is about \emph{data representation} and not about the
  data itself. For example, \code{is_sparse()} always returns \code{FALSE}
  on an ordinary matrix, even if the matrix contains 99\% zeros, because
  it uses a dense representation of the data.
  OTOH \code{is_sparse()} always returns \code{TRUE} on a
  \linkS4class{dgCMatrix} object from the \pkg{Matrix} package,
  even if the data contains no zeros, because the object uses a sparse
  representation of the data.
}

\usage{
is_sparse(x)
}

\arguments{
  \item{x}{
    An array-like object.
    This can be an ordinary array, a \link[SparseArray]{SparseArray}
    object from the \pkg{SparseArray} package, a \linkS4class{dgCMatrix}
    object from the \pkg{Matrix} package, a \link[DelayedArray]{DelayedArray}
    object from the \pkg{DelayedArray} package, or any object with an array
    semantic (i.e. an object for which \code{dim(x)} is not NULL).
  }
}

\value{
  \code{TRUE} or \code{FALSE}
}

\seealso{
  \itemize{
    \item \link[base]{array} objects in base R.

    \item \linkS4class{dgCMatrix} objects in the \pkg{Matrix} package.
  }
}

\examples{
m <- matrix(0L, nrow=50, ncol=20)
stopifnot(identical(is_sparse(m), FALSE))

dgc <- as(m + runif(1000), "dgCMatrix")
stopifnot(identical(is_sparse(dgc), TRUE))
}
\keyword{methods}
