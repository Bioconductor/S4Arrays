\name{is_sparse}

\alias{is_sparse}
\alias{is_sparse<-}
\alias{is_sparse,ANY-method}
\alias{is_sparse,dgCMatrix-method}
\alias{is_sparse,lgCMatrix-method}
\alias{is_sparse,dgRMatrix-method}
\alias{is_sparse,lgRMatrix-method}

\title{Check for sparse representation}

\description{
  \code{is_sparse} indicates whether an object (typically array-like)
  uses a sparse representation of the data or not.

  Note that this is about \emph{data representation} and not about the
  data itself. For example, \code{is_sparse()} always returns \code{FALSE}
  on an \emph{ordinary} matrix, even if the matrix contains 99\% zeros,
  because the data in such a matrix is always stored in a dense form.
  OTOH \code{is_sparse()} always returns \code{TRUE} on a
  \link[SparseArray]{SparseArray} derivative from the \pkg{SparseArray}
  package, or on a \linkS4class{dgCMatrix} object from the \pkg{Matrix}
  package, even if the data contains no zeros, because these objects
  use a sparse representation of the data.
}

\usage{
is_sparse(x)
}

\arguments{
  \item{x}{
    Any object, but will typically be an array-like object.

    Examples of array-like objects: ordinary arrays,
    \link[SparseArray]{SparseArray} objects from the \pkg{SparseArray} package,
    \linkS4class{dgCMatrix} objects from the \pkg{Matrix} package,
    \link[DelayedArray]{DelayedArray} objects from the \pkg{DelayedArray}
    package, or any object with an array semantic (i.e. an object for
    which \code{dim(x)} is not NULL).
  }
}

\value{
  \code{TRUE} or \code{FALSE}
}

\seealso{
  \itemize{
    \item \code{\link{read_block}} to read a block of data from an
          array-like object.

    \item \link[base]{array} and \link[base]{matrix} objects in base R.

    \item \linkS4class{dgCMatrix} objects implemented in the \pkg{Matrix}
          package.
  }
}

\examples{
m <- matrix(0L, nrow=50, ncol=20)
stopifnot(identical(is_sparse(m), FALSE))

dgc <- as(m + runif(1000), "dgCMatrix")
stopifnot(identical(is_sparse(dgc), TRUE))
}
\keyword{methods}
