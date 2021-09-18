\name{SparseArray}
\docType{class}

\alias{class:SparseArray}
\alias{SparseArray-class}
\alias{SparseArray}

\alias{dimnames,SparseArray-method}
\alias{dimnames<-,SparseArray-method}
\alias{is_sparse}
\alias{is_sparse<-}
\alias{is_sparse,ANY-method}
\alias{is_sparse,SparseArray-method}
\alias{sparsity}

\alias{show,SparseArray-method}

\title{SparseArray objects}

\description{
  SparseArray is a virtual class to be extended by classes that aim at
  representing in-memory multidimensional sparse arrays.

  It has currently two concrete subclasses: \link{COO_SparseArray}
  and \link{SVT_SparseArray}.

  Here is where SparseArray and its two concrete subclasses fit in the
  \link{Array} class hierarchy:
  \preformatted{
                          Array
                            ^
                            |
                       SparseArray
                       ^         ^
                       |         |
            COO_SparseArray   SVT_SparseArray}
}

\details{
  Any object that belongs to a class that extends SparseArray is called
  a \emph{SparseArray derivative}.

  Most of the \emph{standard array API} defined in base R should work
  on SparseArray derivatives, including \code{dim()}, \code{length()},
  \code{dimnames()}, \code{`dimnames<-`()}, \code{[}, \code{drop()},
  \code{t()}, \code{rbind()}, \code{cbind()}, etc...

  SparseArray derivatives also support \code{type()}, \code{`type<-`()},
  \code{is_sparse()}, \code{sparsity()}, \code{arbind()} and \code{acbind()},
  six generic functions defined in this package (the \pkg{S4Arrays} package).

  \code{sparsity(x)} returns the ratio of number of zero-valued elements
  over total number of elements (\code{length(x)}) in array-like object
  \code{x}.
}

\seealso{
  \itemize{
    \item The \link{COO_SparseArray} and \link{SVT_SparseArray} classes.

    \item S4 classes \linkS4class{dgCMatrix} and \linkS4class{lgCMatrix}
          defined in the \pkg{Matrix} package, for the de facto standard
          of sparse matrix representations in R.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
showClass("SparseArray")
}

\keyword{methods}
\keyword{classes}
