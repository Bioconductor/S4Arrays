\name{SparseArray}
\docType{class}

\alias{class:SparseArray}
\alias{SparseArray-class}
\alias{SparseArray}

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

  Users should be able to operate on SparseArray derivatives via the
  \emph{standard array API} defined in base R, that is, using things
  like \code{dim()}, \code{length()}, \code{dimnames()}, \code{[},
  \code{drop()}, \code{t()}, \code{rbind()}, \code{cbind()}, etc...

  SparseArray derivatives also support \code{type()}, \code{`type<-`()},
  and \code{is_sparse()}, three generic functions defined in this
  package (the \pkg{S4Arrays} package).
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
