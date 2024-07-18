\name{rowsum}

\alias{rowsum}
\alias{colsum}
\alias{colsum,ANY-method}
\alias{colsum,matrix-method}

\title{Compute column/row sums of a matrix-like object, for groups
       of rows/columns}

\description{
  \code{rowsum()} computes column sums across rows of a numeric matrix-like
  object for each level of a grouping variable.

  \code{colsum()} computes row sums across columns of a numeric matrix-like
  object for each level of a grouping variable.

  NOTE: This man page is for the \code{rowsum} and \code{colsum}
  \emph{S4 generic functions} defined in the \pkg{S4Arrays} package.
  See \code{?base::\link[base]{rowsum}} for the default \code{rowsum()}
  method (defined in the \pkg{base} package).
  Bioconductor packages can define specific methods for objects
  (typically matrix-like) not supported by the default methods.
}

\usage{
rowsum(x, group, reorder=TRUE, ...)

colsum(x, group, reorder=TRUE, ...)
}

\arguments{
  \item{x}{
    A numeric matrix-like object.
  }
  \item{group, reorder, ...}{
    See \code{?base::\link[base]{rowsum}} for a description of
    these arguments.
  }
}

\value{
  See \code{?base::\link[base]{rowsum}} for the value returned
  by the default method.

  The default \code{colsum()} method simply does
  \code{t(rowsum(t(x), group, reorder=reorder, ...))}.

  Specific methods defined in Bioconductor packages should
  behave as consistently as possible with the default methods.
}

\seealso{
  \itemize{
    \item \code{base::\link[base]{rowsum}} for the default
          \code{rowsum} method.

    \item \code{\link[methods]{showMethods}} for displaying a summary of the
          methods defined for a given generic function.

    \item \code{\link[methods]{selectMethod}} for getting the definition of
          a specific method.

    \item \link[DelayedArray]{rowsum,DelayedMatrix-method} in the
          \pkg{DelayedArray} package for an example of a specific
          \code{rowsum} method (defined for \link[DelayedArray]{DelayedMatrix}
          objects).
  }
}

\examples{
rowsum  # note the dispatch on the 'x' arg only
showMethods("rowsum")
selectMethod("rowsum", "ANY")     # the default rowsum() method

colsum  # note the dispatch on the 'x' arg only
showMethods("colsum")
selectMethod("colsum", "ANY")     # the default colsum() method
selectMethod("colsum", "matrix")  # colsum() method for ordinary matrices
}

\keyword{methods}
\keyword{manip}
