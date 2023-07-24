\name{bind-arrays}

\alias{bind-arrays}
\alias{bind arrays}

\alias{abind}
\alias{abind,ANY-method}

\alias{arbind}
\alias{acbind}
\alias{arbind,array-method}
\alias{acbind,array-method}


\title{Combine multidimensional array-like objects}

\description{
  Combine multidimensional array-like objects.

  NOTE: This man page is for the \code{abind}
  \emph{S4 generic function} defined in the \pkg{S4Arrays} package.
  See \code{?abind::\link[abind]{abind}} for the default method
  (defined in the \pkg{abind} package).
  Bioconductor packages can define specific methods for objects
  not supported by the default method.
}

\usage{
## Bind array-like objects along an arbitrary dimension:
abind(..., along=NULL, rev.along=NULL)

## Bind array-like objects along their first or second dimension:
arbind(...)
acbind(...)
}

\arguments{
  \item{...}{
    The array-like objects to bind.
  }
  \item{along, rev.along}{
    See \code{?abind::\link[abind]{abind}} for a description of
    these arguments.
  }
}

\value{
  An array-like object, typically of the same class as the input objects if
  they all have the same class.
}

\seealso{
  \itemize{
    \item \code{abind::\link[abind]{abind}} in the \pkg{abind} package
          for the default \code{abind} method.

    \item \code{\link[base]{rbind}} and \code{\link[base]{cbind}} in the
          \pkg{base} package for the corresponding operations on matrix-like
          objects.
  }
}

\examples{
a1 <- array(1:60, c(3, 5, 4),
            dimnames=list(NULL, paste0("A1y", 1:5), NULL))
a2 <- array(101:240, c(7, 5, 4),
            dimnames=list(paste0("A2x", 1:7), paste0("A2y", 1:5), NULL))
a3 <- array(10001:10100, c(5, 5, 4),
            dimnames=list(paste0("A3x", 1:5), NULL, paste0("A3z", 1:4)))

abind(a1, a2, a3, along=1)  # same as 'arbind(a1, a2, a3)'

m2 <- matrix(1:35, nrow=7)
abind(m2, a2, along=3)
abind(a2, m2, along=3)

abind(m2, m2+0.5, rev.along=0)  # same as 'abind(m2, m2+0.5, along=3)'
}
\keyword{array}
\keyword{manip}
