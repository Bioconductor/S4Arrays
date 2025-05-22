\name{array_recycling}

\alias{array_recycling}

\alias{class:tile}
\alias{tile-class}
\alias{tile}

\alias{as_tile}
\alias{Ops,array,tile-method}
\alias{Ops,tile,array-method}
\alias{Ops,tile,tile-method}

\title{Multidimensional array recycling}

\description{
  In base R, arithmetic and other binary operations between a matrix or
  array and a vector don't let the user control how the latter should
  be recycled: it's always recycled along the first dimension of the
  matrix or array. This has led users to use various tricks when they
  need recycling along the second dimension (a.k.a. "horizontal recycling"),
  like the popular "double-transposition" trick: \code{t(t(m) / colSums(m))}.
  However this is not only inelegant, it's also inefficient.

  \code{as_tile()} is meant to address that.

  It also allows arithmetic and other binary operations between two arrays
  of distinct dimensions (rejected as "non-conformable" by arithmetic
  operations in base R), typically between a small one (the tile) and
  a bigger one, as long as their geometries are compatible.
}

\usage{
as_tile(x, along=1L, dim=NULL)
}

\arguments{
  \item{x}{
    An array-like object or a vector.
  }
  \item{along}{
    Can only be used when \code{x} is a vector.
    Must be a single positive integer indicating the "orientation" of the
    tile to be created, that is, the dimension along which \code{x} will be
    recycled by arithmetic operations and other binary operations.
  }
  \item{dim}{
    \code{NULL} or the dimensions (supplied as an integer vector) of
    the tile to be created.
  }
}

\value{
  A tile object (tile is an extension of array).
}

\seealso{
  \itemize{
    \item \code{base::\link{colSums}} in the \pkg{base} package.

    \item \link[base]{array} and \link[base]{matrix} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## 2D EXAMPLES
## ---------------------------------------------------------------------
m0 <- matrix(1:54, nrow=6)
x <- c(-1, 0, 100)

## Arithmetic operations in base R recycle 'x' along the first dimension
## of the matrix ("vertical recycling"):
m0 * x

## To recycle 'x' along the second dimension of the matrix ("horizontal
## recycling"), we turn it into an "horizontal" tile:
t <- as_tile(x, along=2)
t

m0 * t

## The above produces the same result as the double-transposition trick
## but is more efficient (and also makes code easier to read):
stopifnot(identical(m0 * t, t(t(m0) * x)))

## A less artificial example:
cs0 <- colSums(m0)
m <- m0 / as_tile(cs0, along=2)

stopifnot(all.equal(colSums(m), rep(1, ncol(m))))  # sanity check

## Using an arbitrary 2D tile:

t <- m0[1:2, 1:3]
t
## Unfortunately arithmetic operations in base R refuse to operate on
## arrays that don't have the same dimensions:
\dontrun{
  m0 / t  # ERROR! (non-conformable arrays)
}

## Wrapping 't' in a tile object makes this work:
m0 / as_tile(t)

## ---------------------------------------------------------------------
## 3D EXAMPLES
## ---------------------------------------------------------------------
## Note that colSums() supports multidimensional arrays. In this case
## the user can use the 'dims' argument to control how the array should
## be sliced into "columns". See '?base::colSums' for the details.

a <- array(runif(300), dim=c(10, 5, 6))

## Using 'dims=2' indicates that the columns are the 2D slices obtained
## by slicing the array along its 3rd dimension. With this slicing, each
## column is a 10 x 5 matrix.
cs2 <- colSums(a, dims=2)
cs2  # vector of length 6 (one value per 2D slice)

t <- as_tile(cs2, along=3)
a / t

stopifnot(all.equal(colSums(a / t, dims=2), rep(1, dim(a)[3])))  # sanity check

## By default (i.e. when 'dims=1') the array is considered to be made
## of 5*6 columns of length 10.
cs1 <- colSums(a)
cs1  # 5 x 6 matrix

t <- as_tile(cs1, dim=c(1L, dim(cs1)))
a / t

stopifnot(all.equal(colSums(a / t), matrix(1, dim(a)[2], dim(a)[3])))  # sanity check
}
\keyword{array}
\keyword{methods}
