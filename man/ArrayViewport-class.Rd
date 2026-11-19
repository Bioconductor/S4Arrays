\name{ArrayViewport-class}
\docType{class}

\alias{class:ArrayViewport}
\alias{ArrayViewport-class}
\alias{ArrayViewport}

\alias{class:DummyArrayViewport}
\alias{DummyArrayViewport-class}
\alias{DummyArrayViewport}

\alias{class:SafeArrayViewport}
\alias{SafeArrayViewport-class}
\alias{SafeArrayViewport}

\alias{refdim}
\alias{refdim,ArrayViewport-method}
\alias{ranges,DummyArrayViewport-method}
\alias{ranges,SafeArrayViewport-method}
\alias{start,ArrayViewport-method}
\alias{width,ArrayViewport-method}
\alias{end,ArrayViewport-method}
\alias{dim,ArrayViewport-method}

\alias{classNameForDisplay,ArrayViewport-method}
\alias{show,ArrayViewport-method}

\alias{makeNindexFromArrayViewport}

\title{ArrayViewport objects}

\description{
  An ArrayViewport object represents a block-shaped view on a multidimensional
  array, that is, a view made of consecutive indices along each dimension
  of the array.

  In this context, the array on which the view is defined is referred to as
  the reference array.

  ArrayViewport and \link{ArrayGrid} objects are primarily used in the
  context of block processing of array-like objects.
}

\usage{
## Constructor function:
ArrayViewport(refdim, ranges=NULL)
}

\arguments{
  \item{refdim}{
    An integer vector containing the dimensions of the reference array.
  }
  \item{ranges}{
    An \link[IRanges]{IRanges} object with one range per dimension in
    the reference array. If omitted, an ArrayViewport object that represents
    a view on the full reference array is returned.
  }
}

\details{
  Note that an ArrayViewport object only stores the geometry of the view,
  which can be described by a collection of integer intervals, one along
  each dimension in the reference array.

  In particular, an ArrayViewport object doesn't store any array data, or
  any details about the class or data type of the reference array.
}

\value{
  An ArrayViewport object.
}

\seealso{
  \itemize{
    \item \code{\link{read_block}} to read a block of data from an
          array-like object.

    \item \link{ArrayGrid} objects.

    \item \code{\link[DelayedArray]{blockApply}} and family, in the
          \pkg{DelayedArray} package, for convenient block processing
          of an array-like object.

    \item \link[base]{array} and \link[base]{matrix} objects in base R.
  }
}

\examples{
## Create ArrayViewport objects on a reference array of dimensions
## 15 x 10 x 8:

ranges <- IRanges(c("11-15", "1-10", "3-4"))
viewport1 <- ArrayViewport(c(15L, 10L, 8L), ranges)
viewport1

refdim(viewport1)
dim(viewport1)     # the dimensions of the subarray in the view
length(viewport1)  # the number of array elements in the view
ranges(viewport1)

ranges <- IRanges(start=c(11, 8, 6), width=c(5, 0, 2))  # empty view!
viewport2 <- ArrayViewport(c(15L, 10L, 8L), ranges)
viewport2

refdim(viewport2)
dim(viewport2)
length(viewport2)  # no array elements in the view!
ranges(viewport2)

## Create ArrayViewport objects on a reference matrix of dimensions
## 150 x 20:

ranges <- IRanges(end=c(150, 20), width=10)  # bottom right 10x10 tile
viewport3 <- ArrayViewport(c(150L, 20L), ranges)
viewport3

refdim(viewport3)
dim(viewport3)
length(viewport3)
ranges(viewport3)

viewport4 <- ArrayViewport(c(150L, 20L))  # view on the full matrix
viewport4

refdim(viewport4)
dim(viewport4)     # same as 'refdim(viewport4)'
length(viewport4)  # same as length of full matrix
ranges(viewport4)
}
\keyword{classes}
\keyword{methods}
