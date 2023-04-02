\name{write_block}

\alias{write_block}
\alias{write_block,ANY-method}

\title{Write array blocks}

\description{
  Use \code{write_block} to write a block of data to an array-like object.

  Note that this function is typically used in the context of block processing
  of on-disk objects (e.g. \link[DelayedArray]{DelayedArray} objects), often
  in combination with \code{\link{read_block}}.
}

\usage{
write_block(sink, viewport, block)
}

\arguments{
  \item{sink}{
    A **writable** array-like object. This is typically a
    \link[DelayedArray]{RealizationSink} derivative (RealizationSink
    is a virtual class defined in the \pkg{DelayedArray} package),
    but not necessarily. See \code{?\link[DelayedArray]{RealizationSink}}
    in the \pkg{DelayedArray} package for more information about
    RealizationSink objects.

    Although \code{write_block()} will typically be used on a RealizationSink
    derivative, it can also be used on an ordinary array or other in-memory
    array-like object that supports subassignment (\code{`[<-`}), like
    a \link[SparseArray]{SparseArray} object from the \pkg{SparseArray}
    package, or a dgCMatrix object from the \pkg{Matrix} package.
  }
  \item{viewport}{
    An \link{ArrayViewport} object compatible with \code{sink}, that is,
    such that \code{refdim(viewport)} is identical to \code{dim(sink)}.
  }
  \item{block}{
    An array-like object of the same dimensions as \code{viewport}.
  }
}

\value{
  The modified array-like object \code{sink}.
}

\seealso{
  \itemize{
    \item \link{ArrayGrid} for ArrayGrid and ArrayViewport objects.

    \item \link[SparseArray]{SparseArray} objects implemented in the
          \pkg{SparseArray} package.

    \item \code{\link{read_block}} to read a block of data from an
          array-like object.

    \item \code{\link[DelayedArray]{blockApply}} and family, in the
          \pkg{DelayedArray} package, for convenient block processing
          of an array-like object.

    \item \link[DelayedArray]{RealizationSink} objects implemented in the
          \pkg{DelayedArray} package for more realistic \code{write_block}
          examples.

    \item \link[base]{array} and \link[base]{matrix} objects in base R.
  }
}

\examples{
## Please note that, although educative, the examples below are somewhat
## artificial and do not illustrate real-world usage of write_block().
## See '?RealizationSink' in the DelayedArray package for more realistic
## read_block/write_block examples.

## ---------------------------------------------------------------------
## BASIC EXAMPLE 1: WRITE A BLOCK TO AN ORDINARY MATRIX (DENSE)
## ---------------------------------------------------------------------
m1 <- matrix(1:30, ncol=5)
m1

## Define the viewport on 'm1' to write the data to:
block1_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m1), IRanges(c(3, 2), width=block1_dim))
viewport1

## Data to write:
block1 <- read_block(m1, viewport1) + 1000L

## Write the block:
m1A <- write_block(m1, viewport1, block1)
m1A

## Sanity checks:
stopifnot(identical(`[<-`(m1, 3:6, 2:4, value=block1), m1A))
m1B <- write_block(m1, viewport1, as(block1, "dgCMatrix"))
stopifnot(identical(m1A, m1B))

## ---------------------------------------------------------------------
## BASIC EXAMPLE 2: WRITE A BLOCK TO A SPARSE MATRIX
## ---------------------------------------------------------------------
m2 <- rsparsematrix(12, 20, density=0.2,
                    rand.x=function(n) sample(25, n, replace=TRUE))
m2

## Define the viewport on 'm2' to write the data to:
block2_dim <- c(2, 20)
viewport2 <- ArrayViewport(dim(m2), IRanges(c(1, 1), width=block2_dim))
viewport2

## Data to write:
block2 <- matrix(1001:1040, nrow=2)

## Write the block:
m2A <- write_block(m2, viewport2, block2)
m2A

## Sanity checks:
stopifnot(identical(`[<-`(m2, 1:2, , value=block2), m2A))
m2B <- write_block(m2, viewport2, as(block2, "dgCMatrix"))
stopifnot(identical(m2A, m2B))

## ---------------------------------------------------------------------
## BASIC EXAMPLE 3: WRITE A BLOCK TO A 3D ARRAY
## ---------------------------------------------------------------------
a3 <- array(1:60, dim=5:3)

## Define the viewport on 'a3' to write the data to:
block3_dim <- c(2, 4, 1)
viewport3 <- ArrayViewport(dim(a3), IRanges(c(1, 1, 3), width=block3_dim))
viewport3

## Data to write:
block3 <- array(-(1:8), dim=block3_dim)

## Write the block:
a3A <- write_block(a3, viewport3, block3)
a3A

## Sanity checks:
stopifnot(identical(`[<-`(a3, 1:2, , 3, value=block3), a3A))
a3B <- write_block(a3, viewport3, as(block3, "SparseArray"))
stopifnot(identical(a3A, a3B))

## ---------------------------------------------------------------------
## BASIC EXAMPLE 4: WRITE BLOCKS DEFINED BY A GRID
## ---------------------------------------------------------------------
a4 <- array(NA_real_, dim=6:4)

## Define a grid of 2x3x2 blocks on 'a4':
grid4 <- RegularArrayGrid(dim(a4), spacings=c(2,3,2))
grid4
nblock <- length(grid4)  # number of blocks

## Walk on the grid and write blocks of random data:
for (bid in seq_len(nblock)) {
    viewport <- grid4[[bid]]
    block <- array(runif(length(viewport)), dim=dim(viewport))
    cat("====== Write block ", bid, "/", nblock, " ======\n", sep="")
    a4 <- write_block(a4, viewport, block)
}
a4

## ---------------------------------------------------------------------
## BASIC EXAMPLE 5: READ, PROCESS, AND WRITE BLOCKS DEFINED BY TWO GRIDS
## ---------------------------------------------------------------------
## Say we have a 3D array and want to collapse its 3rd dimension by
## summing the array elements that are stacked vertically, that is, we
## want to compute the matrix 'm' obtained by doing 'sum(a[i, j, ])' for
## all valid i and j. There are several ways to do this.

## 1. Here is a solution based on apply():

collapse_3rd_dim <- function(a) apply(a, MARGIN=1:2, sum)

## 2. Here is a slightly more efficient solution:

collapse_3rd_dim <- function(a) {
    m <- matrix(0, nrow=nrow(a), ncol=ncol(a))
    for (z in seq_len(dim(a)[[3]]))
        m <- m + a[ , , z]
    m
}

## 3. And here is a block-processing solution that involves two grids,
##    one for the sink, and one for the input:

a5 <- array(runif(8000), dim=c(25, 40, 8))  # input
m  <- array(NA_real_, dim=dim(a5)[1:2])     # sink

## Since we're going to walk on the two grids simultaneously, read a
## block from 'a5' and write it to 'm', we need to make sure that we
## define grids that are "aligned". More precisely, the two grids must
## have the same number of viewports, and the viewports in one must
## correspond to the viewports in the other one:
m_grid  <- RegularArrayGrid(dim(m),  spacings=c(10, 10))
a5_grid <- RegularArrayGrid(dim(a5), spacings=c(10, 10, dim(a5)[[3]]))

## Let's check that our two grids are actually "aligned":
stopifnot(identical(length(m_grid), length(a5_grid)))
stopifnot(identical(dims(m_grid), dims(a5_grid)[ , 1:2, drop=FALSE]))

## Walk on the two grids simultaneously, and read/collapse/write blocks:
for (bid in seq_along(m_grid)) {
    ## Read block from 'a5'.
    a5_viewport <- a5_grid[[bid]]
    block <- read_block(a5, a5_viewport)
    ## Collapse it.
    block <- collapse_3rd_dim(block)
    ## Write the collapsed block to 'm'.
    m_viewport <- m_grid[[bid]]
    m <- write_block(m, m_viewport, block)
}

## Sanity checks:
stopifnot(identical(dim(a5)[1:2], dim(m)))
stopifnot(identical(sum(a5), sum(m)))
stopifnot(identical(collapse_3rd_dim(a5), m))

## See '?RealizationSink' in the DelayedArray package for a more
## realistic "array collapse" example where the blocks are written
## to a RealizationSink object.
}
\keyword{methods}
