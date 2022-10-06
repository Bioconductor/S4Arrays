\name{read_block}

\alias{read_block}

\alias{read_block_as_dense}
\alias{read_block_as_dense,ANY-method}

\title{Read array blocks}

\description{
  Use \code{read_block} to read a block of data from an array-like object.

  Note that this function is typically used in the context of block processing
  of on-disk objects (e.g. \link[DelayedArray]{DelayedArray} objects), often
  in combination with \code{\link{write_block}}.
}

\usage{
read_block(x, viewport, as.sparse=NA)
}

\arguments{
  \item{x}{
    An array-like object. This can be an ordinary array, a
    \link[SparseArray]{SparseArray} object from the \pkg{SparseArray} package,
    a \linkS4class{dgCMatrix} object from the \pkg{Matrix} package,
    a \link[DelayedArray]{DelayedArray} object from the \pkg{DelayedArray}
    package, or any object with an array semantic (i.e. an object for
    which \code{dim(x)} is not NULL).
  }
  \item{viewport}{
    An \link{ArrayViewport} object compatible with \code{x}, that is,
    such that \code{refdim(viewport)} is identical to \code{dim(x)}.
  }
  \item{as.sparse}{
    Can be \code{FALSE}, \code{TRUE}, or \code{NA}.

    If \code{FALSE}, the block is returned as an ordinary
    array (a.k.a. dense array).

    If \code{TRUE}, it's returned as a \link[SparseArray]{SparseArray} object.

    If \code{NA} (the default), the block is returned as an
    ordinary array if \code{is_sparse(x)} is \code{FALSE} and
    as a \link[SparseArray]{SparseArray} object otherwise.
    In other words, using \code{as.sparse=NA} is equivalent
    to using \code{as.sparse=is_sparse(x)}.
    This preserves sparsity and is the most efficient way to read a block.

    Note that when returned as a 2D \link[SparseArray]{SparseArray} object
    with numeric or logical data, a block can easily and efficiently
    be coerced to a \link[Matrix]{sparseMatrix} derivative from the
    \pkg{Matrix} package with \code{as(block, "sparseMatrix")}.
    This will return a dgCMatrix object if \code{type(block)}
    is \code{"double"} or \code{"integer"}, and a lgCMatrix
    object if it's \code{"logical"}.
  }
}

\value{
  A block of data. More precisely, the data from \code{x} that belongs to the
  block delimited by the specified viewport.

  The block of data is returned either as an ordinary (dense) array
  or as a \link[SparseArray]{SparseArray} object from the \pkg{SparseArray}
  package.

  Note that the returned block of data is guaranteed to have the same type
  as \code{x} and the same dimensions as the \code{viewport}.
  More formally, if \code{block} is the value returned by
  \code{read_block(x, viewport)}, then:
  \preformatted{    identical(type(block), type(x))}
  and
  \preformatted{    identical(dim(block), dim(viewport))}
  are always \code{TRUE}.
}

\seealso{
  \itemize{
    \item \link{ArrayGrid} for ArrayGrid and ArrayViewport objects.

    \item \link[SparseArray]{SparseArray} objects in the \pkg{SparseArray}
          package.

    \item \code{\link{write_block}}.

    \item \code{\link[DelayedArray]{blockApply}} and family, in the
          \pkg{DelayedArray} package, for convenient block processing
          of an array-like object.

    \item \link[Matrix]{dgCMatrix-class} and \link[Matrix]{lgCMatrix-class}
          objects in the \pkg{Matrix} package.

    \item \link[DelayedArray]{DelayedArray} objects in the \pkg{DelayedArray}
          package.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## Please note that, although educative, the examples below are somewhat
## artificial and do not illustrate real-world usage of read_block().
## See '?RealizationSink' in the DelayedArray package for more realistic
## read_block/write_block examples.

## ---------------------------------------------------------------------
## BASIC EXAMPLE 1: READ A BLOCK FROM AN ORDINARY MATRIX (DENSE)
## ---------------------------------------------------------------------
m1 <- matrix(1:30, ncol=5)
m1

## Define the viewport on 'm1' to read the data from:
block1_dim <- c(4, 3)
viewport1 <- ArrayViewport(dim(m1), IRanges(c(3, 2), width=block1_dim))
viewport1

## Read the block:
block1 <- read_block(m1, viewport1)  # same as m1[3:6, 2:4, drop=FALSE]
block1

## Use 'as.sparse=TRUE' to read the block as sparse object:
block1b <- read_block(m1, viewport1, as.sparse=TRUE)
block1b
is_sparse(block1b)  # TRUE
class(block1b)  # a SparseArray object

## Sanity checks:
stopifnot(identical(type(m1), type(block1)))
stopifnot(identical(dim(viewport1), dim(block1)))
stopifnot(identical(m1[3:6, 2:4, drop=FALSE], block1))
stopifnot(is(block1b, "SparseArray"))
stopifnot(identical(type(m1), type(block1b)))
stopifnot(identical(dim(viewport1), dim(block1b)))
stopifnot(identical(block1, as.array(block1b)))

## ---------------------------------------------------------------------
## BASIC EXAMPLE 2: READ A BLOCK FROM A SPARSE MATRIX
## ---------------------------------------------------------------------
m2 <- rsparsematrix(12, 20, density=0.2,
                    rand.x=function(n) sample(25, n, replace=TRUE))
m2
is_sparse(m2)  # TRUE

## Define the viewport on 'm2' to read the data from:
block2_dim <- c(2, 20)
viewport2 <- ArrayViewport(dim(m2), IRanges(c(1, 1), width=block2_dim))
viewport2

## By default, read_block() preserves sparsity:
block2 <- read_block(m2, viewport2)
block2
is_sparse(block2)  # TRUE
class(block2)  # a SparseArray object

## Use 'as.sparse=FALSE' to force read_block() to return an ordinary
## matrix or array:
block2b <- read_block(m2, viewport2, as.sparse=FALSE)
block2b
as(block2b, "sparseMatrix")

## Sanity checks:
stopifnot(is(block2, "SparseArray"))
stopifnot(identical(type(m2), type(block2)))
stopifnot(identical(dim(viewport2), dim(block2)))
stopifnot(identical(type(m2), type(block2b)))
stopifnot(identical(dim(viewport2), dim(block2b)))
stopifnot(identical(block2b, as.array(block2)))

## ---------------------------------------------------------------------
## BASIC EXAMPLE 3: READ A BLOCK FROM A 3D ARRAY
## ---------------------------------------------------------------------
a3 <- array(1:60, dim=5:3)

## Define the viewport on 'a3' to read the data from:
block3_dim <- c(2, 4, 1)
viewport3 <- ArrayViewport(dim(a3), IRanges(c(1, 1, 3), width=block3_dim))
viewport3

## Read the block:
block3 <- read_block(a3, viewport3)  # same as a3[1:2, 1:4, 3, drop=FALSE]
block3

## Note that unlike [, read_block() never drops dimensions.

## Sanity checks:
stopifnot(identical(type(a3), type(block3)))
stopifnot(identical(dim(viewport3), dim(block3)))
stopifnot(identical(a3[1:2, 1:4, 3, drop=FALSE], block3))

## ---------------------------------------------------------------------
## BASIC EXAMPLE 4: READ AND PROCESS BLOCKS DEFINED BY A GRID
## ---------------------------------------------------------------------
a4 <- array(runif(120), dim=6:4)

## Define a grid of 2x3x2 blocks on 'a4':
grid4 <- RegularArrayGrid(dim(a4), spacings=c(2,3,2))
grid4
nblock <- length(grid4)  # number of blocks

## Walk on the grid and print the corresponding blocks:
for (bid in seq_len(nblock)) {
    viewport <- grid4[[bid]]
    block <- read_block(a4, viewport)
    cat("====== Block ", bid, "/", nblock, " ======\n", sep="")
    print(block)
}

## Walk on the grid and compute the sum of each block:
block_sums <- sapply(grid4,
    function(viewport) sum(read_block(a4, viewport))
)
block_sums

## Sanity checks:
stopifnot(identical(length(block_sums), nblock))
stopifnot(all.equal(sum(block_sums), sum(a4)))

## ---------------------------------------------------------------------
## THE read_block/write_block COMBO
## ---------------------------------------------------------------------
## See '?write_block' for examples that use the read_block/write_block
## combo.
}
\keyword{methods}
