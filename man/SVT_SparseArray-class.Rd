\name{SVT_SparseArray-class}
\docType{class}

\alias{class:NULL_OR_list}
\alias{NULL_OR_list-class}
\alias{NULL_OR_list}

\alias{class:SVT_SparseArray}
\alias{SVT_SparseArray-class}
\alias{SVT_SparseArray}

\alias{type,SVT_SparseArray}
\alias{type<-,SVT_SparseArray}

\alias{as.array.SVT_SparseArray}
\alias{as.array,SVT_SparseArray-method}
\alias{coerce,array,SVT_SparseArray-method}

\alias{coerce,SVT_SparseArray,dgCMatrix-method}
\alias{coerce,SVT_SparseArray,lgCMatrix-method}
\alias{coerce,CsparseMatrix,SVT_SparseArray-method}

\alias{coerce,SVT_SparseArray,COO_SparseArray-method}
\alias{coerce,COO_SparseArray,SVT_SparseArray-method}

\alias{coerce,ANY,SVT_SparseArray-method}
\alias{coerce,RsparseMatrix,SVT_SparseArray-method}

\alias{t.SVT_SparseArray}
\alias{t,SVT_SparseArray-method}

\title{SVT_SparseArray objects}

\description{
  SVT_SparseArray is a new class for efficient in-memory representation
  of multidimensional sparse arrays.

  It replaces the older and less efficient \link{COO_SparseArray} class.
}

\usage{
## Constructor function:
SVT_SparseArray(x, type=NA)
}

\arguments{
  \item{x}{
    An ordinary array, or a dgCMatrix/lgCMatrix object, or any array-like
    object that supports coercion to SVT_SparseArray.
  }
  \item{type}{
    A single string specifying the requested type of the object.

    Normally, the SVT_SparseArray object returned by the constructor
    function has the same \code{type()} as \code{x} but the user can use
    the \code{type} argument to request a different type. Note that doing:
    \preformatted{    svt <- SVT_SparseArray(x, type=type)}
    is equivalent to doing:
    \preformatted{    svt <- SVT_SparseArray(x)
    type(svt) <- type}
    but the former is more convenient and will be generally more efficient.

    Supported types are all R atomic types plus \code{"list"}.
  }
}

\details{
  SVT_SparseArray is a concrete subclass of the \link{SparseArray}
  virtual class. This makes SVT_SparseArray objects SparseArray derivatives.

  The nonzero data in a SVT_SparseArray object is stored in a \emph{Sparse
  Vector Tree}. We'll refer to this internal data representation as
  the \emph{SVT layout}. See the "SVT layout" section below for more
  information.

  The SVT layout is similar to the CSC layout (compressed, sparse,
  column-oriented format) used by dgCMatrix/lgCMatrix objects from
  the \pkg{Matrix} package, but with the following improvements:
  \itemize{
      \item The SVT layout supports sparse arrays of arbitrary dimensions.
      \item With the SVT layout, the sparse data can be of any type.
            With dgCMatrix/lgCMatrix objects, it must be either
            \code{"double"} or \code{"logical"}.
      \item The SVT layout imposes no limit on the number of nonzero elements
            that can be stored. With dgCMatrix/lgCMatrix objects, this number
            must be < 2^31.
      \item Overall, the SVT layout allows more efficient operations on
            SVT_SparseArray objects.
  }
}

\value{
  An SVT_SparseArray object.
}

\section{SVT layout}{
  An SVT (Sparse Vector Tree) is a tree of depth N - 1 where N is the number
  of dimensions of the sparse array.

  The leaves in the tree can only be of two kinds: NULL or \emph{leaf vector}.
  Leaves that are leaf vectors can only be found at the deepest level in the
  tree (i.e. at depth N - 1). All leaves found at a lower depth must be NULLs.

  A leaf vector represents a sparse vector of length equal to the first
  dimension of the sparse array. This is done using a set of offset/value
  pairs sorted by strictly ascending offset.
  More precisely, a leaf vector is represented by an ordinary list of 2
  parallel vectors:
  \enumerate{
    \item an integer vector of offsets (i.e. 0-based positions);
    \item a vector (atomic or list) of nonzero values.
  }
  The 2nd vector determines the type of the leaf vector i.e. \code{"double"},
  \code{"integer"}, \code{"logical"}, etc...
  All the leaf vectors in the SVT have the type of the sparse array.

  Examples:
  \itemize{
    \item An SVT_SparseArray object with 1 dimension has its nonzero data
          stored in an SVT of depth 0. Such SVT is represented by a
          single "leaf vector".

    \item An SVT_SparseArray object with 2 dimensions has its nonzero data
          stored in an SVT of depth 1. Such SVT is represented by a list of
          length the extend of the 2nd dimension (nb of columns). Each list
          element is an SVT of depth 0 (as described above), or a NULL if
          the corresponding column is empty (i.e. has no nonzero data).

          For example, the nonzero data of an 8-column sparse matrix will be
          stored in an SVT that looks like this:
          \preformatted{
    .------------------list-of-length-8-----------------.
   /       /       /      |       |      \       \       \
  |       |       |       |       |       |       |       |
 leaf    leaf    NULL    leaf    leaf    leaf    leaf    NULL
vector  vector          vector  vector  vector  vector}

          The NULL leaves represent the empty columns (i.e. the columns
          with no nonzero elements).

    \item An SVT_SparseArray object with 3 dimensions has its nonzero data
          stored in an SVT of depth 2. Such SVT is represented by a list of
          length the extend of the 3rd dimension. Each list element must be
          an SVT of depth 1 (as described above) that stores the nonzero data
          of the corresponding 2D slice, or a NULL if the 2D slice is empty
          (i.e. has no nonzero data).
  }
}

\seealso{
  \itemize{
    \item The \link{SparseArray} class for the virtual parent class of
          COO_SparseArray and SVT_SparseArray.

    \item S4 classes \linkS4class{dgCMatrix} and \linkS4class{lgCMatrix}
          defined in the \pkg{Matrix} package, for the de facto standard
          of sparse matrix representations in R.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
m <- matrix(rpois(3e6, lambda=0.4), ncol=1200)

## Note that 'SparseArray(m)' can also be used for this:
svt <- SVT_SparseArray(m)
svt
type(svt)
object.size(svt)

dgcm <- as(m, "dgCMatrix")
type(dgcm)
object.size(dgcm)
}
\keyword{classes}
\keyword{methods}
