\name{readSparseCSV}

\alias{writeSparseCSV}
\alias{readSparseCSV}
\alias{readSparseTable}

\title{Write/read a sparse matrix from/to a CSV file}

\description{
  Write or read a sparse matrix from or to a CSV (comma-separated values) file.
}

\usage{
writeSparseCSV(x, filepath, sep=",", transpose=FALSE, write.zeros=FALSE)

readSparseCSV(filepath, sep=",", transpose=FALSE,
              as=c("SparseMatrix", "dgCMatrix"))
}

\arguments{
  \item{x}{
    A matrix-like object, typically sparse. IMPORTANT: The object must have
    rownames and colnames! These will be written to the file.

    Another requirement is that the object must be subsettable. More
    precisely: it must support 2D-style subsetting of the kind \code{x[i, ]}
    and \code{x[ , j]} where \code{i} and \code{j} are valid row and column
    indices.
  }
  \item{filepath}{
    The path (as a single string) to the file where to write the
    matrix-like object or to read it from.
    If \code{""}, \code{writeSparseCSV()} will write the data to the
    standard output connection.

    Note that \code{filepath} can also be a connection.
  }
  \item{sep}{
    The field separator character. Values on each line of the
    file are separated by this character.
  }
  \item{transpose}{
    \code{TRUE} or \code{FALSE}. By default, rows in the matrix-like object
    correspond to lines in the CSV file. Set \code{transpose} to \code{TRUE}
    to transpose the matrix-like object on-the-fly, that is, to have its
    columns written to or read from the lines in the CSV file.

    Note that using \code{transpose=TRUE} is semantically equivalent to
    calling \code{t()} on the object before writting it or after reading
    it, but it will tend to be more efficient. Also it will work even if
    \code{x} does not support \code{t()} (not all matrix-like objects are
    guaranteed to be transposables).
  }
  \item{write.zeros}{
    \code{TRUE} or \code{FALSE}. By default, the zero values in \code{x}
    are not written to the file. Set \code{write.zeros} to \code{TRUE} to
    write them.
  }
  \item{as}{
    What to return the sparse data as. By default it's returned as a
    \link{SparseMatrix} object.
  }
}

\value{
  \code{writeSparseCSV} returns an invisible \code{NULL}.

  \code{readSparseCSV} returns a \link{SparseMatrix} object by default.
}

\seealso{
  \itemize{
    \item \link{SparseArray} objects.

    \item \linkS4class{dgCMatrix} objects in the \pkg{Matrix} package.
  }
}

\examples{
## Prepare toy matrix 'm0':
rownames0 <- LETTERS[1:6]
colnames0 <- letters[1:4]
m0 <- matrix(0L, nrow=length(rownames0), ncol=length(colnames0),
                 dimnames=list(rownames0, colnames0))
m0[c(1:2, 8, 10, 15:17, 24)] <- (1:8)*10L
m0

## writeSparseCSV():
writeSparseCSV(m0, filepath="", sep="\t")
writeSparseCSV(m0, filepath="", sep="\t", write.zeros=TRUE)

outfile <- tempfile()
writeSparseCSV(m0, outfile)

## readSparseCSV():
svt1 <- readSparseCSV(outfile)
svt1

svt2 <- readSparseCSV(outfile, transpose=TRUE)
svt2

readSparseCSV(outfile, as="dgCMatrix")

## Sanity checks:
stopifnot(identical(m0, as.matrix(svt1)))
stopifnot(identical(t(m0), as.matrix(svt2)))
}
\keyword{utilities}
