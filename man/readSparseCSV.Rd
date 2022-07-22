\name{readSparseCSV}
  
\alias{readSparseCSV}
\alias{readSparseTable}

\title{An alternative to \code{scuttle::readSparseCounts}}

\description{
  THIS IS A WORK-IN-PROGRESS!
}

\usage{
readSparseCSV(filepath, sep=",", transpose=FALSE,
              as=c("SparseMatrix", "dgCMatrix"))
}

\arguments{
  \item{filepath}{
    COMING SOON...
  }
  \item{sep}{
    COMING SOON...
  }
  \item{transpose}{
    COMING SOON...
  }
  \item{as}{
    COMING SOON...
  }
}

\details{
  COMING SOON...
}

\value{
  COMING SOON...
}

\seealso{
  \itemize{
    \item \code{\link[scuttle]{readSparseCounts}} in the \pkg{scuttle} package.
  }
}

\examples{
library(scuttle)
example(readSparseCounts)  # to create 'outfile'

readSparseCSV(outfile, sep="\t")
readSparseCSV(outfile, sep="\t", transpose=TRUE)

readSparseCSV(outfile, sep="\t", as="dgCMatrix")
}
\keyword{utilities}
