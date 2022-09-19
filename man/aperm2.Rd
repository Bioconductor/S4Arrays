\name{aperm2}

\alias{aperm2}

\title{Generalized permutation of the dimensions of an array}

\description{
  \code{aperm2()} extends the functionality of
  \code{base::\link[base]{aperm}()} by allowing dropping and/or
  adding \emph{ineffective dimensions} (i.e. dimensions with an extent of 1)
  from/to the input array, in addition to permuting its dimensions.

  Note that, like \code{base::\link[base]{aperm}()}, \code{aperm2()}
  always preserves the length of the input array. However, unlike with
  \code{base::\link[base]{aperm}()}, the array returned by \code{aperm2()}
  doesn't necessarily have the same number of dimensions as the input array.
}

\usage{
aperm2(a, perm)
}

\arguments{
  \item{a}{
    An array.
  }
  \item{perm}{
    An integer vector, possibly containing \code{NA}s, indicating how
    the dimensions of the returned array should be mapped to those of
    the input array.

    More precisely, \code{perm} can be one of the following:
    \itemize{
      \item A permutation of the \code{seq_along(dim(a))} vector,
            like for \code{base::\link[base]{aperm}()}.
            Note that if the identity permutation is used (i.e.
            \code{perm=seq_along(dim(a))}), then \code{aperm2()}
            is a no-op (like with \code{base::\link[base]{aperm}()}).

      \item A permutation of a \emph{subset} of the \code{seq_along(dim(a))}
            vector. In this case the dimensions that are excluded
            must be \emph{ineffective dimensions} i.e. each of them
            must have an extent of 1. In other words, only integers
            that belong to \code{which(dim(a) == 1)} can be missing
            from \code{perm}.

            In this case, the ineffective dimensions that are excluded
            will be dropped i.e. they won't be carried over to the
            returned array.

      \item Additionally, any number of \code{NA}s can be inserted
            anywhere in a \code{perm} vector like one described above.

            In this case, ineffective dimensions will be added to the
            returned array. These added dimensions will materialize
            as additional \code{1}'s in the \code{dim()} vector of the
            returned array, at positions that match the positions of
            the \code{NA}s in \code{perm}.
    }

    Note that if \code{perm} is missing, then \code{aperm2(a)} reverses
    the order of \code{a}'s dimensions (i.e. \code{perm} gets set to
    \code{rev(seq_along(dim(a)))}), like \code{base::\link[base]{aperm}(a)}
    does.
  }
}

\value{
An array with one dimension per element in the \code{perm} argument.
The length of the returned array will always be the same as the length
of the input array. (Note that for an array \code{a}, \code{length(a)}
is \code{prod(dim(a))}.)
}

\note{
The \code{aperm()} method for \link[DelayedArray]{DelayedArray} objects
defined in the \pkg{DelayedArray} package implements the "\code{aperm2}
semantic", that is, it allows dropping and/or adding \emph{ineffective
dimensions} from/to the input \link[DelayedArray]{DelayedArray} object.
}

\seealso{
  \itemize{
    \item \code{\link[base]{aperm}} in the \pkg{base} package for
          the function that \code{aperm2} is based on.

    \item \code{\link[BiocGenerics]{aperm}} in the \pkg{BiocGenerics}
          package for the \code{aperm} \emph{S4 generic function}.

    \item \link[DelayedArray]{aperm,DelayedArray-method} in the
          \pkg{DelayedArray} package for an \code{aperm()} method
          that implements the "\code{aperm2} semantic".
  }
}

\examples{
## ---------------------------------------------------------------------
## SOME EXAMPLES WITH A 4D ARRAY
## ---------------------------------------------------------------------

a <- array(1:72, c(3, 6, 1, 4),
           dimnames=list(NULL, letters[1:6], NULL, LETTERS[1:4]))
a

## Permute first two dimensions:
aperm2(a, perm=c(2,1,3,4))

## Permute first and last dimensions:
aperm2(a, perm=c(4,2,3,1))

## Drop 3rd dimension:
aperm2(a, perm=c(1,2,4))

## Drop 3rd dimension and permute 2nd and last:
aperm2(a, perm=c(1,4,2))

## Drop 3rd dimension and cycle the order of the remaining ones:
aperm2(a, perm=c(2,4,1))

## Add one ineffective dimension:
aperm2(a, perm=c(NA,1,2,3,4))
aperm2(a, perm=c(1,NA,2,3,4))
aperm2(a, perm=c(1,2,NA,3,4))
aperm2(a, perm=c(1,2,3,NA,4))
aperm2(a, perm=c(1,2,3,4,NA))

## Add four ineffective dimensions:
aperm2(a, perm=c(NA,1,2,3,NA,NA,4,NA))

## Permute first and last dimensions and add one ineffective dimension:
aperm2(a, perm=c(4,2,3,NA,1))

## Drop 3rd dimension, cycle the order of the remaining ones, and add
## two ineffective dimensions:
aperm2(a, perm=c(2,4,NA,1,NA))

## No-op:
aperm2(a, perm=seq_along(dim(a)))

## Reverse the order of the dimensions (multidimensional transposition):
aperm2(a)  # same as 'aperm2(a, perm=rev(seq_along(dim(a))))'

## ---------------------------------------------------------------------
## COMPOSING aperm2() TRANSFORMATIONS
## ---------------------------------------------------------------------

## Applying two successive aperm() transformations, first with 'perm'
## set to 'perm1' then set to 'perm2', is equivalent to applying a
## single aperm() transformation with 'perm' set to 'perm1[perm2]'.
##
## More formally:
##     aperm(aperm(a, perm=perm1), perm=perm2)
## is equivalent to:
##     aperm(a, perm=perm1[perm2])
##
## Note that this also applies to aperm2()!

## Examples with aperm():

perm1 <- c(2,4,3,1)
perm2 <- c(4,3,2,1)
perm3 <- c(2,1,4,3)

a12 <- aperm(aperm(a, perm=perm1), perm=perm2)
stopifnot(identical(a12, aperm(a, perm=perm1[perm2])))

a13 <- aperm(aperm(a, perm=perm1), perm=perm3)
stopifnot(identical(a13, aperm(a, perm=perm1[perm3])))

a23 <- aperm(aperm(a, perm=perm2), perm=perm3)
stopifnot(identical(a23, aperm(a, perm=perm2[perm3])))

a123 <- aperm(aperm(aperm(a, perm=perm1), perm=perm2), perm=perm3)
stopifnot(identical(a123, aperm(a, perm=perm1[perm2][perm3])))
stopifnot(identical(a123, aperm(a, perm=perm1[perm2[perm3]])))

## Examples with aperm2():

perm1 <- c(2,4,1)
perm2 <- c(1,3,NA,2,NA)
perm3 <- c(5,4,2,1)

a12 <- aperm2(aperm2(a, perm=perm1), perm=perm2)
stopifnot(identical(a12, aperm2(a, perm=perm1[perm2])))

a123 <- aperm2(aperm2(aperm2(a, perm=perm1), perm=perm2), perm=perm3)
stopifnot(identical(a123, aperm2(a, perm=perm1[perm2][perm3])))
stopifnot(identical(a123, aperm2(a, perm=perm1[perm2[perm3]])))

## ---------------------------------------------------------------------
## REVERSIBILITY OF THE aperm2() TRANSFORMATION
## ---------------------------------------------------------------------

## An aperm() or aperm2() transformation is always reversible.
## The 'perm' vector to use to achieve the reverse transformation
## can be inferred from the initial 'perm' vector using the following
## helper function ('n' must be the number of dimensions of
## the original array):
build_rev_perm <- function(perm, n=length(perm)) {
    rev_perm <- rep.int(NA_integer_, n)
    na_idx <- which(!is.na(perm))
    rev_perm[perm[na_idx]] <- na_idx
    rev_perm
}

## Examples:

perm <- c(2,4,NA,1,NA)
rev_perm <- build_rev_perm(perm, n=length(dim(a)))
stopifnot(identical(aperm2(aperm2(a, perm=perm), perm=rev_perm), a))

## The "composed" 'perm' vector achieves identity:
perm[rev_perm]

## Sanity checks:

perm <- seq_len(10)
stopifnot(identical(build_rev_perm(perm), perm))

perm <- c(2:5,1L)
rev_perm <- build_rev_perm(perm)
stopifnot(identical(perm[rev_perm], seq_along(perm)))

perm <- c(5L,NA,2:4,NA,NA,1L)
rev_perm <- build_rev_perm(perm, n=6)
stopifnot(identical(perm[rev_perm], c(1:5,NA)))
}
\keyword{array}
\keyword{manip}
