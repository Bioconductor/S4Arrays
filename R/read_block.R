### =========================================================================
### read_block()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The read_block_as_dense() generic
###
### Not meant to be called directly by the end user. They should call
### higher-level user-facing read_block() function instead, with
### the 'as.sparse' argument set to FALSE.
### Must return an ordinary array.
### Note that the read_block() frontend will take care of propagating the
### dimnames, so, for the sake of efficiency, individual methods should not
### try to do it.

setGeneric("read_block_as_dense", signature="x",
    function(x, viewport) standardGeneric("read_block_as_dense")
)

### This default read_block_as_dense() method will work on any object 'x'
### that supports extract_array() e.g. an ordinary array, a sparseMatrix
### derivative from the Matrix package, a SparseArray derivative,
### a DelayedArray object, a DelayedOp object, an HDF5ArraySeed object, etc...
### Does NOT propagate the dimnames.
setMethod("read_block_as_dense", "ANY",
    function(x, viewport)
    {
        Nindex <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
        extract_array(x, Nindex)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_block()
###

.load_SparseArray <- function(...)
{
    if (!requireNamespace("SparseArray", quietly=TRUE))
        stop("Couldn't load the SparseArray package.\n\n  ",
             wmsg("Note that calling read_block() ", ..., " requires ",
                  "the SparseArray package. Please install it with:"),
             "\n\n    BiocManager::install(\"SparseArray\")")
}

### A user-facing frontend for read_block_as_dense() and read_block_as_sparse().
### Reads a block of data from array-like object 'x'. Depending on the value
### of argument 'as.sparse', the block is returned either as an ordinary
### array (dense representation) or a SparseArray object (sparse
### representation).
### 'as.sparse' can be TRUE, FALSE, or NA. If FALSE, the block is returned
### as an ordinary array. If TRUE, it's returned as a SparseArray object.
### Using 'as.sparse=NA' (the default) is equivalent to
### using 'as.sparse=is_sparse(x)'. This is the most efficient
### way to read a block.
### Must propagate the dimnames.
read_block <- function(x, viewport, as.sparse=NA)
{
    stopifnot(is(viewport, "ArrayViewport"),
              identical(refdim(viewport), dim(x)),
              is.logical(as.sparse),
              length(as.sparse) == 1L)
    if (is_sparse(x)) {
        .load_SparseArray("on a ", class(x), " object ")
        ans <- SparseArray::read_block_as_sparse(x, viewport)
        SparseArray:::check_returned_SparseArray(
                             ans, dim(viewport),
                             "read_block_as_sparse", class(x))
        if (isFALSE(as.sparse))
            ans <- as.array(ans)
    } else {
        ans <- read_block_as_dense(x, viewport)
        check_returned_array(ans, dim(viewport),
                             "read_block_as_dense", class(x))
        if (isTRUE(as.sparse)) {
            .load_SparseArray("with 'as.sparse=TRUE'")
            ans <- as(ans, "SparseArray")
        }
    }
    ## Individual read_block_as_dense() and read_block_as_sparse() methods
    ## are not supposed to propagate the dimnames. We take care of this now.
    Nindex <- makeNindexFromArrayViewport(viewport)
    ans_dimnames <- subset_dimnames_by_Nindex(dimnames(x), Nindex)
    set_dimnames(ans, ans_dimnames)
}

