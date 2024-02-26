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

### --- OLD read_block() behavior (BioC <= 3.17) ---

.load_DelayedArray_for_read_block <- function(...)
    load_package_with_graceful_failure("DelayedArray",
                                       "calling read_block() ", ...)

### Provides the old read_block() behavior (used in BioC <= 3.17) where
### a sparse block gets returned as a SparseArraySeed object from the
### DelayedArray package.
.OLD_read_block <- function(x, viewport, as.sparse=NA)
{
    if (is_sparse(x)) {
        .load_DelayedArray_for_read_block("on a ", class(x), " object ")
        ans <- DelayedArray::read_sparse_block(x, viewport)
        if (isFALSE(as.sparse))
            ans <- DelayedArray::sparse2dense(ans)
    } else {
        ans <- read_block_as_dense(x, viewport)
        check_returned_array(ans, dim(viewport),
                             "read_block_as_dense", class(x))
        if (isTRUE(as.sparse)) {
            .load_DelayedArray_for_read_block("with 'as.sparse=TRUE'")
            ans <- DelayedArray::dense2sparse(ans)
        }
    }
    ans
}

### --- NEW read_block() behavior (BioC >= 3.18) ---

.load_SparseArray_for_read_block <- function(...)
    load_package_with_graceful_failure("SparseArray",
                                       "calling read_block() ", ...)

### Provides the new read_block() behavior (used in BioC >= 3.18) where
### a sparse block gets returned as a SparseArray object from the
### new SparseArray package. Note that this new behavior makes use of
### the new SparseArray::read_block_as_sparse() generic (replaces
### DelayedArray::read_sparse_block()).
.NEW_read_block <- function(x, viewport, as.sparse=NA)
{
    if (is_sparse(x)) {
        .load_SparseArray_for_read_block("on a ", class(x), " object ")
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
            .load_SparseArray_for_read_block("with 'as.sparse=TRUE'")
            ans <- as(ans, "SparseArray")
        }
    }
    ans
}

### A user-facing frontend for read_block_as_dense() and
### SparseArray::read_block_as_sparse().
### Reads a block of data from array-like object 'x'. Depending on the value
### of argument 'as.sparse', the block is returned either as an ordinary
### array (dense representation) or a SparseArray object (sparse
### representation).
### 'as.sparse' can be TRUE, FALSE, or NA. If FALSE, the block is returned
### as an ordinary array. If TRUE, it's returned as a SparseArray object.
### Using 'as.sparse=NA' (the default) is equivalent to
### using 'as.sparse=is_sparse(x)'. This is the most efficient way to read
### a block.
### Propagate the dimnames.
read_block <- function(x, viewport, as.sparse=NA)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("the first argument to read_block() must be an ",
                  "array-like object (i.e. it must have dimensions)"))
    stopifnot(is(viewport, "ArrayViewport"),
              identical(refdim(viewport), x_dim),
              is.logical(as.sparse),
              length(as.sparse) == 1L)

    ## IMPORTANT NOTE: We temporarily preserve the old read_block() behavior
    ## for backward compatibility. See comments for .OLD_read_block()
    ## and .NEW_read_block() above for additional details.
    ## TODO: In BioC 3.18, the plan is to switch to the new behavior, and
    ## to update man/read_block.Rd accordingly. See TODO file in DelayedArray
    ## for the details.
    ans <- .OLD_read_block(x, viewport, as.sparse=as.sparse)
    #ans <- .NEW_read_block(x, viewport, as.sparse=as.sparse)

    ## Individual read_block_as_dense() and read_block_as_sparse() methods
    ## are not expected to propagate the dimnames so we take care of this
    ## now.
    Nindex <- makeNindexFromArrayViewport(viewport)
    ans_dimnames <- subset_dimnames_by_Nindex(dimnames(x), Nindex)
    set_dimnames(ans, ans_dimnames)
}

