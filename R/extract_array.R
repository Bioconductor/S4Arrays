### =========================================================================
### extract_array()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

### Return the slice as a list.
.extract_data_frame_slice <- function(x, index)
{
    slice <- subset_by_Nindex(x, index)
    ## Turn into a list and replace factors with character vectors.
    lapply(slice, as.vector)
}

.extract_DataFrame_slice0 <- function(x)
{
    ## Make sure that this remains consistent with .get_DataFrame_type()
    ## defined in R/type.R
    x0 <- x[0L, , drop=FALSE]
    df0 <- as.data.frame(x0)
    if (ncol(df0) != ncol(x))
        stop(wmsg("DataFrame object 'x' can be used as the seed of ",
                  "a DelayedArray object only if as.data.frame(x) ",
                  "preserves the number of columns"))
    BiocGenerics:::extract_data_frame_slice0(df0)
}

.extract_DataFrame_slice <- function(x, index)
{
    slice <- subset_by_Nindex(x, index)
    slice <- as.data.frame(slice)
    ## Turn into a list and replace factors with character vectors.
    lapply(slice, as.vector)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array() generic and methods
###
### Note that extract_array() is part of the "seed contract" as defined in
### the "Implementing A DelayedArray Backend" vignette from the DelayedArray
### package.

### Similar to SparseArray:::.contact_author_msg2().
.contact_author_msg1 <- function(.Generic, x_class)
{
    msg <- c("Please contact the authors/maintainers of the ",
             x_class, " class")
    class_package <- attr(x_class, "package")
    if (!is.null(class_package))
        msg <- c(msg, " (defined in the ", class_package, " package)")
    c(msg, " about this, and point them to the man page for the ",
           .Generic, "() generic function defined in the S4Arrays ",
           "package ('?S4Arrays::", .Generic, "').")
}

check_returned_array <- function(ans, expected_dim, .Generic, x_class)
{
    if (!is.array(ans))
        stop(wmsg("The ", .Generic, "() method for ", x_class, " ",
                  "objects didn't return an ordinary array. ",
                  .Generic, "() methods should **always** return an ",
                  "ordinary array. ",
                  .contact_author_msg1(.Generic, x_class)))
    if (!identical(dim(ans), expected_dim))
        stop(wmsg("The ", .Generic, "() method for ", x_class, " ",
                  "objects returned an array with incorrect ",
                  "dimensions. ", .contact_author_msg1(.Generic, x_class)))
    ans
}

### 'index' is expected to be an unnamed list of subscripts as positive
### integer vectors, one vector per dimension in 'x'. *Missing* list elements
### are allowed and represented by NULLs.
### The extract_array() methods don't need to support anything else.
### They must return an ordinary array. No need to propagate the dimnames.
setGeneric("extract_array", signature="x",
    function(x, index)
    {
        x_dim <- dim(x)
        if (is.null(x_dim))
            stop(wmsg("the first argument to extract_array() must be an ",
                      "array-like object (i.e. it must have dimensions)"))
        ans <- standardGeneric("extract_array")
        expected_dim <- get_Nindex_lengths(index, x_dim)
        check_returned_array(ans, expected_dim, "extract_array", class(x))
    }
)

### subset_by_Nindex() uses `[` internally to perform the subsetting, so
### this default extract_sparse_array() method will work on any object 'x'
### that supports `[` and as.array().
setMethod("extract_array", "ANY",
    function(x, index)
    {
        slice <- subset_by_Nindex(x, index)
        as.array(slice)
    }
)

setMethod("extract_array", "array",
    function(x, index) subset_by_Nindex(x, index)
)

### Equivalent to
###
###     subset_by_Nindex(as.matrix(x), index)
###
### but avoids the cost of turning the full data frame 'x' into a matrix so
### memory footprint stays small when 'index' is small.
setMethod("extract_array", "data.frame",
    function(x, index)
    {
        slice0 <- BiocGenerics:::extract_data_frame_slice0(x)
        slice <- .extract_data_frame_slice(x, index)
        data <- unlist(c(slice0, slice), recursive=FALSE, use.names=FALSE)
        array(data, dim=get_Nindex_lengths(index, dim(x)))
    }
)

### Equivalent to
###
###     subset_by_Nindex(as.matrix(as.data.frame(x)), index)
###
### but avoids the cost of turning the full DataFrame 'x' first into a data
### frame then into a matrix so memory footprint stays small when 'index' is
### small.
setMethod("extract_array", "DataFrame",
    function(x, index)
    {
        slice0 <- .extract_DataFrame_slice0(x)
        slice <- .extract_DataFrame_slice(x, index)
        data <- unlist(c(slice0, slice), recursive=FALSE, use.names=FALSE)
        array(data, dim=get_Nindex_lengths(index, dim(x)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A convenience wrapper around extract_array()
###

### An enhanced version of extract_array() that accepts an Nindex (see
### Nindex-utils.R) and propagates the dimnames.
### WARNING: The list elements in 'Nindex' can only be NULLs, integer
### vectors, or RangeNSBS objects at the moment. extract_array_by_Nindex()
### will break if they are not! See FIXME below.
extract_array_by_Nindex <- function(x, Nindex)
{
    ## TODO: Once we have a full Nindex normalization mechanism in place
    ## (see FIXME below), use it to normalize the supplied 'Nindex' in 2
    ## steps: (1) by normalizing with something like 'as.NSBSlist=TRUE'
    ## to produce an Nindex with NSBS list elements, then (2) by doing
    ## something like:
    ##
    ##   lapply( , function(i) if (is.null(i)) NULL else as.integer(i))
    ##
    ## on the Nindex obtained at (1).
    ## Pass the Nindex obtained at (1) to subset_dimnames_by_Nindex() and
    ## the Nindex obtained at (2) to extract_array().
    ans_dimnames <- subset_dimnames_by_Nindex(dimnames(x), Nindex)

    ## FIXME: The list elements of an Nindex can be anything (see
    ## Nindex-utils.R) so it's not enough to expand only those list elements
    ## that are RangeNSBS objects. For example the call to extract_array()
    ## below will fail if some subscripts in 'Nindex' are character vectors
    ## or Rle objects. We need to perform a full normalization of 'Nindex'
    ## like we do in new_DelayedSubset() (see DelayedOp-class.R). Note that
    ## we're good for now because extract_array_by_Nindex() is only used
    ## in the context of show_compact_array() and the default "read_block"
    ## method where the supplied 'Nindex' is guaranteed to contain only
    ## NULLs, integer vectors, or RangeNSBS objects.
    ans <- extract_array(x, expand_Nindex_RangeNSBS(Nindex))

    set_dimnames(ans, ans_dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### as.array(x) (in-memory realization of an array-like object)
###

### Realize the object i.e. execute all the delayed operations and turn the
### object back into an ordinary array.
.from_Array_to_array <- function(x, drop=FALSE)
{
    if (!isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    index <- vector("list", length=length(dim(x)))
    ans <- extract_array(x, index)
    ans <- set_dimnames(ans, dimnames(x))
    if (drop)
        ans <- drop(ans)
    ans
}

### S3/S4 combo for as.array.Array
as.array.Array <- function(x, ...) .from_Array_to_array(x, ...)
setMethod("as.array", "Array", .from_Array_to_array)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other coercions to in-memory representations
###
### All these coercions are based on as.array().
###

### S3/S4 combo for as.data.frame.Array
as.data.frame.Array <- function(x, row.names=NULL, optional=FALSE, ...)
    as.data.frame(as.array(x, drop=TRUE),
                  row.names=row.names, optional=optional, ...)
setMethod("as.data.frame", "Array", as.data.frame.Array)

### S3/S4 combo for as.vector.Array
as.vector.Array <- function(x, mode="any")
{
    ans <- as.array(x, drop=TRUE)
    as.vector(ans, mode=mode)
}
setMethod("as.vector", "Array", as.vector.Array)

### S3/S4 combo for as.logical.Array
as.logical.Array <- function(x, ...) as.vector(x, mode="logical", ...)
setMethod("as.logical", "Array", as.logical.Array)

### S3/S4 combo for as.integer.Array
as.integer.Array <- function(x, ...) as.vector(x, mode="integer", ...)
setMethod("as.integer", "Array", as.integer.Array)

### S3/S4 combo for as.numeric.Array
as.numeric.Array <- function(x, ...) as.vector(x, mode="numeric", ...)
setMethod("as.numeric", "Array", as.numeric.Array)

### S3/S4 combo for as.complex.Array
as.complex.Array <- function(x, ...) as.vector(x, mode="complex", ...)
setMethod("as.complex", "Array", as.complex.Array)

### S3/S4 combo for as.character.Array
as.character.Array <- function(x, ...) as.vector(x, mode="character", ...)
setMethod("as.character", "Array", as.character.Array)

### S3/S4 combo for as.raw.Array
as.raw.Array <- function(x) as.vector(x, mode="raw")
setMethod("as.raw", "Array", as.raw.Array)

