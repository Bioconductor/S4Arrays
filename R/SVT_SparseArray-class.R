### =========================================================================
### SVT_SparseArray objects
### -------------------------------------------------------------------------
###
### Use SVT layout to store the sparse data.
###
### An SVT_SparseArray object stores its nonzero data in a "Sparse Vector
### Tree" (SVT).
###
### If the sparse array is empty (i.e. has no nonzero data), the 'SVT' slot
### must be set to NULL.
###
### IMPORTANT NOTES:
### - All the "leaf vectors" in the SVT are guaranteed to have a
###   length <= the first dimension of the SVT_SparseArray object, which
###   itself is guaranteed to be <= INT_MAX (2^31 - 1).
### - The cumulated length of the "leaf vectors" in the SVT is the number
###   of nonzero values (i.e. nzdata length) in the SVT_SparseArray object.
###   There is no upper limit to this number.
###   In other words, unlike dgCMatrix objects where this number is
###   limited to INT_MAX, an SVT_SparseArray can store an arbitrary number
###   of nonzero values.
###

setClassUnion("NULL_OR_list", c("NULL", "list"))

setClass("SVT_SparseArray",
    contains="SparseArray",
    representation(
        type="character",
        SVT="NULL_OR_list"  # NULL or Sparse Vector Tree (SVT)
    ),
    prototype(
        type="logical"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor
###

.new_SVT_SparseArray <- function(dim, dimnames=NULL, type="logical",
                                 SVT=NULL, check=TRUE)
{
    dimnames <- normarg_dimnames(dimnames, dim)
    new2("SVT_SparseArray", dim=dim, dimnames=dimnames, type=type,
                            SVT=SVT, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("type", "SVT_SparseArray", function(x) x@type)

### Note that like for the length of atomic vectors in base R, the returned
### length will be a double if it's > .Machine$integer.max
.get_SVT_SparseArray_nzdata_length <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    .Call2("C_get_SVT_SparseArray_nzdata_length",
           x@dim, x@SVT, PACKAGE="S4Arrays")
}

#setMethod("nzdataLength", "SVT_SparseArray",
#    .get_SVT_SparseArray_nzdata_length
#)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() setter
###

.set_SVT_SparseArray_type <- function(x, value)
{
    stopifnot(is(x, "SVT_SparseArray"))

    value <- normarg_array_type(value, "the supplied type")
    x_type <- type(x)
    if (value == x_type)
        return(x)

    new_SVT <- .Call2("C_set_SVT_SparseArray_type",
                      x@dim, x@type, x@SVT, value, PACKAGE="S4Arrays")

    BiocGenerics:::replaceSlots(x, type=value, SVT=new_SVT, check=FALSE)
}

setReplaceMethod("type", "SVT_SparseArray", .set_SVT_SparseArray_type)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray objects and ordinary arrays
###

.from_SVT_SparseArray_to_array <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    .Call2("C_from_SVT_SparseArray_to_Rarray",
           from@dim, dimnames(from), from@type, from@SVT, PACKAGE="S4Arrays")
}

### S3/S4 combo for as.array.SVT_SparseArray
as.array.SVT_SparseArray <- function(x, ...) .from_SVT_SparseArray_to_array(x)
setMethod("as.array", "SVT_SparseArray", as.array.SVT_SparseArray)

.build_SVT_SparseArray_from_array <- function(x, type=NA)
{
    stopifnot(is.array(x))
    if (identical(type, NA))
        type <- type(x)
    ans_SVT <- .Call2("C_build_SVT_from_Rarray",
                      x, type, PACKAGE="S4Arrays")
    .new_SVT_SparseArray(dim(x), dimnames(x), type, ans_SVT, check=FALSE)
}

setAs("array", "SVT_SparseArray",
    function(from) .build_SVT_SparseArray_from_array(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray and [d|l]gCMatrix objects
###

.make_CsparseMatrix_from_SVT_SparseArray <- function(from, to_type)
{
    stopifnot(is(from, "SVT_SparseArray"))
    if (length(from@dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgCMatrix ",
                  "or lgCMatrix must have exactly 2 dimensions"))

    ## Coercion to dgCMatrix (i.e. to_type="double"):
    ## If 'from@type' is "logical", "integer", or "raw", we'll coerce 'ans_x'
    ## to "double" right before passing it to new_CsparseMatrix() below. This
    ## is ok because it won't introduce zeros in 'ans_x'. Also it should be
    ## slightly more efficient than switching the type of 'from' now.
    ## However, if the coercion to "double" can potentially introduce zeros
    ## (e.g. if 'from@type' is "complex"), then we need to switch the type now.
    ## Otherwise we will end up with zeros in the "x" slot of the resulting
    ## dgCMatrix object.

    ## Coercion to lgCMatrix (i.e. to_type="logical"):
    ## If 'from@type' is "integer", "double", "complex", or "raw", we'll
    ## coerce 'ans_x' to "logical" right before passing it to
    ## new_CsparseMatrix() below. This is ok because it won't introduce
    ## logical zeros (i.e. FALSEs) in 'ans_x'. Also it should be slightly
    ## more efficient than switching the type of 'from' now.
    ## However, if the coercion to "logical" can potentially introduce zeros
    ## (e.g. if 'from@type' is "character"), then we need to switch the type
    ## now. Otherwise we will end up with zeros in the "x" slot of the
    ## resulting lgCMatrix object.

    postpone <- coercion_can_introduce_zeros(from@type, to_type)
    if (!postpone)
        type(from) <- to_type  # early type switching

    ## Returns 'ans_p', 'ans_i', and 'ans_x', in a list of length 3.
    C_ans <- .Call2("C_from_SVT_SparseArray_to_CsparseMatrix",
                    from@dim, from@type, from@SVT, PACKAGE="S4Arrays")
    ans_p <- C_ans[[1L]]
    ans_i <- C_ans[[2L]]
    ans_x <- C_ans[[3L]]  # same type as 'from'

    ## This type switching is safe only if it does not introduce zeros.
    if (postpone)
        storage.mode(ans_x) <- to_type  # late type switching

    new_CsparseMatrix(from@dim, ans_p, ans_i, ans_x, dimnames=from@dimnames)
}

.from_SVT_SparseArray_to_dgCMatrix <- function(from)
    .make_CsparseMatrix_from_SVT_SparseArray(from, "double")

.from_SVT_SparseArray_to_lgCMatrix <- function(from)
    .make_CsparseMatrix_from_SVT_SparseArray(from, "logical")

setAs("SVT_SparseArray", "dgCMatrix", .from_SVT_SparseArray_to_dgCMatrix)
setAs("SVT_SparseArray", "lgCMatrix", .from_SVT_SparseArray_to_lgCMatrix)

.build_SVT_SparseArray_from_CsparseMatrix <- function(x, type=NA)
{
    stopifnot(is(x, "CsparseMatrix"))
    if (identical(type, NA))
        type <- type(x)
    ans_SVT <- .Call2("C_build_SVT_from_CsparseMatrix",
                      x, type, PACKAGE="S4Arrays")
    .new_SVT_SparseArray(dim(x), dimnames(x), type, ans_SVT, check=FALSE)
}

setAs("CsparseMatrix", "SVT_SparseArray",
    function(from) .build_SVT_SparseArray_from_CsparseMatrix(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray and COO_SparseArray objects
###

.from_SVT_SparseArray_to_COO_SparseArray <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    ## Returns 'ans_nzcoo' and 'ans_nzvals' in a list of length 2.
    C_ans <- .Call2("C_from_SVT_SparseArray_to_COO_SparseArray",
                    from@dim, from@type, from@SVT, PACKAGE="S4Arrays")
    ans_nzcoo <- C_ans[[1L]]
    ans_nzvals <- C_ans[[2L]]
    new2("COO_SparseArray", dim=from@dim, dimnames=from@dimnames,
                            nzcoo=ans_nzcoo, nzvals=ans_nzvals,
                            check=FALSE)
}

setAs("SVT_SparseArray", "COO_SparseArray",
    .from_SVT_SparseArray_to_COO_SparseArray
)

.build_SVT_SparseArray_from_COO_SparseArray <- function(x, type=NA)
{
    stopifnot(is(x, "COO_SparseArray"))
    if (identical(type, NA)) {
        type <- type(x)
    } else {
        ## Some quick testing/benchmarking seemed to suggest that it's
        ## slightly more efficient to change the type of the input
        ## COO_SparseArray object than that of the output SVT_SparseArray
        ## object.
        type(x) <- type
    }
    ans_SVT <- .Call2("C_build_SVT_from_COO_SparseArray",
                      x@dim, x@nzcoo, x@nzvals, type,
                      PACKAGE="S4Arrays")
    .new_SVT_SparseArray(x@dim, x@dimnames, type, ans_SVT, check=FALSE)
}

setAs("COO_SparseArray", "SVT_SparseArray",
    function(from) .build_SVT_SparseArray_from_COO_SparseArray(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SVT_SparseArray() constructor
###

.new_empty_SVT_SparseArray <- function(type=NA)
{
    if (identical(type, NA))
        type <- "logical"
    new2("SVT_SparseArray", type=type, check=FALSE)
}

.SVT_SparseArray <- function(x, type=NA)
{
    if (is.array(x))
        return(.build_SVT_SparseArray_from_array(x, type=type))

    if (is(x, "CsparseMatrix"))
        return(.build_SVT_SparseArray_from_CsparseMatrix(x, type=type))

    if (is(x, "COO_SparseArray"))
        return(.build_SVT_SparseArray_from_COO_SparseArray(x, type=type))

    ans <- as(x, "SVT_SparseArray")
    if (!identical(type, NA))
        type(ans) <- type
    ans
}

SVT_SparseArray <- function(x, type=NA)
{
    if (!identical(type, NA))
        type <- normarg_array_type(type, "the requested type")

    if (missing(x))
        return(.new_empty_SVT_SparseArray(type))

    .SVT_SparseArray(x, type=type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to SparseArray, and the SparseArray() constructor
###
### The SVT_SparseArray representation is preferred over COO_SparseArray when
### coercing to SparseArray or when using the SparseArray() constructor. With
### one exception: when the object to coerce is an RsparseMatrix derivative!
### This is because we don't have an efficient way to coerce these objects to
### SVT_SparseArray at the moment.

setAs("ANY", "SparseArray", function(from) as(from, "SVT_SparseArray"))
setAs("RsparseMatrix", "SparseArray", function(from) as(from, "COO_SparseArray"))

SparseArray <- function(x, type=NA)
{
    if (!identical(type, NA))
        type <- normarg_array_type(type, "the requested type")

    if (missing(x))
        return(.new_empty_SVT_SparseArray(type))

    if (is(x, "RsparseMatrix")) {
        ans <- as(x, "COO_SparseArray")
        if (!identical(type, NA))
            type(ans) <- type
        return(ans)
    }

    .SVT_SparseArray(x, type=type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transposition
###

### S3/S4 combo for t.SVT_SparseArray
t.SVT_SparseArray <- function(x)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("the ", class(x), " object to transpose ",
                  "must have exactly 2 dimensions"))

    new_SVT <- .Call2("C_transpose_SVT_SparseArray",
                      x_dim, x@type, x@SVT, PACKAGE="S4Arrays")

    BiocGenerics:::replaceSlots(x, dim=rev(x_dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}
setMethod("t", "SVT_SparseArray", t.SVT_SparseArray)

