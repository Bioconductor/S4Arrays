### =========================================================================
### SVT_SparseArray objects
### -------------------------------------------------------------------------
###
### Use SVT layout to store the sparse data.
###
### An SVT_SparseArray object stores its nonzero data in a "Sparse Vector
### Tree" (SVT).
###
### An SVT is a tree of depth the number of dimensions in the array where
### the leaves are sparse vectors, also called "leaf vectors".
### A "leaf vector" is a vector of offset/value pairs sorted by strictly
### ascending offset. It is represented by a list of 2 parallel vectors:
### an integer vector of offsets and a vector (atomic or list) of nonzero
### values. The 2nd vector determines the type of the leaf vector i.e.
### double, integer, logical, etc... All the leaf vectors in the SVT must
### have the same type, which should match the type specified in the 'type'
### slot of the SVT_SparseArray object.
###
### More precisely:
###
### - An SVT_SparseArray object with 1 dimension stores its nonzero data in an
###   SVT of depth 1. Such SVT is represented by a single "leaf vector".
###
### - An SVT_SparseArray object with 2 dimensions stores its nonzero data in an
###   SVT of depth 2. Such SVT is represented by a list of length the extend
###   of the 2nd dimension (nb of columns). Each list element is an SVT of
###   depth 1 (as described above), or a NULL if the corresponding column is
###   empty (i.e. has no nonzero data).
###
### - An SVT_SparseArray object with 3 dimensions stores its nonzero data in an
###   SVT of depth 3. Such SVT is represented by a list of length the extend
###   of the 3rd dimension. Each list element must be an SVT of depth 2 (as
###   described above) that stores the nonzero data of the corresponding 2D
###   slice, or a NULL if the 2D slice is empty (i.e. has no nonzero data).
###
### - etc...
###
### If the sparse array is empty (i.e. has no nonzero data), the 'SVT' slot
### is set to NULL.
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

### Supports SVT_SparseArray objects of any atomic type except "character".
.from_SVT_SparseArray_to_dgCMatrix <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    if (length(from@dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgCMatrix ",
                  "must have exactly 2 dimensions"))

    ## If 'from@type' is "integer" or "raw", we'll coerce 'ans_x' to "double"
    ## right before passing it to new_CsparseMatrix() below. This is ok
    ## because it won't introduce zeros in 'ans_x'. Also it should be slightly
    ## more efficient than switching the type of 'from' now.
    ## However, if 'from@type' is "complex", we need to do the switch now,
    ## because coercing 'ans_x' to "double" would potentially introduce
    ## zeros in it.
    if (from@type == "complex")
        type(from) <- "double"

    ## Returns 'ans_p', 'ans_i', and 'ans_x', in a list of length 3.
    C_ans <- .Call2("C_from_SVT_SparseArray_to_CsparseMatrix",
                    from@dim, from@type, from@SVT, PACKAGE="S4Arrays")
    ans_p <- C_ans[[1L]]
    ans_i <- C_ans[[2L]]
    ans_x <- C_ans[[3L]]  # same type as 'from'
    if (storage.mode(ans_x) != "double")
        storage.mode(ans_x) <- "double"  # won't introduce zeros
    new_CsparseMatrix(from@dim, ans_p, ans_i, ans_x, dimnames=from@dimnames)
}

### Supports SVT_SparseArray objects of any atomic type except "character".
.from_SVT_SparseArray_to_lgCMatrix <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    if (length(from@dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to lgCMatrix ",
                  "must have exactly 2 dimensions"))

    ## Returns 'ans_p', 'ans_i', and 'ans_x', in a list of length 3.
    C_ans <- .Call2("C_from_SVT_SparseArray_to_CsparseMatrix",
                    from@dim, from@type, from@SVT, PACKAGE="S4Arrays")
    ans_p <- C_ans[[1L]]
    ans_i <- C_ans[[2L]]
    ans_x <- C_ans[[3L]]  # same type as 'from'
    if (storage.mode(ans_x) != "logical")
        storage.mode(ans_x) <- "logical"  # won't introduce FALSEs
    new_CsparseMatrix(from@dim, ans_p, ans_i, ans_x, dimnames=from@dimnames)
}

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
### Constructor
###

SVT_SparseArray <- function(x, type=NA)
{
    if (!identical(type, NA))
        type <- normarg_array_type(type, "the requested type")

    if (missing(x)) {
        if (identical(type, NA))
            type <- "logical"
        return(new2("SVT_SparseArray", type=type, check=FALSE))
    }

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

