### =========================================================================
### SVT_SparseArray objects
### -------------------------------------------------------------------------
###
### SparseArray objects using the SVT layout.
###
### An SVT_SparseArray object stores its nonzero data in a "Sparse Vector
### Tree" (SVT).
###
### An SVT is a tree of depth the number of dimensions in the array where
### the leaves are sparse vectors, also called "leaf vectors".
### A leaf vector is represented by a list of 2 parallel vectors: an integer
### vector of positions and a vector (atomic or list) of nonzero values.
### The 2nd vector determines the type of the leaf vector. All the leaf
### vectors in the SVT must have the same type, which should match the type
### specified in the 'type' slot of the SVT_SparseArray object.
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
### Going back and forth between SVT_SparseArray and COO_SparseArray objects
###

.from_SVT_SparseArray_to_COO_SparseArray <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    ## Returns 'ans_nzindex' and 'ans_nzdata' in a list of length 2.
    C_ans <- .Call2("C_from_SVT_SparseArray_to_COO_SparseArray",
                    from@dim, from@type, from@SVT, PACKAGE="S4Arrays")
    ans_nzindex <- C_ans[[1L]]
    ans_nzdata  <- C_ans[[2L]]
    new2("COO_SparseArray", dim=from@dim, dimnames=from@dimnames,
                            nzindex=ans_nzindex, nzdata=ans_nzdata,
                            check=FALSE)
}

.from_COO_SparseArray_to_SVT_SparseArray <- function(from)
{
    stopifnot(is(from, "COO_SparseArray"))
    ans_SVT <- .Call2("C_build_SVT_from_COO_SparseArray",
                      from@dim, from@nzindex, from@nzdata, PACKAGE="S4Arrays")
    .new_SVT_SparseArray(from@dim, from@dimnames, type(from), ans_SVT,
                         check=FALSE)
}

setAs("SVT_SparseArray", "COO_SparseArray",
    .from_SVT_SparseArray_to_COO_SparseArray
)
setAs("COO_SparseArray", "SVT_SparseArray",
    .from_COO_SparseArray_to_SVT_SparseArray
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SVT_SparseArray constructor
###

.make_SVT_SparseArray_from_dgCMatrix <- function(x, as.integer=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(as.integer))
        stop(wmsg("'as.integer' must be TRUE or FALSE"))
    ans_SVT <- .Call2("C_build_SVT_from_dgCMatrix",
                      x, as.integer, PACKAGE="S4Arrays")
    .new_SVT_SparseArray(dim(x), dimnames(x), type(x), ans_SVT,
                         check=FALSE)
}

SVT_SparseArray <- function(x, as.integer=FALSE)
{
    if (missing(x))
        return(new2("SVT_SparseArray", check=FALSE))
    if (!is(x, "dgCMatrix"))
        return(as(x, "SVT_SparseArray"))
    .make_SVT_SparseArray_from_dgCMatrix(x, as.integer=as.integer)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray and [d|l]gCMatrix objects
###

setAs("dgCMatrix", "SVT_SparseArray",
    function(from) .make_SVT_SparseArray_from_dgCMatrix(from)
)

setAs("lgCMatrix", "SVT_SparseArray",
    function(from)
    {
        stop("not ready yet")
    }
)

.from_SVT_SparseArray_to_CsparseMatrix <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    if (length(from@dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgCMatrix ",
                  "or lgCMatrix must have exactly 2 dimensions"))
    ans_class <- switch(from@type,
                        'integer'=, 'double'="dgCMatrix",
                        'logical'="lgCMatrix",
                        stop(wmsg("unsupported data type: ", from@type)))
    ## Returns 'ans_p', 'ans_i', and 'ans_x', in a list of length 3.
    C_ans <- .Call2("C_from_SVT_SparseArray_to_CsparseMatrix",
                    from@dim, from@type, from@SVT, PACKAGE="S4Arrays")
    ans_p <- C_ans[[1L]]
    ans_i <- C_ans[[2L]]
    ans_x <- C_ans[[3L]]
    new(ans_class, Dim=from@dim, p=ans_p, i=ans_i, x=ans_x,
                   Dimnames=from@dimnames)
}

setAs("SVT_SparseArray", "dgCMatrix", .from_SVT_SparseArray_to_CsparseMatrix)
setAs("SVT_SparseArray", "lgCMatrix", .from_SVT_SparseArray_to_CsparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray objects and ordinary arrays
###

.from_SVT_SparseArray_to_array <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    .Call2("C_from_SVT_SparseArray_to_Rarray",
           from@dim, dimnames(from), from@type, from@SVT, PACKAGE="S4Arrays")
}

.from_array_to_SVT_SparseArray <- function(from)
{
    stopifnot(is.array(from))
    ans_SVT <- .Call2("C_build_SVT_from_Rarray", from, PACKAGE="S4Arrays")
    .new_SVT_SparseArray(dim(from), dimnames(from), type(from), ans_SVT,
                         check=FALSE)
}

### S3/S4 combo for as.array.SVT_SparseArray
as.array.SVT_SparseArray <- function(x, ...) .from_SVT_SparseArray_to_array(x)
setMethod("as.array", "SVT_SparseArray", as.array.SVT_SparseArray)

setAs("array", "SVT_SparseArray", .from_array_to_SVT_SparseArray)

