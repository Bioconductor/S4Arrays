### =========================================================================
### SVTSparseArray objects
### -------------------------------------------------------------------------
###
### An SVTSparseArray object stores its nonzero data in a "Sparse Vector
### Tree" (SVT).
###
### An SVT is a tree of depth the number of dimensions in the array where
### the leaves are sparse vectors, also called "leaf vectors".
### A leaf vector is represented by a list of 2 parallel vectors: an integer
### vector of positions and a vector (atomic or list) of nonzero values.
### The 2nd vector determines the type of the leaf vector. All the leaf
### vectors in the SVT must have the same type, which should match the type
### specified in the 'type' slot of the SVTSparseArray object.
###
### More precisely:
###
### - An SVTSparseArray object with 1 dimension stores its nonzero data in an
###   SVT of depth 1. Such SVT is represented by a single "leaf vector".
###
### - An SVTSparseArray object with 2 dimensions stores its nonzero data in an
###   SVT of depth 2. Such SVT is represented by a list of length the extend
###   of the 2nd dimension (nb of columns). Each list element is an SVT of
###   depth 1 (as described above), or a NULL if the corresponding column is
###   empty (i.e. has no nonzero data).
###
### - An SVTSparseArray object with 3 dimensions stores its nonzero data in an
###   SVT of depth 3. Such SVT is represented by a list of length the extend
###   of the 3rd dimension. Each list element must be an SVT of depth 2 (as
###   described above) that stores the nonzero data of the corresponding 2D
###   slice, or a NULL if the 2D slice is empty (i.e. has no nonzero data).
###
### - etc...
###
### If the sparse array is empty (i.e. has no nonzero data), the 'svtree'
### slot is set to NULL.
###
### IMPORTANT NOTES:
### - All the "leaf vectors" in the SVT are guaranteed to have a
###   length <= the first dimension of the SVTSparseArray object, which
###   itself is guaranteed to be <= INT_MAX (2^31 - 1).
### - The cumulated length of the "leaf vectors" in the SVT is the number
###   of nonzero values (i.e. nzdata length) in the SVTSparseArray object.
###   There is no upper limit to this number.
###   In other words, unlike dgCMatrix objects where this number is
###   limited to INT_MAX, an SVTSparseArray can store an arbitrary number
###   of nonzero values.
###

setClassUnion("NULL_OR_list", c("NULL", "list"))

setClass("SVTSparseArray",
    contains="SparseArray",
    representation(
        type="character",
        svtree="NULL_OR_list"  # NULL or Sparse Vector Tree (SVT)
    ),
    prototype(
        type="logical"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("type", "SVTSparseArray", function(x) x@type)

### Note that like for the length of atomic vectors in base R, the returned
### length will be a double if it's > .Machine$integer.max
.get_SVTSparseArray_nzdata_length <- function(x)
    .Call2("C_get_SVTSparseArray_nzdata_length", x@dim, x@svtree,
                                                 PACKAGE="S4Arrays")

#setMethod("nzdataLength", "SVTSparseArray",
#    .get_SVTSparseArray_nzdata_length
#)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVTSparseArray and COOSparseArray objects
###

.from_SVTSparseArray_to_COOSparseArray <- function(from)
{
    ## Returns 'ans_nzindex' and 'ans_nzdata' in a list of length 2.
    C_ans <- .Call2("C_from_SVTSparseArray_to_COOSparseArray",
                    from@dim, from@type, from@svtree, PACKAGE="S4Arrays")
    ans_nzindex <- C_ans[[1L]]
    ans_nzdata  <- C_ans[[2L]]
    new2("COOSparseArray", dim=from@dim, dimnames=from@dimnames,
                     nzindex=ans_nzindex, nzdata=ans_nzdata, check=FALSE)
}

.from_COOSparseArray_to_SVTSparseArray <- function(from)
{
    ans_svtree <- .Call2("C_from_COOSparseArray_to_SVTSparseArray",
                         from@dim, from@nzindex, from@nzdata,
                         PACKAGE="S4Arrays")
    new2("SVTSparseArray", dim=from@dim, dimnames=from@dimnames,
                           type=type(from), svtree=ans_svtree, check=FALSE)
}

setAs("SVTSparseArray", "COOSparseArray",
    .from_SVTSparseArray_to_COOSparseArray
)
setAs("COOSparseArray", "SVTSparseArray",
    .from_COOSparseArray_to_SVTSparseArray
)

