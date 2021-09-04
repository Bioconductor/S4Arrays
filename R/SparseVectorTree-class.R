### =========================================================================
### SparseVectorTree objects
### -------------------------------------------------------------------------


setClassUnion("NULL_OR_list", c("NULL", "list"))

setClass("SparseVectorTree",
    contains="SparseArray",
    representation(
        type="character",

        # - A tree of depth 1 if the reference array has 1 dimension. A tree
        #   of depth 1 is represented by a NULL or a "leaf vector".
        #   A leaf vector is a sparse vector represented by a list of 2
        #   parallel vectors: an integer vector of positions and a vector
        #   (atomic or list, exact type should match 'type' slot) of nonzero
        #   values.
        #
        # - A tree of depth 2 if the reference array has 2 dimensions. A tree
        #   of depth 2 is represented by a NULL or a list with one list element
        #   per column. If the latter, each list element must be a "tree of
        #   depth 1" i.e. a NULL or a "leaf vector".
        #
        # - A tree of depth 3 if the reference array has 3 dimensions. A tree
        #   of depth 3 is represented by a NULL or a list of length the extend
        #   of the 3rd dimension. If the latter, each list element must be
        #   a "tree of depth 2" (like described above) that represents a
        #   selection on the corresponding 2D slice.
        #
        # - etc...
        tree="NULL_OR_list"
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

setMethod("type", "SparseVectorTree", function(x) x@type)

### Note that like for the length of atomic vectors in base R, the returned
### length will be a double if it's > .Machine$integer.max
.get_SparseVectorTree_nzdata_length <- function(x)
    .Call2("C_get_SparseVectorTree_nzdata_length", x@dim, x@tree,
                                                   PACKAGE="S4Arrays")

#setMethod("nzdataLength", "SparseVectorTree",
#    .get_SparseVectorTree_nzdata_length
#)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SparseVectorTree and COOSparseArray objects
###

.from_SparseVectorTree_to_COOSparseArray <- function(from)
{
    ## Returns 'ans_nzindex' and 'ans_nzdata' in a list of length 2.
    C_ans <- .Call2("C_from_SparseVectorTree_to_COOSparseArray",
                    from@dim, from@type, from@tree, PACKAGE="S4Arrays")
    ans_nzindex <- C_ans[[1L]]
    ans_nzdata  <- C_ans[[2L]]
    new2("COOSparseArray", dim=from@dim, dimnames=from@dimnames,
                     nzindex=ans_nzindex, nzdata=ans_nzdata, check=FALSE)
}

.from_COOSparseArray_to_SparseVectorTree <- function(from)
{
    ans_tree <- .Call2("C_from_COOSparseArray_to_SparseVectorTree",
                       from@dim, from@nzindex, from@nzdata, PACKAGE="S4Arrays")
    new2("SparseVectorTree", dim=from@dim, dimnames=from@dimnames,
                             type=type(from), tree=ans_tree, check=FALSE)
}

setAs("SparseVectorTree", "COOSparseArray",
    .from_SparseVectorTree_to_COOSparseArray
)
setAs("COOSparseArray", "SparseVectorTree",
    .from_COOSparseArray_to_SparseVectorTree
)

