### =========================================================================
### ArraySelection objects
### -------------------------------------------------------------------------


setClass("ArraySelection",
    representation(
        "VIRTUAL",
        # Dimensions of "the reference array" i.e. the array on top of which
        # the selection is defined.
        refdim="integer",
        selection="ANY"
    ),
    prototype(refdim=0L)
)

.validate_ArraySelection <- function(x)
{
    msg <- validate_dim_slot(x, "refdim")
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("ArraySelection", .validate_ArraySelection)

### IMPORTANT NOTES:
### - Because ordinary matrices in base R must have less than 2^31 rows, an
###   Mindex object cannot represent a "long selection" i.e. a selection of
###   2^31 array elements or more!
### - Mindex not a good name as it can easily be confused with MIndex from
###   Biostrings!
setClass("Mindex",
    contains="ArraySelection",
    representation(
        # An integer matrix with one column per dimension in the reference
        # array and one row per selected element in the reference array.
        # Only non-NA positive integers are allowed in the matrix.
        selection="matrix"
    )
)

### A workaround for the limitation of Mindex objects (see above).
setClass("LongMindex",
    contains="ArraySelection",
    representation(
        # A list of integer vectors of the same length. One integer vector
        # per dimension in the reference array. Only non-NA positive integers
        # are allowed in the individual vectors.
        selection="list"
    )
)

setClassUnion("NULL_OR_integer_OR_list", c("NULL", "integer", "list"))

### Does not have the limitation of Mindex objects (see above) and much
### compact than Mindex or LongMindex.
setClass("SelectionTree",
    contains="ArraySelection",
    representation(
        # - A tree of depth 1 if the reference array has 1 dimension. A tree
        #   of depth 1 is represented by a NULL or an integer vector.
        #
        # - A tree of depth 2 if the reference array has 2 dimensions. A tree
        #   of depth 2 is represented by a NULL or a list with one list element
        #   per column. If the latter, each list element must be a "tree of
        #   depth 1" i.e. a NULL or an integer vector of row indices.
        #
        # - A tree of depth 3 if the reference array has 3 dimensions. A tree
        #   of depth 3 is represented by a NULL or a list of length the extend
        #   of the 3rd dimension. If the latter, each list element must be
        #   a "tree of depth 2" (like described above) that represents a
        #   selection on the corresponding 2D slice.
        #
        # - etc...
        selection="NULL_OR_integer_OR_list"
    )
)

### Equivalent to 'length(unlist(x@selection, recursive=TRUE))' but much
### faster and more memory efficient.
### Note that like for atomic vectors, the returned length will be a double
### if it's > .Machine$integer.max
.get_SelectionTree_length <- function(x)
    .Call2("C_get_SelectionTree_length", x@refdim, x@selection,
                                         PACKAGE="S4Arrays")

setMethod("length", "SelectionTree", .get_SelectionTree_length)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SelectionTree and matrix objects
###

.from_SelectionTree_to_matrix <- function(x)
    .Call2("C_from_SelectionTree_to_matrix", x@refdim, x@selection,
                                             PACKAGE="S4Arrays")

.from_matrix_to_SelectionTree <- function(m, refdim)
{
    stopifnot(is.matrix(m), ncol(m) >= 1L, is.integer(m))
    refdim <- normarg_dim(refdim, "refdim")
    stopifnot(length(refdim) == ncol(m))

    ans_selection <- .Call2("C_from_matrix_to_SelectionTree", m, refdim,
                            PACKAGE="S4Arrays")
    new2("SelectionTree", refdim=refdim, selection=ans_selection)
}

### S3/S4 combo for as.matrix.SelectionTree
as.matrix.SelectionTree <-
    function(x, ...) .from_SelectionTree_to_matrix(x, ...)
setMethod("as.matrix", "SelectionTree", .from_SelectionTree_to_matrix)

### SelectionTree constructor.
SelectionTree <- function(refdim, m) .from_matrix_to_SelectionTree(m, refdim)

