### =========================================================================
### SparseArray objects
### -------------------------------------------------------------------------


setClass("SparseArray",
    contains="Array",
    representation(
        "VIRTUAL",
        dim="integer",     # This gives us dim() for free!
        dimnames="list"    # List with one list element per dimension. Each
                           # list element must be NULL or a character vector.
    ),
    prototype(
        dim=0L,
        dimnames=list(NULL)
    )
)

.validate_SparseArray <- function(x)
{
    msg <- validate_dim_slot(x, "dim")
    if (!isTRUE(msg))
        return(msg)
    msg <- validate_dimnames_slot(x, x@dim)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("SparseArray", .validate_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseArray API:
### - Getters: dim(), length(), dimnames(), type().
### - Setters: `dimnames<-`(), `type<-`().
### - is_sparse().
###

setMethod("dimnames", "SparseArray",
    function(x) simplify_NULL_dimnames(x@dimnames)
)

setReplaceMethod("dimnames", "SparseArray",
    function(x, value)
    {
        x@dimnames <- normarg_dimnames(value, dim(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### which_is_nonzero() and coercion_can_introduce_zeros()
###

which_is_nonzero <- function(x, arr.ind=FALSE)
{
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    is_nonzero <- x != zero
    which(is_nonzero | is.na(is_nonzero), arr.ind=arr.ind)
}

coercion_can_introduce_zeros <- function(from_type, to_type)
{
    if (identical(to_type, "double"))
        return(from_type %in% c("logical", "integer", "raw"))
    if (identical(to_type, "logical"))
        return(from_type %in% c("integer", "double", "complex", "raw"))
    stop(wmsg("'to_type' must be \"double\" or \"logical\""))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics is_sparse(), is_sparse<-(), and sparsity()
###

### is_sparse() detects **structural** sparsity which is a global qualitative
### property of array-like object 'x' rather than a quantitative one.
### In other words it doesn't look at the data in 'x' to decide whether 'x'
### should be considered sparse or not. Said otherwise, it is NOT about
### quantitative sparsity measured by sparsity().
### IMPORTANT: Seeds for which is_sparse() returns TRUE **must** support
### extract_sparse_array().
setGeneric("is_sparse", function(x) standardGeneric("is_sparse"))

setGeneric("is_sparse<-", signature="x",
    function(x, value) standardGeneric("is_sparse<-")
)

### By default, nothing is considered sparse.
setMethod("is_sparse", "ANY", function(x) FALSE)

setMethod("is_sparse", "SparseArray", function(x) TRUE)

setGeneric("sparsity", function(x) standardGeneric("sparsity"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The extract_sparse_array() generics
###

### This is the workhorse behind read_sparse_block().
### Similar to extract_array() except that:
###   (1) The extracted array data must be returned in a COO_SparseArray
###       object. Methods should always operate on the sparse representation
###       of the data and never "expand" it, that is, never turn it into a
###       dense representation for example by doing something like
###       'dense2sparse(extract_array(x, index))'. This would defeat the
###       purpose of read_sparse_block().
###   (2) It should be called only on an array-like object 'x' for which
###       'is_sparse(x)' is TRUE.
###   (3) The subscripts in 'index' should NOT contain duplicates.
### IMPORTANT NOTE: For the sake of efficiency, (2) and (3) are NOT checked
### and are the responsibility of the user. We'll refer to (2) and (3) as
### the "extract_sparse_array() Terms of Use".
setGeneric("extract_sparse_array",
    function(x, index)
    {
        x_dim <- dim(x)
        if (is.null(x_dim))
            stop(wmsg("first argument to extract_sparse_array() ",
                      "must be an array-like object"))
        ans <- standardGeneric("extract_sparse_array")
        expected_dim <- get_Nindex_lengths(index, x_dim)
        ## TODO: Display a more user/developper-friendly error by
        ## doing something like the extract_array() generic where
        ## check_returned_array() is used to display a long and
        ## detailed error message.
        stopifnot(is(ans, "COO_SparseArray"),
                  identical(dim(ans), expected_dim))
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "SparseArray",
    function(object) show_compact_array(object)
)

