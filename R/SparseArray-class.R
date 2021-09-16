### =========================================================================
### SparseArray objects
### -------------------------------------------------------------------------


### Virtual class with 2 concrete subclasses: COO_SparseArray and
### SVT_SparseArray.
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
### extract_sparse_array(). More info about this in extract_sparse_array.R
setGeneric("is_sparse", function(x) standardGeneric("is_sparse"))

setGeneric("is_sparse<-", signature="x",
    function(x, value) standardGeneric("is_sparse<-")
)

### By default, nothing is considered sparse.
setMethod("is_sparse", "ANY", function(x) FALSE)

setMethod("is_sparse", "SparseArray", function(x) TRUE)

setGeneric("sparsity", function(x) standardGeneric("sparsity"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

### show_compact_array() relies on extract_array() so show() will work
### out-of-the-box on any SparseArray derivative that supports aperm().
setMethod("show", "SparseArray",
    function(object) show_compact_array(object)
)

