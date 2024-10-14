### =========================================================================
### The is_sparse() and `is_sparse<-`() generics
### -------------------------------------------------------------------------
###

### IMPORTANT: Array-like objects for which is_sparse() is TRUE **must**
### support SparseArray::extract_sparse_array(). More about this in the
### extract_sparse_array.R file from the SparseArray package where the
### extract_sparse_array() generic is defined.

setGeneric("is_sparse", function(x) standardGeneric("is_sparse"))

### By default, nothing is considered sparse.
setMethod("is_sparse", "ANY", function(x) FALSE)

### All [C|R|T]sparseMatrix derivatives are considered sparse. Note that
### they all support SparseArray::extract_sparse_array() via the default
### extract_sparse_array() method defined in the SparseArray package.
setMethod("is_sparse", "CsparseMatrix", function(x) TRUE)
setMethod("is_sparse", "RsparseMatrix", function(x) TRUE)
setMethod("is_sparse", "TsparseMatrix", function(x) TRUE)


setGeneric("is_sparse<-", signature="x",
    function(x, value) standardGeneric("is_sparse<-")
)

