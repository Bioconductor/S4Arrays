### =========================================================================
### The is_sparse() and `is_sparse<-`() generics
### -------------------------------------------------------------------------
###

### IMPORTANT: Array-like objects for which is_sparse() is TRUE **must**
### support SparseArray::extract_sparse_array(). More about this in the
### SparseArray-subsetting.R file from the SparseArray package where the
### extract_sparse_array() generic is defined.

setGeneric("is_sparse", function(x) standardGeneric("is_sparse"))

### By default, nothing is considered sparse.
setMethod("is_sparse", "ANY", function(x) FALSE)

### The corresponding (and **mandatory**, see IMPORTANT note above)
### SparseArray::extract_sparse_array() methods are defined in the
### SparseArray package.
setMethod("is_sparse", "dgCMatrix", function(x) TRUE)
setMethod("is_sparse", "lgCMatrix", function(x) TRUE)
setMethod("is_sparse", "dgRMatrix", function(x) TRUE)
setMethod("is_sparse", "lgRMatrix", function(x) TRUE)

setGeneric("is_sparse<-", signature="x",
    function(x, value) standardGeneric("is_sparse<-")
)

