### =========================================================================
### The rowsum() and colsum() S4 generics
### -------------------------------------------------------------------------
###


### The base package provides rowsum() only (as an S3 generic).
### TODO: Maybe move these generics to the BiocGenerics or MatrixGenerics
### package?
setGeneric("rowsum", signature="x")

setGeneric("colsum", signature="x",
    function(x, group, reorder=TRUE, ...)
        standardGeneric("colsum")
)

setMethod("colsum", "ANY",
    function(x, group, reorder=TRUE, ...)
    {
        t(rowsum(t(x), group, reorder=reorder, ...))
    }
)

### NOT exported but used in the SparseArray and DelayedArray packages!
compute_ugroup <- function(group, expected_group_len, reorder=TRUE)
{
    if (!(is.vector(group) || is.factor(group)))
        stop(wmsg("'group' must be a vector or factor"))
    if (length(group) != expected_group_len)
        stop(wmsg("incorrect length for 'group'"))
    if (!isTRUEorFALSE(reorder))
        stop(wmsg("'reorder' must be TRUE or FALSE"))
    ## Taken from base::rowsum.default().
    ugroup <- unique(group)
    if (anyNA(ugroup))
        warning(wmsg("missing values for 'group'"))
    if (reorder)
        ugroup <- sort(ugroup, na.last=TRUE, method="quick")
    ugroup
}

