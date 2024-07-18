### =========================================================================
### The rowsum() and colsum() S4 generics
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .fast_rowsum() and .fast_colsum()
###

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

### A fast re-implementation of base::rowsum().
.fast_rowsum <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is.matrix(x))
    ugroup <- compute_ugroup(group, nrow(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    ans <- .Call2("C_rowsum", x, group, length(ugroup), na.rm,
                              PACKAGE="S4Arrays")
    set_dimnames(ans, list(as.character(ugroup), colnames(x)))
}

.fast_colsum <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is.matrix(x))
    ugroup <- compute_ugroup(group, ncol(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    ans <- .Call2("C_colsum", x, group, length(ugroup), na.rm,
                              PACKAGE="S4Arrays")
    set_dimnames(ans, list(rownames(x), as.character(ugroup)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The rowsum() and colsum() S4 generics and methods
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
        t(rowsum(t(x), group, reorder=reorder, ...))
)

setMethod("colsum", "matrix",
    function(x, group, reorder=TRUE, ...)
        .fast_colsum(x, group, reorder=reorder, ...)
)

