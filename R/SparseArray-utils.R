### =========================================================================
### Operate natively on SparseArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "unary isometric" array transformations
###
### A "unary isometric" array transformation is a transformation that returns
### an array-like object with the same dimensions as the input and where each
### element is the result of applying a function to the corresponding element
### in the input.
###
### Note that some "unary isometric" transformations preserve sparsity (e.g.
### is.na(), nchar(), round(), sqrt(), log1p(), etc...) and others don't
### (e.g. is.finite(), !, log(), etc..). We only implement the former.
###
### All the "unary isometric" array transformations implemented in this
### section return a COO_SparseArray object of the same dimensions as the
### input COO_SparseArray object.
###

.UNARY_ISO_OPS <- c("is.na", "is.infinite", "is.nan", "tolower", "toupper")

for (.Generic in .UNARY_ISO_OPS) {
    setMethod(.Generic, "COO_SparseArray",
        function(x)
        {
            GENERIC <- match.fun(.Generic)
            new_nzvals <- GENERIC(x@nzvals)
            BiocGenerics:::replaceSlots(x, nzvals=new_nzvals, check=FALSE)
        }
    )
}

setMethod("nchar", "COO_SparseArray",
    function(x, type="chars", allowNA=FALSE, keepNA=NA)
    {
        new_nzvals <- nchar(x@nzvals, type=type, allowNA=allowNA, keepNA=keepNA)
        BiocGenerics:::replaceSlots(x, nzvals=new_nzvals, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

setMethod("anyNA", "COO_SparseArray",
    function(x, recursive=FALSE) anyNA(x@nzvals, recursive=recursive)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### which()
###

.nzcoo_order <- function(nzcoo)
    do.call(order, lapply(ncol(nzcoo):1L, function(along) nzcoo[ , along]))

setMethod("which", "COO_SparseArray",
    function(x, arr.ind=FALSE, useNames=TRUE)
    {
        if (!identical(useNames, TRUE))
            warning(wmsg("'useNames' is ignored when 'x' is ",
                         "a COO_SparseArray object or derivative"))
        if (!isTRUEorFALSE(arr.ind))
            stop(wmsg("'arr.ind' must be TRUE or FALSE"))
        idx1 <- which(x@nzvals)
        nzcoo1 <- x@nzcoo[idx1, , drop=FALSE]
        oo <- .nzcoo_order(nzcoo1)
        ans <- nzcoo1[oo, , drop=FALSE]
        if (arr.ind)
            return(ans)
        Mindex2Lindex(ans, dim=dim(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Summary" group generic
###
### Members: max(), min(), range(), sum(), prod(), any(), all()
###

setMethod("Summary", "COO_SparseArray",
    function(x, ..., na.rm=FALSE)
    {
        GENERIC <- match.fun(.Generic)
        if (length(list(...)) != 0L)
            stop(wmsg(.Generic, "() method for COO_SparseArray objects ",
                      "only accepts a single object"))
        ## Whether 'x' contains zeros or not doesn't make a difference for
        ## sum() and any().
        if (.Generic %in% c("sum", "any"))
            return(GENERIC(x@nzvals, na.rm=na.rm))
        ## Of course a typical COO_SparseArray object "contains" zeros
        ## (i.e. it would contain zeros if we converted it to a dense
        ## representation with sparse2dense()). However, this is not
        ## guaranteed so we need to make sure to properly handle the case
        ## where it doesn't (admittedly unusual and definitely an inefficient
        ## way to represent dense data!)
        x_has_zeros <- length(x@nzvals) < length(x)
        if (!x_has_zeros)
            return(GENERIC(x@nzvals, na.rm=na.rm))
        x_type <- typeof(x@nzvals)
        if (.Generic == "all") {
            ## Mimic what 'all(sparse2dense(x))' would do.
            if (x_type == "double")
                warning("coercing argument of type 'double' to logical")
            return(FALSE)
        }
        zero <- vector(x_type, length=1L)
        GENERIC(zero, x@nzvals, na.rm=na.rm)
    }
)

### We override the range() method defined above via the Summary() method
### because we want to support the 'finite' argument like S3 method
### base::range.default() does.

### S3/S4 combo for range.COO_SparseArray
range.COO_SparseArray <- function(..., na.rm=FALSE, finite=FALSE)
{
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("range() method for COO_SparseArray objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    x_has_zeros <- length(x@nzvals) < length(x)
    if (!x_has_zeros)
        return(range(x@nzvals, na.rm=na.rm, finite=finite))
    zero <- vector(typeof(x@nzvals), length=1L)
    range(zero, x@nzvals, na.rm=na.rm, finite=finite)
}
### The signature of all the members of the S4 "Summary" group generic is
### 'x, ..., na.rm' (see getGeneric("range")) which means that the S4 methods
### cannot add arguments after 'na.rm'. So we add the 'finite' argument before.
setMethod("range", "COO_SparseArray",
    function(x, ..., finite=FALSE, na.rm=FALSE)
        range.COO_SparseArray(x, ..., na.rm=na.rm, finite=finite)

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

.mean_COO_SparseArray <- function(x, na.rm=FALSE)
{
    s <- sum(x@nzvals, na.rm=na.rm)
    nval <- length(x)
    if (na.rm)
        nval <- nval - sum(is.na(x@nzvals))
    s / nval
}

### S3/S4 combo for mean.COO_SparseArray
mean.COO_SparseArray <- function(x, na.rm=FALSE, ...)
    .mean_COO_SparseArray(x, na.rm=na.rm, ...)
setMethod("mean", "COO_SparseArray", .mean_COO_SparseArray)

