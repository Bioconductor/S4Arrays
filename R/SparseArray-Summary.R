### =========================================================================
### "Summary" methods for SparseArray objects
### -------------------------------------------------------------------------
###
### "Summary" is a group generic with members: max(), min(), range(), sum(),
### prod(), any(), all(). In this file we defines methods to make all these
### functions work on SparseArray objects.
###
### We also define methods to make mean() and anyNA() work on these objects.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### max(), min(), range(), sum(), prod(), any(), all()
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

### S3/S4 combo for range.SVT_SparseArray
range.SVT_SparseArray <- function(..., na.rm=FALSE, finite=FALSE)
{
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("range() method for SVT_SparseArray objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    stop(wmsg("not ready yet!"))
    #x_has_zeros <- length(x@nzvals) < length(x)
    #if (!x_has_zeros)
    #    return(range(x@nzvals, na.rm=na.rm, finite=finite))
    #zero <- vector(typeof(x@nzvals), length=1L)
    #range(zero, x@nzvals, na.rm=na.rm, finite=finite)
}
### The signature of all the members of the S4 "Summary" group generic is
### 'x, ..., na.rm' (see getGeneric("range")) which means that the S4 methods
### cannot add arguments after 'na.rm'. So we add the 'finite' argument before.
setMethod("range", "SVT_SparseArray",
    function(x, ..., finite=FALSE, na.rm=FALSE)
        range.SVT_SparseArray(x, ..., na.rm=na.rm, finite=finite)

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### TODO: Maybe introduce a new generic for this e.g. countNAs()?
.count_SparseArray_NAs <- function(x)
{
    if (is(x, "COO_SparseArray"))
        return(sum(is.na(x@nzvals)))

    if (is(x, "SVT_SparseArray"))
        return(.Call2("C_count_SVT_SparseArray_NAs",
                      x@dim, x@type, x@SVT, PACKAGE="S4Arrays"))

    stop(wmsg(class(x)[[1L]], " objects are not supported"))
}

.mean_SparseArray <- function(x, na.rm=FALSE)
{
    s <- as.double(sum(x, na.rm=na.rm))
    nval <- length(x)
    if (na.rm)
        nval <- nval - .count_SparseArray_NAs(x)
    s / nval
}

### S3/S4 combo for mean.SparseArray
mean.SparseArray <- function(x, na.rm=FALSE, ...)
    .mean_SparseArray(x, na.rm=na.rm, ...)
setMethod("mean", "SparseArray", .mean_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

setMethod("anyNA", "COO_SparseArray",
    function(x, recursive=FALSE) anyNA(x@nzvals, recursive=recursive)
)

setMethod("anyNA", "SVT_SparseArray",
    function(x, recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop(wmsg("anyNA() method for SVT_SparseArray objects ",
                      "does not support the 'recursive' argument"))
        .Call2("C_anyNA_SVT_SparseArray",
               x@dim, x@type, x@SVT, PACKAGE="S4Arrays")
    }
)

