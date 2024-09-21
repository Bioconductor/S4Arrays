### =========================================================================
### Array subsetting
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level generics to support subsetting of Array derivatives
###
### We define 4 low-level generics that are called by the "[" method for
### Array objects defined below in this file. The aim is to falicitate the
### implementation of subsetting operations on array-like S4 objects.
###

setGeneric("subset_Array_by_logical_array", signature="x",
    function(x, y) standardGeneric("subset_Array_by_logical_array")
)
setGeneric("subset_Array_by_Lindex", signature="x",
    function(x, Lindex) standardGeneric("subset_Array_by_Lindex")
)
setGeneric("subset_Array_by_Mindex", signature="x",
    function(x, Mindex) standardGeneric("subset_Array_by_Mindex")
)
setGeneric("subset_Array_by_Nindex", signature="x",
    function(x, Nindex) standardGeneric("subset_Array_by_Nindex")
)

setMethod("subset_Array_by_logical_array", "Array",
    function(x, y)
        stop(wmsg(class(x)[[1L]], " objects don't support this ",
                  "form of subsetting at the moment"))
)

setMethod("subset_Array_by_Lindex", "Array",
    function(x, Lindex)
        stop(wmsg(class(x)[[1L]], " objects don't support this ",
                  "form of subsetting at the moment"))
)

### Simply delegates to subset_Array_by_Lindex().
.subset_Array_by_Mindex <- function(x, Mindex)
{
    stopifnot(is.matrix(Mindex), is.numeric(Mindex))
    subset_Array_by_Lindex(x, Mindex2Lindex(Mindex, dim(x)))
}

setMethod("subset_Array_by_Mindex", "Array", .subset_Array_by_Mindex)

setMethod("subset_Array_by_Nindex", "Array",
    function(x, Nindex)
        stop(wmsg(class(x)[[1L]], " objects don't support this ",
                  "form of subsetting at the moment"))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "[" method for Array objects
###

### Works on any array-like object that supports the 4 low-level generics
### above.
.subset_Array <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop(wmsg("'x' is missing"))
    if (!isTRUEorFALSE(drop))
        stop(wmsg("'drop' must be TRUE or FALSE"))
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    if (nsubscript == 0L)
        return(x)  # no-op
    x_dim <- dim(x)
    if (nsubscript == 1L && drop) {
        i <- Nindex[[1L]]
        if (type(i) == "logical" && identical(x_dim, dim(i)))
            return(subset_Array_by_logical_array(x, i))
        if (is.matrix(i))
            return(subset_Array_by_Mindex(x, i))
        if (is.numeric(i))
            return(subset_Array_by_Lindex(x, i))
    }
    if (nsubscript != length(x_dim))
        stop(wmsg("incorrect number of subscripts"))
    ## Note that this normalization will coerce the numeric subscripts
    ## in 'Nindex' to integer. However this coercion is no longer necessary
    ## because subset_Array_by_Nindex() should be able to handle subscripts
    ## of type "double".
    ## TODO: Consider using a normalization process here that preserves the
    ## numeric subscripts.
    Nindex <- normalize_Nindex(Nindex, x)
    ans <- subset_Array_by_Nindex(x, Nindex)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "Array", .subset_Array)

