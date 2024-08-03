### =========================================================================
### Array subassignment
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level generics to support subassignment of Array derivatives
###
### We define 4 low-level generics that are called by the "[<-" method for
### Array objects defined below in this file. The aim is to falicitate the
### implementation of subassignment operations on array-like S4 objects.
###

setGeneric("subassign_Array_by_logical_array", signature="x",
    function(x, y, value) standardGeneric("subassign_Array_by_logical_array")
)
setGeneric("subassign_Array_by_Lindex", signature="x",
    function(x, Lindex, value) standardGeneric("subassign_Array_by_Lindex")
)
setGeneric("subassign_Array_by_Mindex", signature="x",
    function(x, Mindex, value) standardGeneric("subassign_Array_by_Mindex")
)
setGeneric("subassign_Array_by_Nindex", signature="x",
    function(x, Nindex, value) standardGeneric("subassign_Array_by_Nindex")
)

setMethod("subassign_Array_by_logical_array", "Array",
    function(x, y, value)
        stop(wmsg(class(x)[[1L]], " objects don't support this ",
                  "form of subassignment at the moment"))
)

setMethod("subassign_Array_by_Lindex", "Array",
    function(x, Lindex, value)
        stop(wmsg(class(x)[[1L]], " objects don't support this ",
                  "form of subassignment at the moment"))
)

### Simply delegates to subassign_Array_by_Lindex().
.subassign_Array_by_Mindex <- function(x, Mindex, value)
{
    stopifnot(is.matrix(Mindex), is.numeric(Mindex))
    subassign_Array_by_Lindex(x, Mindex2Lindex(Mindex, dim(x)), value)
}

setMethod("subassign_Array_by_Mindex", "Array", .subassign_Array_by_Mindex)

setMethod("subassign_Array_by_Nindex", "Array",
    function(x, Nindex, value)
        stop(wmsg(class(x)[[1L]], " objects don't support this ",
                  "form of subassignment at the moment"))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "[<-" method for Array objects
###

### Works on any array-like object that supports the 4 low-level generics
### above.
.subassign_Array <- function(x, i, j, ..., value)
{
    if (missing(x))
        stop(wmsg("'x' is missing"))
    Nindex <- extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (nsubscript == 1L) {
        i <- Nindex[[1L]]
        if (type(i) == "logical" && identical(x_dim, dim(i)))
            return(subassign_Array_by_logical_array(x, i, value))
        if (is.matrix(i) && is.numeric(i))
            return(subassign_Array_by_Mindex(x, i, value))
        ## Linear single bracket subassignment e.g. x[5:2] <- 4.
        return(subassign_Array_by_Lindex(x, i, value))
    }
    if (nsubscript != x_ndim)
        stop(wmsg("incorrect number of subscripts"))
    Nindex <- normalize_Nindex(Nindex, x)
    subassign_Array_by_Nindex(x, Nindex, value)
}

setReplaceMethod("[", "Array", .subassign_Array)

