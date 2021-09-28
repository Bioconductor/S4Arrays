### =========================================================================
### SparseArray subassignment
### -------------------------------------------------------------------------
###


.subassign_SVT_SparseArray_by_logical_array <- function(x, y, value)
    stop("subassignment operation not supported yet")

.subassign_SVT_SparseArray_by_Mindex <- function(x, Mindex, value)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.matrix(Mindex), is.numeric(Mindex))
    if (!is.vector(value))
        stop(wmsg("the supplied value must be a vector for this form ",
                  "of subassignment to an SVT_SparseArray object"))

    ## Change 'x' type if necessary. */
    new_type <- type(c(vector(type(x)), vector(type(value))))
    type(x) <- new_type
    if (nrow(Mindex) == 0L)
        return(x)  # nothing else to do

    ## Normalize 'value'. */
    storage.mode(value) <- new_type
    if (length(value) == 0L)
        stop(wmsg("replacement has length zero"))
    value <- S4Vectors:::recycleVector(value, nrow(Mindex))

    if (storage.mode(Mindex) != "integer")
        storage.mode(Mindex) <- "integer"
    new_SVT <- .Call2("C_subassign_SVT_by_Mindex",
                      x@dim, x@type, x@SVT, Mindex, value, PACKAGE="S4Arrays")
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

.subassign_SVT_SparseArray_by_Lindex <- function(x, Lindex, value)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.vector(Lindex), is.numeric(Lindex))
    if (!is.vector(value))
        stop(wmsg("the supplied value must be a vector for this form ",
                  "of subassignment to an SVT_SparseArray object"))

    ## Change 'x' type if necessary. */
    new_type <- type(c(vector(type(x)), vector(type(value))))
    type(x) <- new_type
    if (length(Lindex) == 0L)
        return(x)  # nothing else to do

    ## Normalize 'value'. */
    storage.mode(value) <- new_type
    if (length(value) == 0L)
        stop(wmsg("replacement has length zero"))
    value <- S4Vectors:::recycleVector(value, length(Lindex))

    new_SVT <- .Call2("C_subassign_SVT_by_Lindex",
                      x@dim, x@type, x@SVT, Lindex, value, PACKAGE="S4Arrays")
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

.subassign_SVT_SparseArray <- function(x, i, j, ..., value)
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
            return(.subassign_SVT_SparseArray_by_logical_array(x, i, value))
        if (is.matrix(i) && is.numeric(i))
            return(.subassign_SVT_SparseArray_by_Mindex(x, i, value))
        ## Linear single bracket subassignment e.g. x[5:2] <- 4.
        return(.subassign_SVT_SparseArray_by_Lindex(x, i, value))
    }
    if (nsubscript != x_ndim)
        stop(wmsg("incorrect number of subscripts"))
    index <- normalize_Nindex(Nindex, x)
    .subassign_SVT_SparseArray_by_Nindex(x, index, value)
}

setReplaceMethod("[", "SVT_SparseArray", .subassign_SVT_SparseArray)

