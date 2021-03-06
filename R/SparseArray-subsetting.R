### =========================================================================
### SparseArray subsetting
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The extract_sparse_array() generic
###
### extract_sparse_array() is the workhorse behind read_sparse_block(), and,
### more generally, behind efficient subsetting of sparse array objects.
### Similar to extract_array() except that:
###   (1) The extracted array data must be returned as a SparseArray object.
###       Methods should always operate on the sparse representation
###       of the data and never "expand" it, that is, never turn it into a
###       dense representation for example by doing something like
###       'dense2sparse(extract_array(x, index))'. This would defeat the
###       purpose of read_sparse_block().
###   (2) It should be called only on an array-like object 'x' for which
###       'is_sparse(x)' is TRUE.
###   (3) The subscripts in 'index' should NOT contain duplicates.
### IMPORTANT NOTE: For the sake of efficiency, (2) and (3) are NOT checked
### and are the responsibility of the user. We'll refer to (2) and (3) as
### the "extract_sparse_array() Terms of Use".

setGeneric("extract_sparse_array",
    function(x, index)
    {
        x_dim <- dim(x)
        if (is.null(x_dim))
            stop(wmsg("first argument to extract_sparse_array() ",
                      "must be an array-like object"))
        ans <- standardGeneric("extract_sparse_array")
        expected_dim <- get_Nindex_lengths(index, x_dim)
        ## TODO: Display a more user/developper-friendly error by
        ## doing something like the extract_array() generic where
        ## check_returned_array() is used to display a long and
        ## detailed error message.
        stopifnot(is(ans, "SparseArray"),
                  identical(dim(ans), expected_dim))
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_sparse_array() and extract_array() methods for COO_SparseArray
### objects
###

### IMPORTANT NOTE: The returned COO_SparseArray object is guaranteed to be
### **correct** ONLY if the subscripts in 'index' do NOT contain duplicates!
### If they contain duplicates, the correct COO_SparseArray object to return
### should contain repeated nonzero data. However, in order to keep it as
### efficient as possible, the code below does NOT repeat the nonzero data
### that corresponds to duplicates subscripts. It does not check for
### duplicates in 'index' either because this check could have a
### non-negligible cost.
### All this is OK because .extract_COO_SparseArray_subset() should
### always be used in a context where 'index' does NOT contain duplicates.
### The only situation where 'index' CAN contain duplicates is when
### .extract_COO_SparseArray_subset() is called by
### .extract_array_from_COO_SparseArray(), in which case the
### missing nonzero data are added later.
.extract_COO_SparseArray_subset <- function(x, index)
{
    stopifnot(is(x, "COO_SparseArray"))
    ans_dim <- get_Nindex_lengths(index, dim(x))
    x_nzcoo <- x@nzcoo
    for (along in seq_along(ans_dim)) {
        i <- index[[along]]
        if (is.null(i))
            next
        x_nzcoo[ , along] <- match(x_nzcoo[ , along], i)
    }
    keep_idx <- which(!rowAnyNAs(x_nzcoo))
    ans_nzcoo <- x_nzcoo[keep_idx, , drop=FALSE]
    ans_nzvals <- x@nzvals[keep_idx]
    COO_SparseArray(ans_dim, ans_nzcoo, ans_nzvals, check=FALSE)
}
setMethod("extract_sparse_array", "COO_SparseArray",
    .extract_COO_SparseArray_subset
)

.extract_array_from_COO_SparseArray <- function(x, index)
{
    coo0 <- .extract_COO_SparseArray_subset(x, index)
    ## If the subscripts in 'index' contain duplicates, 'coo0' is
    ## "incomplete" in the sense that it does not contain the nonzero data
    ## that should have been repeated according to the duplicates in the
    ## subscripts (see IMPORTANT NOTE above).
    ans0 <- sparse2dense(coo0)
    ## We "complete" 'ans0' by repeating the nonzero data according to the
    ## duplicates present in 'index'. Note that this is easy and cheap to
    ## do now because 'ans0' uses a dense representation (it's an ordinary
    ## array). This would be harder to do **natively** on the
    ## COO_SparseArray form (i.e. without converting to dense first
    ## then back to sparse).
    sm_index <- lapply(index,
        function(i) {
            if (is.null(i))
                return(NULL)
            sm <- match(i, i)
            if (isSequence(sm))
                return(NULL)
            sm
        })
    if (all(S4Vectors:::sapply_isNULL(sm_index)))
        return(ans0)
    subset_by_Nindex(ans0, sm_index)
}
setMethod("extract_array", "COO_SparseArray",
    .extract_array_from_COO_SparseArray
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_sparse_array() method for SVT_SparseArray objects
###

### This is the workhorse behind the extract_sparse_array(), extract_array(),
### and `[` methods for SVT_SparseArray objects.
### Like for 'extract_array()', the supplied 'index' must be a list with
### one list element per dimension in 'x'. Each list element must be an
### integer vector of valid indices along the corresponding dimension
### in 'x', or a NULL.
.subset_SVT_SparseArray <- function(x, index, ignore.dimnames=FALSE)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.list(index),
              length(index) == length(x@dim),
              isTRUEorFALSE(ignore.dimnames))

    ## Returns 'new_dim' and 'new_SVT' in a list of length 2.
    C_ans <- .Call2("C_subset_SVT_SparseArray",
                    x@dim, x@type, x@SVT, index, PACKAGE="S4Arrays")
    new_dim <- C_ans[[1L]]
    new_SVT <- C_ans[[2L]]

    ## Compute 'new_dimnames'.
    if (is.null(dimnames(x)) || ignore.dimnames) {
        new_dimnames <- vector("list", length(x@dim))
    } else {
        new_dimnames <- subset_dimnames_by_Nindex(x@dimnames, index)
    }
    BiocGenerics:::replaceSlots(x, dim=new_dim,
                                   dimnames=new_dimnames,
                                   SVT=new_SVT,
                                   check=FALSE)
}

### No need to propagate the dimnames.
setMethod("extract_sparse_array", "SVT_SparseArray",
    function(x, index) .subset_SVT_SparseArray(x, index, ignore.dimnames=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop()
###

### Always returns an SVT_SparseArray object (endomorphism).
.drop_SVT_SparseArray <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    ## Returns 'ans_dim', 'ans_dimnames', and 'ans_SVT', in a list of length 3.
    C_ans <- .Call2("C_drop_SVT_SparseArray_ineffective_dims",
                    x@dim, x@dimnames, x@type, x@SVT, PACKAGE="S4Arrays")
    ans_dim <- C_ans[[1L]]
    ans_dimnames <- C_ans[[2L]]
    ans_SVT <- C_ans[[3L]]
    new_SVT_SparseArray(ans_dim, ans_dimnames, x@type, ans_SVT, check=FALSE)
}

### Returns an SVT_SparseArray object or an ordinary vector of type 'type(x)'.
### The dimnames are propagated except when 'length(x)' is 1 (i.e.
### 'all(dim(x) == 1)' is TRUE).
setMethod("drop", "SVT_SparseArray",
    function(x)
    {
        x <- .drop_SVT_SparseArray(x)
        if (length(dim(x)) != 1L)
            return(x)     # SVT_SparseArray object
        a <- as.array(x)  # 1d ordinary array
        ans <- as.vector(a)  # unfortunately, this drops the names
        ## Restore the names.
        a_dimnames <- dimnames(a)  # NULL or list of length 1
        if (!is.null(a_dimnames))
            names(ans) <- a_dimnames[[1L]]
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SVT_SparseArray single-bracket subsetting (based on extract_sparse_array())
###

.subset_SVT_SparseArray_by_logical_array <- function(x, y, drop=TRUE)
    stop("subsetting operation not supported yet")

.subset_SVT_SparseArray_by_Mindex <- function(x, Mindex, drop=FALSE)
    stop("subsetting operation not supported yet")

.subset_SVT_SparseArray_by_Lindex <- function(x, Lindex)
    stop("subsetting operation not supported yet")

.single_bracket_SVT_SparseArray <- function(x, i, j, ..., drop=TRUE)
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
    x_ndim <- length(x_dim)
    if (nsubscript == 1L) {
        i <- Nindex[[1L]]
        if (type(i) == "logical" && identical(x_dim, dim(i)))
            return(.subset_SVT_SparseArray_by_logical_array(x, i, drop=drop))
        if (is.matrix(i) && is.numeric(i))
            return(.subset_SVT_SparseArray_by_Mindex(x, i, drop=drop))
        ## Linear single bracket subsetting e.g. x[5:2].
        ## If 'x' is monodimensional and 'drop' is FALSE, we fallback
        ## to "multi-dimensional single bracket subsetting" which is an
        ## endomorphism.
        if (x_ndim != 1L || drop)
            return(.subset_SVT_SparseArray_by_Lindex(x, i))
    }
    if (nsubscript != x_ndim)
        stop(wmsg("incorrect number of subscripts"))
    index <- normalize_Nindex(Nindex, x)
    ans <- .subset_SVT_SparseArray(x, index)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "SVT_SparseArray", .single_bracket_SVT_SparseArray)

### Note that the default extract_array() method would do the job but it
### relies on single-bracket subsetting so would needlessly go thru the
### complex .single_bracket_SVT_SparseArray() machinery above. It would
### also propagate the dimnames which extract_array() does not need to do.
### The method below completely bypasses all this complexity.
setMethod("extract_array", "SVT_SparseArray",
    function(x, index)
        as.array(.subset_SVT_SparseArray(x, index, ignore.dimnames=TRUE))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse() and extract_sparse_array() methods for [d|l]g[C|R]Matrix
### objects from the Matrix package
###
### TODO: Support more sparseMatrix derivatives (e.g. dgTMatrix, dgRMatrix)
### as the need arises.
###

setMethod("is_sparse", "dgCMatrix", function(x) TRUE)
setMethod("is_sparse", "lgCMatrix", function(x) TRUE)
setMethod("is_sparse", "dgRMatrix", function(x) TRUE)
setMethod("is_sparse", "lgRMatrix", function(x) TRUE)

### TODO: Return an SVT_SparseArray instead of a COO_SparseArray object when
### extracting from a dgCMatrix or lgCMatrix.
.extract_sparse_array_from_dgCMatrix_or_lgCMatrix <- function(x, index)
{
    sm <- subset_by_Nindex(x, index)  # a dgCMatrix or lgCMatrix object
    make_COO_SparseMatrix_from_dgCMatrix_or_lgCMatrix(sm, use.dimnames=FALSE)
}

.extract_sparse_array_from_dgRMatrix_or_lgRMatrix <- function(x, index)
{
    sm <- subset_by_Nindex(x, index)  # a dgRMatrix or lgRMatrix object
    make_COO_SparseMatrix_from_dgRMatrix_or_lgRMatrix(sm, use.dimnames=FALSE)
}

setMethod("extract_sparse_array", "dgCMatrix",
    .extract_sparse_array_from_dgCMatrix_or_lgCMatrix
)
setMethod("extract_sparse_array", "lgCMatrix",
    .extract_sparse_array_from_dgCMatrix_or_lgCMatrix
)
setMethod("extract_sparse_array", "dgRMatrix",
    .extract_sparse_array_from_dgRMatrix_or_lgRMatrix
)
setMethod("extract_sparse_array", "lgRMatrix",
    .extract_sparse_array_from_dgRMatrix_or_lgRMatrix
)

