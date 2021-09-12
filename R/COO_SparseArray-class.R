### =========================================================================
### COO_SparseArray objects
### -------------------------------------------------------------------------
###
### SparseArray objects using the COO layout to store the nonzero data.
###
### Same as SparseArraySeed objects in the DelayedArray package.
### Extends the Coordinate List (COO) layout used for sparse matrices to
### multiple dimensions.
### See https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
### This layout is also used by https://sparse.pydata.org/
###
### COO_SparseArray API:
### - The SparseArray API.
### - Getters: nzcoo(), nzdata().
### - sparsity().
### - dense2sparse(), sparse2dense().
### - Based on sparse2dense(): extract_array(), as.array(), as.matrix().
### - Based on dense2sparse(): coercion to COO_SparseArray.
### - Back and forth coercion between COO_SparseArray and dg[C|R]Matrix or
###   lg[C|R]Matrix objects from the Matrix package.
###

setClass("COO_SparseArray",
    contains="SparseArray",
    representation(
        nzcoo="matrix",  # M-index containing the coordinates of the
                         # nonzero data.
        nzdata="vector"  # A vector (atomic or list) of length
                         # 'nrow(nzcoo)' containing the nonzero data.
    ),
    prototype(
        nzcoo=matrix(integer(0), ncol=1L),
        nzdata=logical(0)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_nzcoo_slot <- function(x)
{
    x_nzcoo <- x@nzcoo
    if (!(is.matrix(x_nzcoo) && typeof(x_nzcoo) == "integer"))
        return("'nzcoo' slot must be an integer matrix")
    x_dim <- x@dim
    if (ncol(x_nzcoo) != length(x_dim))
        return(paste0("'nzcoo' slot must be a matrix with ",
                      "one column per dimension"))
    for (along in seq_along(x_dim)) {
        not_ok <- S4Vectors:::anyMissingOrOutside(x_nzcoo[ , along],
                                                  1L, x_dim[[along]])
        if (not_ok)
            return(paste0("'nzcoo' slot must contain valid indices, ",
                          "that is, indices that are not NA and are ",
                          ">= 1 and <= their corresponding dimension"))
    }
    TRUE
}

.validate_nzdata_slot <- function(x)
{
    x_nzdata <- x@nzdata
    if (!(is.vector(x_nzdata) && length(x_nzdata) == nrow(x@nzcoo)))
        return(paste0("'nzdata' slot must be a vector of length ",
                      "the number of rows in the 'nzcoo' slot"))
    TRUE
}

.validate_COO_SparseArray <- function(x)
{
    msg <- .validate_nzcoo_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzdata_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("COO_SparseArray", .validate_COO_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("type", "COO_SparseArray", function(x) type(x@nzdata))

setGeneric("nzcoo", function(x) standardGeneric("nzcoo"))
setMethod("nzcoo", "COO_SparseArray", function(x) x@nzcoo)

setGeneric("nzdata", function(x) standardGeneric("nzdata"))
setMethod("nzdata", "COO_SparseArray", function(x) x@nzdata)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() setter
###

.set_COO_SparseArray_type <- function(x, value)
{
    stopifnot(is(x, "COO_SparseArray"))

    value <- normarg_array_type(value, "the supplied type")
    x_type <- type(x)
    if (value == x_type)
        return(x)

    new_nzdata <- x@nzdata
    storage.mode(new_nzdata) <- value
    zero <- vector(value, length=1L)
    is_not_zero <- new_nzdata != zero
    nzidx <- which(is_not_zero | is.na(is_not_zero))
    new_nzcoo <- x@nzcoo[nzidx, , drop=FALSE]
    new_nzdata <- new_nzdata[nzidx]

    BiocGenerics:::replaceSlots(x, nzcoo=new_nzcoo,
                                   nzdata=new_nzdata,
                                   check=FALSE)
}

setReplaceMethod("type", "COO_SparseArray", .set_COO_SparseArray_type)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sparsity()
###

setGeneric("sparsity", function(x) standardGeneric("sparsity"))

setMethod("sparsity", "COO_SparseArray",
    function(x) { 1 - length(nzdata(x)) / length(x) }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.normarg_nzdata <- function(nzdata, length.out)
{
    if (is.null(nzdata))
        stop(wmsg("'nzdata' cannot be NULL when 'nzcoo' is not NULL"))
    if (!is.vector(nzdata))
        stop(wmsg("'nzdata' must be a vector"))
    ## Same logic as S4Vectors:::V_recycle().
    nzdata_len <- length(nzdata)
    if (nzdata_len == length.out)
        return(nzdata)
    if (nzdata_len > length.out && nzdata_len != 1L)
        stop(wmsg("'length(nzdata)' is greater than 'nrow(nzcoo)'"))
    if (nzdata_len == 0L)
        stop(wmsg("'length(nzdata)' is 0 but 'nrow(nzcoo)' is not"))
    if (length.out %% nzdata_len != 0L)
        warning(wmsg("'nrow(nzcoo)' is not a multiple of 'length(nzdata)'"))
    rep(nzdata, length.out=length.out)
}

COO_SparseArray <- function(dim, nzcoo=NULL, nzdata=NULL, dimnames=NULL,
                                 check=TRUE)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(nzcoo)) {
        if (is.null(nzdata)) {
            nzdata <- logical(0)  # vector()
        } else if (!(is.vector(nzdata) && length(nzdata) == 0L)) {
            stop(wmsg("'nzdata' must be NULL or a vector of length 0 ",
                      "when 'nzcoo' is NULL"))
        }
        nzcoo <- matrix(integer(0), ncol=length(dim))
    } else {
        if (!is.matrix(nzcoo))
            stop(wmsg("'nzcoo' must be a matrix"))
        if (storage.mode(nzcoo) == "double")
            storage.mode(nzcoo) <- "integer"
        if (!is.null(dimnames(nzcoo)))
            dimnames(nzcoo) <- NULL
        nzdata <- .normarg_nzdata(nzdata, nrow(nzcoo))
    }
    dimnames <- normarg_dimnames(dimnames, dim)
    new2("COO_SparseArray", dim=dim, dimnames=dimnames,
                            nzcoo=nzcoo, nzdata=nzdata,
                            check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dense2sparse() and sparse2dense()
###

### 'x' must be an array-like object that supports 'type()' and subsetting
### by an M-index subscript.
### Returns a COO_SparseArray object.
dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    is_not_zero <- x != zero
    nzcoo <- which(is_not_zero | is.na(is_not_zero), arr.ind=TRUE)  # M-index
    COO_SparseArray(x_dim, nzcoo, x[nzcoo], dimnames(x), check=FALSE)
}

### 'coosa' must be a COO_SparseArray object.
### Return an ordinary array.
sparse2dense <- function(coosa)
{
    if (!is(coosa, "COO_SparseArray"))
        stop(wmsg("'coosa' must be a COO_SparseArray object"))
    sa_nzdata <- nzdata(coosa)
    zero <- vector(typeof(sa_nzdata), length=1L)
    ans <- array(zero, dim=dim(coosa))
    ans[nzcoo(coosa)] <- sa_nzdata
    set_dimnames(ans, dimnames(coosa))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The extract_sparse_array() generics
###

### This is the workhorse behind read_sparse_block().
### Similar to extract_array() except that:
###   (1) The extracted array data must be returned in a COO_SparseArray
###       object. Methods should always operate on the sparse representation
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
        stopifnot(is(ans, "COO_SparseArray"),
                  identical(dim(ans), expected_dim))
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse(), extract_sparse_array(), and extract_array() methods
###

### IMPORTANT NOTE: The returned COO_SparseArray object is guaranteed to be
### **correct** ONLY if the subscripts in 'index' do NOT contain duplicates!
### If they contain duplicates, the correct COO_SparseArray object to return
### should contain repeated nonzero data. However, in order to keep it as
### efficient as possible, the code below does NOT repeat the nonzero data
### that corresponds to duplicates subscripts. It does not check for
### duplicates in 'index' either because this check could have a
### non-negligible cost.
### All this is OK because .extract_sparse_array_from_COO_SparseArray()
### should always be used in a context where 'index' does NOT contain
### duplicates. The only situation where 'index' CAN contain duplicates
### is when .extract_sparse_array_from_COO_SparseArray() is called by
### .extract_array_from_COO_SparseArray(), in which case the missing
### nonzero data is added later.
.extract_sparse_array_from_COO_SparseArray <- function(x, index)
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
    ans_nzdata <- x@nzdata[keep_idx]
    COO_SparseArray(ans_dim, ans_nzcoo, ans_nzdata, check=FALSE)
}
setMethod("extract_sparse_array", "COO_SparseArray",
    .extract_sparse_array_from_COO_SparseArray
)

.extract_array_from_COO_SparseArray <- function(x, index)
{
    sa0 <- .extract_sparse_array_from_COO_SparseArray(x, index)
    ## If the subscripts in 'index' contain duplicates, 'sa0' is
    ## "incomplete" in the sense that it does not contain the nonzero data
    ## that should have been repeated according to the duplicates in the
    ## subscripts (see IMPORTANT NOTE above).
    ans0 <- sparse2dense(sa0)
    ## We "complete" 'ans0' by repeating the nonzero data according to the
    ## duplicates present in the subscripts. Note that this is easy and cheap
    ## to do now because 'ans0' uses a dense representation (it's an ordinary
    ## array). This would be much harder to do **natively** on the
    ## COO_SparseArray form (i.e. without converting first to dense then
    ## back to sparse in the process).
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
### Show
###

setMethod("show", "COO_SparseArray",
    function(object) show_compact_array(object)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to/from COO_SparseArray
###

### S3/S4 combo for as.array.COO_SparseArray
as.array.COO_SparseArray <- function(x, ...) sparse2dense(x)
setMethod("as.array", "COO_SparseArray", as.array.COO_SparseArray)

### S3/S4 combo for as.matrix.COO_SparseArray
.from_COO_SparseArray_to_matrix <- function(x)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    sparse2dense(x)
}
as.matrix.COO_SparseArray <-
    function(x, ...) .from_COO_SparseArray_to_matrix(x, ...)
setMethod("as.matrix", "COO_SparseArray", .from_COO_SparseArray_to_matrix)

setAs("ANY", "COO_SparseArray", function(from) dense2sparse(from))

### Going back and forth between COO_SparseArray and dg[C|R]Matrix or
### lg[C|R]Matrix objects from the Matrix package:

.from_COO_SparseArray_to_CsparseMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgCMatrix ",
                  "or lgCMatrix must have exactly 2 dimensions"))
    i <- from@nzcoo[ , 1L]
    j <- from@nzcoo[ , 2L]
    CsparseMatrix(from_dim, i, j, from@nzdata, dimnames=dimnames(from))
}

.from_COO_SparseArray_to_RsparseMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgRMatrix ",
                  "or lgRMatrix must have exactly 2 dimensions"))
    i <- from@nzcoo[ , 1L]
    j <- from@nzcoo[ , 2L]
    RsparseMatrix(from_dim, i, j, from@nzdata, dimnames=dimnames(from))
}

setAs("COO_SparseArray", "CsparseMatrix",
    .from_COO_SparseArray_to_CsparseMatrix
)
setAs("COO_SparseArray", "RsparseMatrix",
    .from_COO_SparseArray_to_RsparseMatrix
)
setAs("COO_SparseArray", "sparseMatrix",
    .from_COO_SparseArray_to_CsparseMatrix
)
setAs("COO_SparseArray", "dgCMatrix",
    function(from) as(as(from, "CsparseMatrix"), "dgCMatrix")
)
setAs("COO_SparseArray", "dgRMatrix",
    function(from) as(as(from, "RsparseMatrix"), "dgRMatrix")
)
setAs("COO_SparseArray", "lgCMatrix",
    function(from) as(as(from, "CsparseMatrix"), "lgCMatrix")
)
### Will fail if 'as(from, "RsparseMatrix")' returns a dgRMatrix object
### because the Matrix package doesn't support coercion from dgRMatrix
### to lgRMatrix at the moment:
###   > as(COO_SparseArray(4:3, rbind(c(4, 3)), -2), "lgRMatrix")
###   Error in as(as(from, "RsparseMatrix"), "lgRMatrix") :
###     no method or default for coercing “dgRMatrix” to “lgRMatrix”
setAs("COO_SparseArray", "lgRMatrix",
    function(from) as(as(from, "RsparseMatrix"), "lgRMatrix")
)

.make_COO_SparseArray_from_dgCMatrix_or_lgCMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- from@i + 1L
    j <- rep.int(seq_len(ncol(from)), diff(from@p))
    ans_nzcoo <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    COO_SparseArray(dim(from), ans_nzcoo, from@x, ans_dimnames, check=FALSE)
}

.make_COO_SparseArray_from_dgRMatrix_or_lgRMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- rep.int(seq_len(nrow(from)), diff(from@p))
    j <- from@j + 1L
    ans_nzcoo <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    COO_SparseArray(dim(from), ans_nzcoo, from@x, ans_dimnames, check=FALSE)
}

setAs("dgCMatrix", "COO_SparseArray",
    function(from) .make_COO_SparseArray_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("dgRMatrix", "COO_SparseArray",
    function(from) .make_COO_SparseArray_from_dgRMatrix_or_lgRMatrix(from)
)
setAs("lgCMatrix", "COO_SparseArray",
    function(from) .make_COO_SparseArray_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("lgRMatrix", "COO_SparseArray",
    function(from) .make_COO_SparseArray_from_dgRMatrix_or_lgRMatrix(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse() and extract_sparse_array() methods for dg[C|R]Matrix and
### lg[C|R]Matrix objects from the Matrix package
###
### TODO: Support more sparseMatrix derivatives (e.g. dgTMatrix, dgRMatrix)
### as the need arises.
###

setMethod("is_sparse", "dgCMatrix", function(x) TRUE)
setMethod("is_sparse", "dgRMatrix", function(x) TRUE)
setMethod("is_sparse", "lgCMatrix", function(x) TRUE)
setMethod("is_sparse", "lgRMatrix", function(x) TRUE)

.extract_sparse_array_from_dgCMatrix_or_lgCMatrix <- function(x, index)
{
    sm <- subset_by_Nindex(x, index)  # a dgCMatrix or lgCMatrix object
    .make_COO_SparseArray_from_dgCMatrix_or_lgCMatrix(sm, use.dimnames=FALSE)
}

.extract_sparse_array_from_dgRMatrix_or_lgRMatrix <- function(x, index)
{
    sm <- subset_by_Nindex(x, index)  # a dgRMatrix or lgRMatrix object
    .make_COO_SparseArray_from_dgRMatrix_or_lgRMatrix(sm, use.dimnames=FALSE)
}

setMethod("extract_sparse_array", "dgCMatrix",
    .extract_sparse_array_from_dgCMatrix_or_lgCMatrix
)
setMethod("extract_sparse_array", "dgRMatrix",
    .extract_sparse_array_from_dgRMatrix_or_lgRMatrix
)
setMethod("extract_sparse_array", "lgCMatrix",
    .extract_sparse_array_from_dgCMatrix_or_lgCMatrix
)
setMethod("extract_sparse_array", "lgRMatrix",
    .extract_sparse_array_from_dgRMatrix_or_lgRMatrix
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Extend base::aperm() by allowing dropping and/or adding ineffective
### dimensions. See aperm2.R
###

.aperm.COO_SparseArray <- function(a, perm)
{
    a_dim <- dim(a)
    perm <- normarg_perm(perm, a_dim)
    msg <- validate_perm(perm, a_dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    new_dim <- a_dim[perm]
    new_dim[is.na(perm)] <- 1L
    new_nzcoo <- a@nzcoo[ , perm, drop=FALSE]
    new_nzcoo[ , is.na(perm)] <- 1L
    new_dimnames <- a@dimnames[perm]
    BiocGenerics:::replaceSlots(a, dim=new_dim,
                                   dimnames=new_dimnames,
                                   nzcoo=new_nzcoo,
                                   check=FALSE)
}

### S3/S4 combo for aperm.COO_SparseArray
aperm.COO_SparseArray <-
    function(a, perm, ...) .aperm.COO_SparseArray(a, perm, ...)
setMethod("aperm", "COO_SparseArray", aperm.COO_SparseArray)

