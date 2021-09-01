### =========================================================================
### SparseArray objects
### -------------------------------------------------------------------------


setClass("SparseArray",
    contains="Array",
    representation(
        dim="integer",     # This gives us dim() for free!
        nzindex="matrix",  # M-index of the nonzero data.
        nzdata="vector",   # A vector (atomic or list) of length
                           # 'nrow(nzindex)' containing the nonzero data.
        dimnames="list"    # List with one list element per dimension. Each
                           # list element must be NULL or a character vector.
    ),
    prototype(
        dim=0L,
        nzindex=matrix(integer(0), ncol=1L),
        dimnames=list(NULL)
    )
)

### API:
### - Getters/setters: dim(), length(), nzindex(), nzdata(),
###                    dimnames(), dimnames<-()
### - sparsity().
### - dense2sparse(), sparse2dense().
### - Based on sparse2dense(): extract_array(), as.array(), as.matrix().
### - Based on dense2sparse(): coercion to SparseArray.
### - Back and forth coercion between SparseArray and dg[C|R]Matrix or
###   lg[C|R]Matrix objects from the Matrix package.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_nzindex_slot <- function(x)
{
    x_nzindex <- x@nzindex
    if (!(is.matrix(x_nzindex) && typeof(x_nzindex) == "integer"))
        return("'nzindex' slot must be an integer matrix")
    x_dim <- x@dim
    if (ncol(x_nzindex) != length(x_dim))
        return(paste0("'nzindex' slot must be a matrix with ",
                      "one column per dimension"))
    for (along in seq_along(x_dim)) {
        not_ok <- S4Vectors:::anyMissingOrOutside(x_nzindex[ , along],
                                                  1L, x_dim[[along]])
        if (not_ok)
            return(paste0("'nzindex' slot must contain valid indices, ",
                          "that is, indices that are not NA and are ",
                          ">= 1 and <= their corresponding dimension"))
    }
    TRUE
}

.validate_nzdata_slot <- function(x)
{
    x_nzdata <- x@nzdata
    if (!(is.vector(x_nzdata) && length(x_nzdata) == nrow(x@nzindex)))
        return(paste0("'nzdata' slot must be a vector of length ",
                      "the number of rows in the 'nzindex' slot"))
    TRUE
}

.validate_SparseArray <- function(x)
{
    msg <- validate_dim_slot(x, "dim")
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzindex_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzdata_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- validate_dimnames_slot(x, x@dim)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}

setValidity2("SparseArray", .validate_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters/setters
###

setGeneric("nzindex", function(x) standardGeneric("nzindex"))
setMethod("nzindex", "SparseArray", function(x) x@nzindex)

setGeneric("nzdata", function(x) standardGeneric("nzdata"))
setMethod("nzdata", "SparseArray", function(x) x@nzdata)

setMethod("dimnames", "SparseArray",
    function(x) simplify_NULL_dimnames(x@dimnames)
)

setReplaceMethod("dimnames", "SparseArray",
    function(x, value)
    {
        x@dimnames <- normarg_dimnames(value, dim(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sparsity()
###

setGeneric("sparsity", function(x) standardGeneric("sparsity"))

setMethod("sparsity", "SparseArray",
    function(x) { 1 - length(nzdata(x)) / length(x) }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.normarg_nzdata <- function(nzdata, length.out)
{
    if (is.null(nzdata))
        stop(wmsg("'nzdata' cannot be NULL when 'nzindex' is not NULL"))
    if (!is.vector(nzdata))
        stop(wmsg("'nzdata' must be a vector"))
    ## Same logic as S4Vectors:::V_recycle().
    nzdata_len <- length(nzdata)
    if (nzdata_len == length.out)
        return(nzdata)
    if (nzdata_len > length.out && nzdata_len != 1L)
        stop(wmsg("'length(nzdata)' is greater than 'nrow(nzindex)'"))
    if (nzdata_len == 0L)
        stop(wmsg("'length(nzdata)' is 0 but 'nrow(nzindex)' is not"))
    if (length.out %% nzdata_len != 0L)
        warning(wmsg("'nrow(nzindex)' is not a multiple of 'length(nzdata)'"))
    rep(nzdata, length.out=length.out)
}

SparseArray <- function(dim, nzindex=NULL, nzdata=NULL, dimnames=NULL,
                            check=TRUE)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(nzindex)) {
        if (is.null(nzdata)) {
            nzdata <- logical(0)  # vector()
        } else if (!(is.vector(nzdata) && length(nzdata) == 0L)) {
            stop(wmsg("'nzdata' must be NULL or a vector of length 0 ",
                      "when 'nzindex' is NULL"))
        }
        nzindex <- matrix(integer(0), ncol=length(dim))
    } else {
        if (!is.matrix(nzindex))
            stop(wmsg("'nzindex' must be a matrix"))
        if (storage.mode(nzindex) == "double")
            storage.mode(nzindex) <- "integer"
        if (!is.null(dimnames(nzindex)))
            dimnames(nzindex) <- NULL
        nzdata <- .normarg_nzdata(nzdata, nrow(nzindex))
    }
    dimnames <- normarg_dimnames(dimnames, dim)
    new2("SparseArray", dim=dim, nzindex=nzindex, nzdata=nzdata,
                            dimnames=dimnames,
                            check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dense2sparse() and sparse2dense()
###

### 'x' must be an array-like object that supports 'type()' and subsetting
### by an M-index subscript.
### Return a SparseArray object.
dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    nzindex <- which(x != zero, arr.ind=TRUE)  # M-index
    SparseArray(x_dim, nzindex, x[nzindex], dimnames(x), check=FALSE)
}

### 'sa' must be a SparseArray object.
### Return an ordinary array.
sparse2dense <- function(sa)
{
    if (!is(sa, "SparseArray"))
        stop(wmsg("'sa' must be a SparseArray object"))
    sa_nzdata <- nzdata(sa)
    zero <- vector(typeof(sa_nzdata), length=1L)
    ans <- array(zero, dim=dim(sa))
    ans[nzindex(sa)] <- sa_nzdata
    set_dimnames(ans, dimnames(sa))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The is_sparse() and extract_sparse_array() generics
###

### is_sparse() detects **structural** sparsity which is a global qualitative
### property of array-like object 'x' rather than a quantitative one.
### In other words it doesn't look at the data in 'x' to decide whether 'x'
### should be considered sparse or not. Said otherwise, it is NOT about
### quantitative sparsity measured by sparsity().
### IMPORTANT: Seeds for which is_sparse() returns TRUE **must** support
### extract_sparse_array().
setGeneric("is_sparse", function(x) standardGeneric("is_sparse"))

setGeneric("is_sparse<-", signature="x",
    function(x, value) standardGeneric("is_sparse<-")
)

### By default, nothing is considered sparse.
setMethod("is_sparse", "ANY", function(x) FALSE)

### This is the workhorse behind read_sparse_block().
### Similar to extract_array() except that:
###   (1) The extracted array data must be returned in a SparseArray
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
        stopifnot(is(ans, "SparseArray"),
                  identical(dim(ans), expected_dim))
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse(), extract_sparse_array(), and extract_array() methods for
### SparseArray objects
###

setMethod("is_sparse", "SparseArray", function(x) TRUE)

### IMPORTANT NOTE: The returned SparseArray object is guaranteed to be
### **correct** ONLY if the subscripts in 'index' do NOT contain duplicates!
### If they contain duplicates, the correct SparseArray object to return
### should contain repeated nonzero data. However, in order to keep it as
### efficient as possible, the code below does NOT repeat the nonzero data
### that corresponds to duplicates subscripts. It does not check for
### duplicates in 'index' either because this check could have a
### non-negligible cost.
### All this is OK because .extract_sparse_array_from_SparseArray()
### should always be used in a context where 'index' does NOT contain
### duplicates. The only situation where 'index' CAN contain duplicates
### is when .extract_sparse_array_from_SparseArray() is called by
### .extract_array_from_SparseArray(), in which case the missing
### nonzero data is added later.
.extract_sparse_array_from_SparseArray <- function(x, index)
{
    stopifnot(is(x, "SparseArray"))
    ans_dim <- get_Nindex_lengths(index, dim(x))
    x_nzindex <- x@nzindex
    for (along in seq_along(ans_dim)) {
        i <- index[[along]]
        if (is.null(i))
            next
        x_nzindex[ , along] <- match(x_nzindex[ , along], i)
    }
    keep_idx <- which(!rowAnyNAs(x_nzindex))
    ans_nzindex <- x_nzindex[keep_idx, , drop=FALSE]
    ans_nzdata <- x@nzdata[keep_idx]
    SparseArray(ans_dim, ans_nzindex, ans_nzdata, check=FALSE)
}
setMethod("extract_sparse_array", "SparseArray",
    .extract_sparse_array_from_SparseArray
)

.extract_array_from_SparseArray <- function(x, index)
{
    sa0 <- .extract_sparse_array_from_SparseArray(x, index)
    ## If the subscripts in 'index' contain duplicates, 'sa0' is
    ## "incomplete" in the sense that it does not contain the nonzero data
    ## that should have been repeated according to the duplicates in the
    ## subscripts (see IMPORTANT NOTE above).
    ans0 <- sparse2dense(sa0)
    ## We "complete" 'ans0' by repeating the nonzero data according to the
    ## duplicates present in the subscripts. Note that this is easy and cheap
    ## to do now because 'ans0' uses a dense representation (it's an ordinary
    ## array). This would be much harder to do **natively** on the
    ## SparseArray form (i.e. without converting first to dense then
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
setMethod("extract_array", "SparseArray",
    .extract_array_from_SparseArray
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "SparseArray",
    function(object) show_compact_array(object)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to/from SparseArray
###

### S3/S4 combo for as.array.SparseArray
as.array.SparseArray <- function(x, ...) sparse2dense(x)
setMethod("as.array", "SparseArray", as.array.SparseArray)

### S3/S4 combo for as.matrix.SparseArray
.from_SparseArray_to_matrix <- function(x)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop(wmsg("'x' must have exactly 2 dimensions"))
    sparse2dense(x)
}
as.matrix.SparseArray <-
    function(x, ...) .from_SparseArray_to_matrix(x, ...)
setMethod("as.matrix", "SparseArray", .from_SparseArray_to_matrix)

setAs("ANY", "SparseArray", function(from) dense2sparse(from))

### Going back and forth between SparseArray and dg[C|R]Matrix or
### lg[C|R]Matrix objects from the Matrix package:

.from_SparseArray_to_CsparseMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgCMatrix ",
                  "or lgCMatrix must have exactly 2 dimensions"))
    i <- from@nzindex[ , 1L]
    j <- from@nzindex[ , 2L]
    CsparseMatrix(from_dim, i, j, from@nzdata, dimnames=dimnames(from))
}

.from_SparseArray_to_RsparseMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) != 2L)
        stop(wmsg("the ", class(from), " object to coerce to dgRMatrix ",
                  "or lgRMatrix must have exactly 2 dimensions"))
    i <- from@nzindex[ , 1L]
    j <- from@nzindex[ , 2L]
    RsparseMatrix(from_dim, i, j, from@nzdata, dimnames=dimnames(from))
}

setAs("SparseArray", "CsparseMatrix",
    .from_SparseArray_to_CsparseMatrix
)
setAs("SparseArray", "RsparseMatrix",
    .from_SparseArray_to_RsparseMatrix
)
setAs("SparseArray", "sparseMatrix",
    .from_SparseArray_to_CsparseMatrix
)
setAs("SparseArray", "dgCMatrix",
    function(from) as(as(from, "CsparseMatrix"), "dgCMatrix")
)
setAs("SparseArray", "dgRMatrix",
    function(from) as(as(from, "RsparseMatrix"), "dgRMatrix")
)
setAs("SparseArray", "lgCMatrix",
    function(from) as(as(from, "CsparseMatrix"), "lgCMatrix")
)
### Will fail if 'as(from, "RsparseMatrix")' returns a dgRMatrix object
### because the Matrix package doesn't support coercion from dgRMatrix
### to lgRMatrix at the moment:
###   > as(SparseArray(4:3, rbind(c(4, 3)), -2), "lgRMatrix")
###   Error in as(as(from, "RsparseMatrix"), "lgRMatrix") :
###     no method or default for coercing “dgRMatrix” to “lgRMatrix”
setAs("SparseArray", "lgRMatrix",
    function(from) as(as(from, "RsparseMatrix"), "lgRMatrix")
)

.make_SparseArray_from_dgCMatrix_or_lgCMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- from@i + 1L
    j <- rep.int(seq_len(ncol(from)), diff(from@p))
    ans_nzindex <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    SparseArray(dim(from), ans_nzindex, from@x, ans_dimnames, check=FALSE)
}

.make_SparseArray_from_dgRMatrix_or_lgRMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- rep.int(seq_len(nrow(from)), diff(from@p))
    j <- from@j + 1L
    ans_nzindex <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    SparseArray(dim(from), ans_nzindex, from@x, ans_dimnames, check=FALSE)
}

setAs("dgCMatrix", "SparseArray",
    function(from) .make_SparseArray_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("dgRMatrix", "SparseArray",
    function(from) .make_SparseArray_from_dgRMatrix_or_lgRMatrix(from)
)
setAs("lgCMatrix", "SparseArray",
    function(from) .make_SparseArray_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("lgRMatrix", "SparseArray",
    function(from) .make_SparseArray_from_dgRMatrix_or_lgRMatrix(from)
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
    .make_SparseArray_from_dgCMatrix_or_lgCMatrix(sm, use.dimnames=FALSE)
}

.extract_sparse_array_from_dgRMatrix_or_lgRMatrix <- function(x, index)
{
    sm <- subset_by_Nindex(x, index)  # a dgRMatrix or lgRMatrix object
    .make_SparseArray_from_dgRMatrix_or_lgRMatrix(sm, use.dimnames=FALSE)
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

.aperm.SparseArray <- function(a, perm)
{
    a_dim <- dim(a)
    perm <- normarg_perm(perm, a_dim)
    msg <- validate_perm(perm, a_dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    new_dim <- a_dim[perm]
    new_dim[is.na(perm)] <- 1L
    new_nzindex <- a@nzindex[ , perm, drop=FALSE]
    new_nzindex[ , is.na(perm)] <- 1L
    new_dimnames <- a@dimnames[perm]
    BiocGenerics:::replaceSlots(a, dim=new_dim,
                                   nzindex=new_nzindex,
                                   dimnames=new_dimnames,
                                   check=FALSE)
}

### S3/S4 combo for aperm.SparseArray
aperm.SparseArray <-
    function(a, perm, ...) .aperm.SparseArray(a, perm, ...)
setMethod("aperm", "SparseArray", aperm.SparseArray)

