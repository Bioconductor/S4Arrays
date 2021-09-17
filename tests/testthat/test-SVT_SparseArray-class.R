check_SparseArray_object <- function(object, expected_class, a0)
{
    expect_s4_class(object, expected_class)
    expect_true(validObject(object))
    expect_identical(dim(object), dim(a0))
    expect_identical(dimnames(object), dimnames(a0))
    expect_identical(type(object), type(a0))
    a <- as.array(object)
    expect_identical(a, a0)
}

### 'a0' is expected to be of type "double".
### We only check types "double", "integer", "logical", and "raw" at the
### moment. No "complex", "character", or "list" (even though they are
### supported).
check_coercion_to_SparseArray_with_various_types <- function(a0, to)
{
    object <- as(a0, to)
    check_SparseArray_object(object, to, a0)
    suppressWarnings(storage.mode(a0) <- "integer")
    suppressWarnings(object <- as(a0, to))
    check_SparseArray_object(object, to, a0)
    suppressWarnings(storage.mode(a0) <- "logical")
    suppressWarnings(object <- as(a0, to))
    check_SparseArray_object(object, to, a0)
    suppressWarnings(storage.mode(a0) <- "raw")
    suppressWarnings(object <- as(a0, to))
    check_SparseArray_object(object, to, a0)
}

test_that("array <==> SVT_SparseArray coercions", {
    ## Only zeros.
    a0 <- array(0.0, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    check_coercion_to_SparseArray_with_various_types(a0, "SVT_SparseArray")
    m0 <- a0[ , , 1]  # 2D
    check_coercion_to_SparseArray_with_various_types(m0, "SVT_SparseArray")
    x0 <- m0[1, , drop=FALSE]  # 1D
    check_coercion_to_SparseArray_with_various_types(x0, "SVT_SparseArray")

    ## Add some nonzero values.
    a0[5*(1:42)] <- runif(42, min=-5, max=10)
    a0[2, c(1:4, 7:9), 1] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    check_coercion_to_SparseArray_with_various_types(a0, "SVT_SparseArray")
    m0 <- a0[ , , 1]  # 2D
    check_coercion_to_SparseArray_with_various_types(m0, "SVT_SparseArray")
    x0 <- m0[2, , drop=FALSE]  # 1D
    check_coercion_to_SparseArray_with_various_types(x0, "SVT_SparseArray")

    ## Length zero.
    a <- a0[ , 0, ]
    check_coercion_to_SparseArray_with_various_types(a, "SVT_SparseArray")
})

test_that("dgCMatrix <==> SVT_SparseArray coercions", {
    ## Only zeros.
    m0 <- matrix(0.0, nrow=7, ncol=10, dimnames=list(NULL, letters[1:10]))
    dgcm0 <- as(m0, "dgCMatrix")
    svt <- as(dgcm0, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m0)
    expect_identical(as(svt, "dgCMatrix"), dgcm0)
    expect_identical(as(svt, "CsparseMatrix"), dgcm0)
    expect_identical(as(svt, "sparseMatrix"), dgcm0)

    ## Add some nonzero values.
    m0[5*(1:14)] <- runif(7, min=-5, max=10)
    m0[2, c(1:4, 7:9)] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    dgcm0 <- as(m0, "dgCMatrix")
    svt <- as(dgcm0, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m0)
    expect_identical(as(svt, "dgCMatrix"), dgcm0)
    expect_identical(as(svt, "CsparseMatrix"), dgcm0)
    expect_identical(as(svt, "sparseMatrix"), dgcm0)

    ## Length zero.
    m <- m0[0 , ]
    dgcm <- as(m, "dgCMatrix")
    svt <- as(dgcm, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m)
    expect_identical(as(svt, "dgCMatrix"), dgcm)
    expect_identical(as(svt, "CsparseMatrix"), dgcm)
    expect_identical(as(svt, "sparseMatrix"), dgcm)
    m <- m0[ , 0]  # this sets the dimnames to list(NULL, NULL)
    dimnames(m) <- NULL  # fix the dimnames
    dgcm <- as(m, "dgCMatrix")
    svt <- as(dgcm, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m)
    expect_identical(as(svt, "dgCMatrix"), dgcm)
    expect_identical(as(svt, "CsparseMatrix"), dgcm)
    expect_identical(as(svt, "sparseMatrix"), dgcm)
})

test_that("lgCMatrix <==> SVT_SparseArray coercions", {
    ## Only zeros.
    m0 <- matrix(FALSE, nrow=7, ncol=10, dimnames=list(NULL, letters[1:10]))
    lgcm0 <- as(m0, "lgCMatrix")
    svt <- as(lgcm0, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m0)
    expect_identical(as(svt, "lgCMatrix"), lgcm0)
    expect_identical(as(svt, "CsparseMatrix"), lgcm0)
    expect_identical(as(svt, "sparseMatrix"), lgcm0)

    ## Add some nonzero values.
    m0[5*(1:14)] <- TRUE
    m0[17*(1:4)] <- NA
    m0[3, 8] <- NA
    lgcm0 <- as(m0, "lgCMatrix")
    svt <- as(lgcm0, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m0)
    expect_identical(as(svt, "lgCMatrix"), lgcm0)
    expect_identical(as(svt, "CsparseMatrix"), lgcm0)
    expect_identical(as(svt, "sparseMatrix"), lgcm0)

    ## Length zero.
    m <- m0[0 , ]
    lgcm <- as(m, "lgCMatrix")
    svt <- as(lgcm, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m)
    expect_identical(as(svt, "lgCMatrix"), lgcm)
    expect_identical(as(svt, "CsparseMatrix"), lgcm)
    expect_identical(as(svt, "sparseMatrix"), lgcm)
    m <- m0[ , 0]  # this sets the dimnames to list(NULL, NULL)
    dimnames(m) <- NULL  # fix the dimnames
    lgcm <- as(m, "lgCMatrix")
    svt <- as(lgcm, "SVT_SparseArray")
    check_SparseArray_object(svt, "SVT_SparseArray", m)
    expect_identical(as(svt, "lgCMatrix"), lgcm)
    expect_identical(as(svt, "CsparseMatrix"), lgcm)
    expect_identical(as(svt, "sparseMatrix"), lgcm)
})

