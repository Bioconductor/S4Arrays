### 'a0' is expected to be of type "double".
### We only check types "double", "integer", "logical", and "raw" at the
### moment. No "complex", "character", or "list" (even though they are
### supported).
.test_coercion_to_SparseArray_with_various_types <-
    function(a0, to, expected_class)
{
    object <- as(a0, to)
    check_SparseArray_object(object, expected_class, a0)
    a <- a0
    suppressWarnings(storage.mode(a) <- "integer")
    suppressWarnings(object <- as(a, to))
    check_SparseArray_object(object, expected_class, a)
    a <- a0
    suppressWarnings(storage.mode(a) <- "logical")
    suppressWarnings(object <- as(a, to))
    check_SparseArray_object(object, expected_class, a)
    a <- a0
    suppressWarnings(storage.mode(a) <- "raw")
    suppressWarnings(object <- as(a, to))
    check_SparseArray_object(object, expected_class, a)
}

.test_SparseMatrix_transposition <- function(m0, to)
{
    svt <- as(m0, to)
    tm0 <- t(m0)
    tsvt <- t(svt)
    check_SparseArray_object(tsvt, to, tm0)
    expect_identical(tsvt, as(tm0, to))
}

test_that("array <==> SVT_SparseArray coercions", {
    ## Only zeros.
    a0 <- array(0.0, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    .test_coercion_to_SparseArray_with_various_types(a0, "SVT_SparseArray",
                                                         "SVT_SparseArray")
    .test_coercion_to_SparseArray_with_various_types(a0, "SparseArray",
                                                         "SVT_SparseArray")
    m0 <- a0[ , , 1]  # 2D
    .test_coercion_to_SparseArray_with_various_types(m0, "SVT_SparseMatrix",
                                                         "SVT_SparseMatrix")
    .test_coercion_to_SparseArray_with_various_types(m0, "SVT_SparseArray",
                                                         "SVT_SparseMatrix")
    .test_coercion_to_SparseArray_with_various_types(m0, "SparseArray",
                                                         "SVT_SparseMatrix")
    x0 <- as.array(m0[1, ])  # 1D
    .test_coercion_to_SparseArray_with_various_types(x0, "SVT_SparseArray",
                                                         "SVT_SparseArray")
    .test_coercion_to_SparseArray_with_various_types(x0, "SparseArray",
                                                         "SVT_SparseArray")

    ## Add some nonzero values.
    set.seed(123)
    a0[5*(1:42)] <- runif(42, min=-5, max=10)
    a0[2, c(1:4, 7:9), 1] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    .test_coercion_to_SparseArray_with_various_types(a0, "SVT_SparseArray",
                                                         "SVT_SparseArray")
    m0 <- a0[ , , 1]  # 2D
    .test_coercion_to_SparseArray_with_various_types(m0, "SVT_SparseMatrix",
                                                         "SVT_SparseMatrix")
    x0 <- as.array(m0[2, ])  # 1D
    .test_coercion_to_SparseArray_with_various_types(x0, "SVT_SparseArray",
                                                         "SVT_SparseArray")

    ## Length zero.
    a <- a0[ , 0, ]
    .test_coercion_to_SparseArray_with_various_types(a, "SVT_SparseArray",
                                                        "SVT_SparseArray")
})

test_that("dgCMatrix <==> SVT_SparseMatrix coercions", {
    ## Only zeros.
    m0 <- matrix(0.0, nrow=7, ncol=10, dimnames=list(NULL, letters[1:10]))
    dgcm0 <- as(m0, "dgCMatrix")
    svt <- as(dgcm0, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m0)
    expect_identical(as(m0, "SVT_SparseMatrix"), svt)
    expect_identical(as(dgcm0, "SparseMatrix"), svt)
    expect_identical(as(dgcm0, "SVT_SparseArray"), svt)
    expect_identical(as(dgcm0, "SparseArray"), svt)
    expect_identical(as(svt, "dgCMatrix"), dgcm0)
    expect_identical(as(svt, "CsparseMatrix"), dgcm0)
    expect_identical(as(svt, "sparseMatrix"), dgcm0)

    ## Add some nonzero values.
    set.seed(456)
    m0[5*(1:14)] <- runif(7, min=-5, max=10)
    m0[2, c(1:4, 7:9)] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    dgcm0 <- as(m0, "dgCMatrix")
    svt <- as(dgcm0, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m0)
    expect_identical(as(m0, "SVT_SparseMatrix"), svt)
    expect_identical(as(dgcm0, "SparseMatrix"), svt)
    expect_identical(as(dgcm0, "SVT_SparseArray"), svt)
    expect_identical(as(dgcm0, "SparseArray"), svt)
    expect_identical(as(svt, "dgCMatrix"), dgcm0)
    expect_identical(as(svt, "CsparseMatrix"), dgcm0)
    expect_identical(as(svt, "sparseMatrix"), dgcm0)

    ## Length zero.
    m <- m0[0 , ]
    dgcm <- as(m, "dgCMatrix")
    svt <- as(dgcm, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)
    expect_identical(as(dgcm, "SparseMatrix"), svt)
    expect_identical(as(dgcm, "SVT_SparseArray"), svt)
    expect_identical(as(dgcm, "SparseArray"), svt)
    expect_identical(as(svt, "dgCMatrix"), dgcm)
    expect_identical(as(svt, "CsparseMatrix"), dgcm)
    expect_identical(as(svt, "sparseMatrix"), dgcm)
    m <- m0[ , 0]  # this sets the dimnames to list(NULL, NULL)
    dimnames(m) <- NULL  # fix the dimnames
    dgcm <- as(m, "dgCMatrix")
    svt <- as(dgcm, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)
    expect_identical(as(dgcm, "SparseMatrix"), svt)
    expect_identical(as(dgcm, "SVT_SparseArray"), svt)
    expect_identical(as(dgcm, "SparseArray"), svt)
    expect_identical(as(svt, "dgCMatrix"), dgcm)
    expect_identical(as(svt, "CsparseMatrix"), dgcm)
    expect_identical(as(svt, "sparseMatrix"), dgcm)
})

test_that("lgCMatrix <==> SVT_SparseMatrix coercions", {
    ## Only zeros.
    m0 <- matrix(FALSE, nrow=7, ncol=10, dimnames=list(NULL, letters[1:10]))
    lgcm0 <- as(m0, "lgCMatrix")
    svt <- as(lgcm0, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m0)
    expect_identical(as(m0, "SVT_SparseMatrix"), svt)
    expect_identical(as(lgcm0, "SparseMatrix"), svt)
    expect_identical(as(lgcm0, "SVT_SparseArray"), svt)
    expect_identical(as(lgcm0, "SparseArray"), svt)
    expect_identical(as(svt, "lgCMatrix"), lgcm0)
    expect_identical(as(svt, "CsparseMatrix"), lgcm0)
    expect_identical(as(svt, "sparseMatrix"), lgcm0)

    ## Add some nonzero values.
    m0[5*(1:14)] <- TRUE
    m0[17*(1:4)] <- NA
    m0[3, 8] <- NA
    lgcm0 <- as(m0, "lgCMatrix")
    svt <- as(lgcm0, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m0)
    expect_identical(as(m0, "SVT_SparseMatrix"), svt)
    expect_identical(as(lgcm0, "SparseMatrix"), svt)
    expect_identical(as(lgcm0, "SVT_SparseArray"), svt)
    expect_identical(as(lgcm0, "SparseArray"), svt)
    expect_identical(as(svt, "lgCMatrix"), lgcm0)
    expect_identical(as(svt, "CsparseMatrix"), lgcm0)
    expect_identical(as(svt, "sparseMatrix"), lgcm0)

    ## Length zero.
    m <- m0[0 , ]
    lgcm <- as(m, "lgCMatrix")
    svt <- as(lgcm, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)
    expect_identical(as(lgcm, "SparseMatrix"), svt)
    expect_identical(as(lgcm, "SVT_SparseArray"), svt)
    expect_identical(as(lgcm, "SparseArray"), svt)
    expect_identical(as(svt, "lgCMatrix"), lgcm)
    expect_identical(as(svt, "CsparseMatrix"), lgcm)
    expect_identical(as(svt, "sparseMatrix"), lgcm)
    m <- m0[ , 0]  # this sets the dimnames to list(NULL, NULL)
    dimnames(m) <- NULL  # fix the dimnames
    lgcm <- as(m, "lgCMatrix")
    svt <- as(lgcm, "SVT_SparseMatrix")
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)
    expect_identical(as(lgcm, "SparseMatrix"), svt)
    expect_identical(as(lgcm, "SVT_SparseArray"), svt)
    expect_identical(as(lgcm, "SparseArray"), svt)
    expect_identical(as(svt, "lgCMatrix"), lgcm)
    expect_identical(as(svt, "CsparseMatrix"), lgcm)
    expect_identical(as(svt, "sparseMatrix"), lgcm)
})

test_that("SVT_SparseMatrix transposition", {
    ## Only zeros.
    m0 <- matrix(0.0, nrow=7, ncol=10, dimnames=list(NULL, letters[1:10]))
    .test_SparseMatrix_transposition(m0, "SVT_SparseMatrix")

    ## Add some nonzero values.
    set.seed(789)
    m0[5*(1:14)] <- runif(7, min=-5, max=10)
    m0[2, c(1:4, 7:9)] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    .test_SparseMatrix_transposition(m0, "SVT_SparseMatrix")

    ## Length zero.
    m <- m0[0 , ]
    .test_SparseMatrix_transposition(m, "SVT_SparseMatrix")
    m <- m0[ , 0]  # this sets the dimnames to list(NULL, NULL)
    dimnames(m) <- NULL  # fix the dimnames
    .test_SparseMatrix_transposition(m, "SVT_SparseMatrix")

    ## Other types.
    m <- m0
    suppressWarnings(storage.mode(m) <- "integer")
    .test_SparseMatrix_transposition(m, "SVT_SparseMatrix")
    m <- m0
    suppressWarnings(storage.mode(m) <- "logical")
    .test_SparseMatrix_transposition(m, "SVT_SparseMatrix")
    m <- m0
    suppressWarnings(storage.mode(m) <- "raw")
    .test_SparseMatrix_transposition(m, "SVT_SparseMatrix")
})

