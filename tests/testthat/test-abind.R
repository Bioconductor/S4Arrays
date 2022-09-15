.TEST_matrices <- list(
    matrix(1:15, nrow=3, ncol=5,
           dimnames=list(NULL, paste0("M1y", 1:5))),
    matrix(101:135, nrow=7, ncol=5,
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5))),
    matrix(1001:1025, nrow=5, ncol=5,
           dimnames=list(paste0("M3x", 1:5), NULL))
)

.TEST_arrays <- list(
    array(1:60, c(3, 5, 4),
           dimnames=list(NULL, paste0("M1y", 1:5), NULL)),
    array(101:240, c(7, 5, 4),
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5), NULL)),
    array(10001:10100, c(5, 5, 4),
           dimnames=list(paste0("M3x", 1:5), NULL, paste0("M3z", 1:4)))
)

test_that("arbind() on matrices", {
    ## on non-empty matrices
    matrices <- .TEST_matrices
    expect_identical(do.call(arbind, matrices), do.call(rbind, matrices))

    ## on empty matrices
    m1 <- matrix(nrow=0, ncol=3, dimnames=list(NULL, letters[1:3]))
    m2 <- matrix(1:15, ncol=3, dimnames=list(NULL, LETTERS[1:3]))
    expect_identical(arbind(m1, m2), rbind(m1, m2))
    expect_identical(arbind(m2, m1), rbind(m2, m1))
    expect_identical(arbind(m1, m1), rbind(m1, m1))

    ## on matrices of type "list"
    m3 <- matrix(list(), nrow=0, ncol=5)
    m4 <- matrix(list(10:9, NULL, letters[1:3], TRUE, raw()), nrow=4, ncol=5)
    expect_identical(arbind(m3, m4), rbind(m3, m4))
    expect_identical(arbind(m4, m3), rbind(m4, m3))
    expect_identical(arbind(m3, m3), rbind(m3, m3))
    matrices <- c(matrices, list(m3, m4))
    expect_identical(do.call(arbind, matrices), do.call(rbind, matrices))
})

test_that("acbind() on matrices", {
    ## on non-empty matrices
    matrices <- lapply(.TEST_matrices, t)
    expect_identical(do.call(acbind, matrices), do.call(cbind, matrices))

    ## on empty matrices
    m1 <- matrix(nrow=3, ncol=0, dimnames=list(letters[1:3], NULL))
    m2 <- matrix(1:15, nrow=3, dimnames=list(LETTERS[1:3], NULL))
    expect_identical(acbind(m1, m2), cbind(m1, m2))
    expect_identical(acbind(m2, m1), cbind(m2, m1))
    expect_identical(acbind(m1, m1), cbind(m1, m1))

    ## on matrices of type "list"
    m3 <- matrix(list(), nrow=5, ncol=0)
    m4 <- matrix(list(10:9, NULL, letters[1:3], TRUE, raw()), nrow=5, ncol=4)
    expect_identical(acbind(m3, m4), cbind(m3, m4))
    expect_identical(acbind(m4, m3), cbind(m4, m3))
    expect_identical(acbind(m3, m3), cbind(m3, m3))
    matrices <- c(matrices, list(m3, m4))
    expect_identical(do.call(acbind, matrices), do.call(cbind, matrices))
    m5 <- matrix(list(), nrow=5, ncol=10)
    expect_identical(acbind(m3, m5), cbind(m3, m5))
    expect_identical(acbind(m5, m3), cbind(m5, m3))
    matrices <- c(matrices, list(m3, m4, m5))
    expect_identical(do.call(acbind, matrices), do.call(cbind, matrices))
})

test_that("arbind() on arrays", {
    ## on 3D arrays
    current <- do.call(arbind, .TEST_arrays)
    check_2D_slice <- function(k) {
        slices <- lapply(.TEST_arrays, `[`, , , k)
        expected_slice <- do.call(rbind, slices)
        expect_identical(current[ , , k], expected_slice)
    }
    for (k in seq_len(dim(current)[[3L]])) check_2D_slice(k)

    ## on 1D arrays
    a1 <- array(11:15, 5, dimnames=list(LETTERS[1:5]))
    expect_identical(arbind(a1), a1)                # unary op
    b1 <- array(letters[1:3])
    expected <- array(c(a1, b1), 8, dimnames=list(c(LETTERS[1:5], rep("", 3))))
    expect_identical(arbind(a1, b1), expected)      # binary op
    expected <- array(c(b1, a1), 8, dimnames=list(c(rep("", 3), LETTERS[1:5])))
    expect_identical(arbind(b1, a1), expected)      # binary op
    a1b1a1 <- arbind(a1, b1, a1)                    # ternary op
    expect_identical(arbind(a1, arbind(b1, a1)), a1b1a1)
    expect_identical(arbind(arbind(a1, b1), a1), a1b1a1)
})

test_that("acbind() on arrays", {
    ## transpose first 2 dimensions of arrays in .TEST_arrays
    arrays <- lapply(.TEST_arrays,
        function(a) {
            a_dimnames <- dimnames(a)
            dim(a)[1:2] <- dim(a)[2:1]
            a_dimnames[1:2] <- a_dimnames[2:1]
            dimnames(a) <- a_dimnames
            a
    })

    ## on 3D arrays
    current <- do.call(acbind, arrays)
    check_2D_slice <- function(k) {
        slices <- lapply(arrays, `[`, , , k)
        expected_slice <- do.call(cbind, slices)
        expect_identical(current[ , , k], expected_slice)
    }
    for (k in seq_len(dim(current)[[3L]])) check_2D_slice(k)

    ## acbind() is not supported on 1D arrays
    expected_msg <- "objects to bind must have at least 2 dimensions"
    a1 <- array(11:15, 5, dimnames=list(LETTERS[1:5]))
    expect_error(acbind(a1), expected_msg)          # unary op
    b1 <- array(letters[1:3])
    expect_error(acbind(a1, b1), expected_msg)      # binary op
    expect_error(acbind(b1, a1), expected_msg)      # binary op
    expect_error(acbind(a1, b1, a1), expected_msg)  # ternary op
})

