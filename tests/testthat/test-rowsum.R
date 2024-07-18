
.test_fast_colsum <- function(m, group)
{
    fast_colsum <- S4Arrays:::.fast_colsum
    stopifnot(is.matrix(m))
    tm <- t(m)

    check_result <- function(current, expected) {
        #expect_true(is.matrix(current))
        #expect_identical(typeof(current), typeof(expected))
        #if (typeof(expected) == "double") {
        #    expect_equal(current, expected)
        #} else {
            expect_identical(current, expected)
        #}
    }

    current <- fast_colsum(m, group)
    expected <- t(base::rowsum(tm, group))
    check_result(current, expected)

    current <- fast_colsum(m, group, na.rm=TRUE)
    expected <- t(base::rowsum(tm, group, na.rm=TRUE))
    check_result(current, expected)

    current <- fast_colsum(m, group, reorder=FALSE)
    expected <- t(base::rowsum(tm, group, reorder=FALSE))
    check_result(current, expected)

    current <- fast_colsum(m, group, reorder=FALSE, na.rm=TRUE)
    expected <- t(base::rowsum(tm, group, reorder=FALSE, na.rm=TRUE))
    check_result(current, expected)
}

test_that("S4Arrays:::.fast_colsum()", {
    ## type "double"
    m1 <- matrix(0, nrow=4, ncol=6)
    group <- c("B", "A", "B", "B", "B", "A")
    rownames(m1) <- letters[1:4]
    .test_fast_colsum(m1, group)

    m1[1 , ] <- c(8.55, Inf, NA_real_, 0, NaN, -Inf)
    m1[3 , ] <- c(0.6, -11.99, 0, 4.44, 0, 0)
    m1[4 , ] <- 1:6
    .test_fast_colsum(m1, group)
    .test_fast_colsum(m1[0L,   ], group)
    .test_fast_colsum(m1[  , 0L], integer(0))
    .test_fast_colsum(m1[0L, 0L], integer(0))

    ## type "integer"
    m2 <- matrix(0L, nrow=4, ncol=6)
    dimnames(m2) <- list(letters[1:4], letters[21:26])
    m2[2, 1] <- NA_integer_
    m2[2, 3] <- 99L
    m2[4,  ] <- 1:6
    .test_fast_colsum(m2, group)

    ## with NAs in 'group'
    group2 <- c("B", NA, "B", "B", NA, "A")
    suppressWarnings(.test_fast_colsum(m1, group2))
    suppressWarnings(.test_fast_colsum(m2, group2))
})

