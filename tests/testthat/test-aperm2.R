test_that("aperm2", {
    expect_error(aperm2(letters), "'a' must be an array")

    a <- array(1:6000, c(1, 40, 1, 50, 3))
    dimnames(a) <- list(NULL,
                        paste0("B", 1:40),
                        NULL,
                        paste0("D", 1:50),
                        paste0("E", 1:3))

    expected_msg <- "'perm' must be an integer vector"
    expect_error(aperm2(a, perm=letters), expected_msg, fixed=TRUE)

    expected_msg <- "'perm' cannot be an empty vector"
    expect_error(aperm2(a, perm=integer(0)), expected_msg, fixed=TRUE)

    expected_msg <- paste0("all non-NA values in 'perm' ",
                           "must be >= 1 and <= 'length(dim(a))'")
    expect_error(aperm2(a, perm=0), expected_msg, fixed=TRUE)
    expect_error(aperm2(a, perm=c(6,NA)), expected_msg, fixed=TRUE)

    expected_msg <- "'perm' cannot contain non-NA duplicates"
    expect_error(aperm2(a, perm=c(4:1,4)), expected_msg, fixed=TRUE)

    expected_msg <- "only dimensions with an extent of 1 can be dropped"
    expect_error(aperm2(a, perm=1:4), expected_msg, fixed=TRUE)

    for (perm in list(1:5, 5:1, c(4,5,2,1,3)))
        expect_identical(aperm2(a, perm), base::aperm(a, perm))

    expected <- S4Arrays:::set_dim(a, dim(a)[-1])
    expected <- S4Arrays:::set_dimnames(expected, dimnames(a)[-1])
    expect_identical(aperm2(a, perm=2:5), expected)
    expected <- base::aperm(expected, perm=4:1)
    expect_identical(aperm2(a, perm=5:2), expected)

    expected <- S4Arrays:::set_dim(a, dim(a)[-3])
    expected <- S4Arrays:::set_dimnames(expected, dimnames(a)[-3])
    expect_identical(aperm2(a, perm=c(1:2,4:5)), expected)
    expected <- base::aperm(expected, perm=c(4,2,1,3))
    expect_identical(aperm2(a, perm=c(5,2,1,4)), expected)

    perm <- c(NA,NA,5,2,NA,1,4,NA)
    expected <- S4Arrays:::set_dim(expected, c(1, 1, 3, 40, 1, 1, 50, 1))
    expected <- S4Arrays:::set_dimnames(expected, dimnames(a)[perm])
    expect_identical(aperm2(a, perm=perm), expected)

    expect_identical(aperm2(a, perm=c(2,4:5)), drop(a))
    expected <- base::aperm(drop(a), perm=c(2:3,1))
    expect_identical(aperm2(a, perm=c(4:5,2)), expected)
})

