test_that("subassign an SVT_SparseArray object by an Mindex or Lindex", {
    ## Only zeros.
    a0 <- array(0L, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    svt0 <- as(a0, "SVT_SparseArray")
    Mindex3 <- rbind(c(7,  9, 3), c(7, 10, 3), c(6, 4, 3), c(2, 4, 3),
                     c(1, 10, 3), c(7, 10, 3), c(1, 1, 3), c(5, 4, 3),
                     c(2,  4, 3))
    Lindex3 <- Mindex2Lindex(Mindex3, dim(a0))
    vals <- c(11:18, 0L)
    a <- `[<-`(a0, Mindex3, value=vals)
    svt <- `[<-`(svt0, Mindex3, value=vals)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    svt <- `[<-`(svt0, Lindex3, value=vals)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    svt <- `[<-`(svt0, as.double(Lindex3), value=vals)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    m0 <- a0[ , , 1]  # 2D
    svt0 <- as(m0, "SVT_SparseMatrix")
    Mindex2 <- Mindex3[ , -3]
    Lindex2 <- Mindex2Lindex(Mindex2, dim(m0))
    m <- `[<-`(m0, Mindex2, value=vals)
    svt <- `[<-`(svt0, Mindex2, value=vals)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    svt <- `[<-`(svt0, Lindex2, value=vals)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    svt <- `[<-`(svt0, as.double(Lindex2), value=vals)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)

    x0 <- as.array(m0[1, ])  # 1D
    svt0 <- as(x0, "SVT_SparseArray")
    Mindex1 <- Mindex2[ , -2, drop=FALSE]
    Lindex1 <- Mindex2Lindex(Mindex1, dim(x0))
    a <- `[<-`(a0, Mindex1, value=vals)
    #svt <- `[<-`(svt0, Mindex1, value=vals)  # not ready yet!
    #svt <- `[<-`(svt0, Lindex1, value=vals)  # not ready yet!
    #svt <- `[<-`(svt0, as.double(Lindex1), value=vals)  # not ready yet!

    ## Add some nonzero values.
    a0 <- make_3D_test_array()
    svt0 <- as(a0, "SVT_SparseArray")
    Mindex <- rbind(cbind(Mindex2, 1L), Mindex3)
    Lindex <- Mindex2Lindex(Mindex, dim(a0))
    vals2 <- c(vals, vals)
    a <- `[<-`(a0, Mindex, value=vals2)
    svt <- `[<-`(svt0, Mindex, value=vals2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    svt <- `[<-`(svt0, Lindex, value=vals2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    svt <- `[<-`(svt0, as.double(Lindex), value=vals2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    m0 <- a0[ , , 1]  # 2D
    svt0 <- as(m0, "SVT_SparseMatrix")
    m <- `[<-`(m0, Mindex2, value=vals)
    svt <- `[<-`(svt0, Mindex2, value=vals)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    svt <- `[<-`(svt0, Lindex2, value=vals)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    svt <- `[<-`(svt0, as.double(Lindex2), value=vals)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)

    x0 <- as.array(m0[1, ])  # 1D
    svt0 <- as(x0, "SVT_SparseArray")
    a <- `[<-`(a0, Mindex1, value=vals)
    #svt <- `[<-`(svt0, Mindex1, value=vals)  # not ready yet!
    #svt <- `[<-`(svt0, Lindex1, value=vals)  # not ready yet!
    #svt <- `[<-`(svt0, as.double(Lindex1), value=vals)  # not ready yet!

    ## With type change.
    a0 <- make_3D_test_array()
    svt0 <- as(a0, "SVT_SparseArray")
    vals2 <- complex(real=vals2, imaginary=-0.75)
    a <- `[<-`(a0, Mindex, value=vals2)
    svt <- `[<-`(svt0, Mindex, value=vals2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    svt <- `[<-`(svt0, Lindex, value=vals2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    svt <- `[<-`(svt0, as.double(Lindex), value=vals2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
})

