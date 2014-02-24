## Kadota original code.
kadota_WAD <- function(data = NULL, data.cl = NULL){ 
    x <- data
    cl <- data.cl
    mean1 <- rowMeans(as.matrix(x[, cl == 1]))
    mean2 <- rowMeans(as.matrix(x[, cl == 2]))
    x_ave <- (mean1 + mean2) / 2
    weight <- (x_ave - min(x_ave)) / (max(x_ave) - min(x_ave))
    statistic <- (mean2 - mean1) * weight
    return(statistic)
}

test_WAD_value <- function() {
    x <- matrix(rnorm(100, 10, 2), ncol = 4)
    g.num <- c(1, 1, 2, 2)
    g.str <- c("A", "A", "B", "B")
    ## y is the logged x.
    y <- x
    y[y < 1] <- 1
    y <- log2(y)
    ## execute WAD test.
    kdt <- kadota_WAD(y, g.num)
    wad.x <- WAD(x, g.num)
    wad.y <- WAD(y, g.num, logged = TRUE)
    wad.z <- WAD(y, g.str, logged = TRUE)
    ## check the numeric.
    checkEqualsNumeric(as.matrix(kdt), as.matrix(wad.x[, 1]))
    checkEqualsNumeric(as.matrix(kdt), as.matrix(wad.y[, 1]))
    checkEqualsNumeric(as.matrix(kdt), as.matrix(wad.z[, 1]))
    ## check the WAD in TCC class.
    tcc.num <- new("TCC", x, g.num)
    tcc.num <- estimateDE(tcc.num, test.method = "wad")
    tcc.str <- new("TCC", y, g.str)
    tcc.str <- estimateDE(tcc.str, test.method = "wad", logged = TRUE)
    checkEqualsNumeric(kdt, tcc.num$stat$testStat)
    checkEqualsNumeric(kdt, tcc.str$stat$testStat)
}
