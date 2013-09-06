## Kadota original code
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
    g <- c(1, 1, 2, 2)
    h <- c("A", "A", "B", "B")
    x <- matrix(rnorm(100, 10, 2), ncol = 4)
    ef <- colSums(x)
    x <- sweep(x, 2, mean(ef) / ef, "*")
    y <- x
    y[y < 1] <- 1
    y <- log2(y)

    kdt <- kadota_WAD(y, g)
    wad.x <- WAD(x, g, log.scale = TRUE, floor.value = 1)
    wad.y <- WAD(y, g)

    checkEqualsNumeric(as.matrix(kdt), as.matrix(wad.x[, 1]))
    checkEqualsNumeric(as.matrix(kdt), as.matrix(wad.y[, 1]))

    wad.h <- WAD(y, h)
    checkEqualsNumeric(as.matrix(kdt), as.matrix(wad.h[, 1]))


    tcc.g <- new("TCC", x, g)
    tcc.g <- estimateDE(tcc.g, test.method = "wad")
    tcc.h <- new("TCC", x, h)
    tcc.h <- estimateDE(tcc.h, test.method = "wad")
    checkEqualsNumeric(kdt, tcc.g$stat$testStat)
    checkEqualsNumeric(kdt, tcc.h$stat$testStat)
}


