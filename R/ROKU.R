.outval <- function(y, upper.limit) {
    if (all(is.na(y)))
        y <- rep(0, length = length(y))
    z <- y[!is.na(y)]
    
    N <- length(z)
    fN <- floor(N * upper.limit) + 1

    z.sorted <- sort(z)
    z.ordered <- order(z)
    
    df <- matrix(c(rep(1:fN, times = c(fN:1)), j = sequence(fN:1)),
                 ncol = 2, byrow = FALSE)
    n <- N - df[, 2] - df[, 1] + 2
    s <- N - n
    f <- rep(0, length = N)

    if (sd(z) != 0) {
        ssd <- apply(df, 1, function(d, w = z.sorted, N = N) {
           return(sd(w[d[1]:(N + 1 - d[2])]))
        }, z.sorted, N)
        u <- n * log(ssd * sqrt((n - 1) / n)) + 
             sqrt(2) * s * lfactorial(n) / n
        min.u <- min(u)
        d <- t(df[u == min.u, ])[1:2]
        if (d[1] > 1)
            f[z.ordered[1:(d[1] - 1)]] <- -1
        if (d[2] > 1)
            f[z.ordered[(N + 1 - d[2] + 1):N]] <- 1
    } 
    return(replace(y, !is.na(y), f))
}


.tbw <- function(y) {
    y <- y[!is.na(y)]
    y.m <- median(y)
    y.u <- (y - y.m) / (5 * median(abs(y - y.m)) + 1e-04)
    y.w <- rep(0, length(y))
    y.i <- abs(y.u) <= 1
    y.w[y.i] <- ((1 - y.u^2)^2)[y.i]
    y.b <- sum(y.w * y) / sum(y.w)
}

.entvalmod <- function(y) {
    y <- y[!is.na(y)]
    l <- length(y)
    y <- y[y != 0]
    if (is.na(sd(y))) {
        return (0)
    } else if (sum(y) <= 0 || sd(y) == 0) {
        return (log2(l))
    } else {
        y.m <- median(y)
        y.u <- (y - y.m) / (5 * median(abs(y - y.m)) + 1e-04)
        y.w <- rep(0, length(y))
        y.i <- abs(y.u) <= 1
        y.w[y.i] <- ((1 - y.u^2)^2)[y.i]
        y.b <- sum(y.w * y) / sum(y.w)
        p <- abs(y - y.b)
        p <- p / sum(p)
        e <- - sum(p * log2(p))
        if (is.na(e))
            e <- 0
        return (e)
    }
}

.entval <- function(y) {
    y <- y[!is.na(y)]
    l <- length(y)
    y <- y[y != 0]
    if (is.na(sd(y))) {
        return (0)
    } else if (sum(y) <= 0 || sd(y) == 0) {
        return (log2(l))
    } else {
        p <- y / sum(y)
        e <- - sum(p * log2(p))
        if (is.na(e))
            e <- 0
        return (e)
    }
}

ROKU <- function(data, upper.limit = 0.25, sort = FALSE) {
    rs <- NULL
    if (is.vector(data)) {
        data <- t(matrix(data))
    } else {
        data <- as.matrix(data)
    }
    rs$outliers <- t(apply(t(scale(t(data))), 1,
                          function (y, upper.limit = upper.limit) {
                          .outval(y, upper.limit = upper.limit)
                     }, upper.limit))
    rs$H <- apply(data, 1, .entval)
    rs$modH <- apply(data, 1, .entvalmod)
    rs$rank <- rank(rs$modH)
    rs$Tbw <- apply(data, 1, .tbw)
    if (!is.null(colnames(data))) {
        l <- colnames(data)
    } else {
        l <- paste("tissue", 1:ncol(data), sep = "_")
    }
    if (!is.null(rownames(data))) {
        r <- rownames(data)
    } else {
        r <- 1:nrow(data)
    }
    colnames(rs$outliers) <- l
    rownames(rs$outliers) <- r
    if (sort) {
        reindex <- order(rs$rank)
        rs$outliers <- rs$outliers[reindex, ]
        rs$H <- rs$H[reindex]
        rs$modH <- rs$modH[reindex]
        rs$rank <- rs$rank[reindex]
        rs$Tbw <- rs$Tbw[reindex]
    }
    return (rs)
}

