.wad <- function(x, group, log.scale = TRUE, floor.value = 1) {
    AD <- FALSE
    if (log.scale) {
        x[x < floor.value] <- floor.value
        x <- log2(x)
    }
    ug <- unique(group)
    s <- combn(length(ug), 2, function(ij, x = x, g = group, ug = ug, AD = AD) {
        g1 <- (g == ug[ij[1]])
        g2 <- (g == ug[ij[2]])
        m1 <- rowMeans(as.matrix(x[, g1]))
        m2 <- rowMeans(as.matrix(x[, g2]))
        if (AD) {
            x_ave <- abs(m1 - m2) / 2
        } else {
            x_ave <- (m1 + m2) / 2
        }
        weight <- (x_ave - min(x_ave)) / (max(x_ave) - min(x_ave))
        s <- (m2 - m1) * weight
        return (s)
    }, TRUE, x, group, ug, AD)
    s <- apply(s, 1, function(i) {
             return (i[max(abs(i)) == abs(i)])
    })
    return(s)
}

WAD <- function(data, group, log.scale = FALSE, floor.value = 1, sort = FALSE) {
    data <- as.matrix(data)
    wad <- .wad(x = data, group = group,
                log.scale = log.scale,
                floor.value = floor.value)
    wad <- data.frame(wad = wad,
                      rank = rank(- abs(wad)))
    if (!is.null(rownames(data))) {
        rownames(wad) <- rownames(data)
    } else {
        rownames(wad) <- 1:nrow(data)
    }
    if (sort) {
        wad <- wad[order(wad[, 2]), ]
    }
    return(wad)
}

