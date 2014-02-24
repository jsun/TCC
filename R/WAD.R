.wad <- function(x, group, logged = FALSE, floor = 1) {
    ## If data does not logged, transfromt it with the base 2 logarithm.
    if (!logged) {
        x[x < floor] <- floor
        x <- log2(x)
    }
    ## WAD only can apply to two-group comparison data. If data consists
    ## of more than three groups, WAD will be performed for all two-group
    ## combinations of given all groups.
    ug <- unique(group)
    s <- combn(length(ug), 2, function(ij, x = x, g = group, ug = ug) {
        g1 <- (g == ug[ij[1]])
        g2 <- (g == ug[ij[2]])
        m1 <- rowMeans(as.matrix(x[, g1]))
        m2 <- rowMeans(as.matrix(x[, g2]))
        x_ave <- (m1 + m2) / 2
        weight <- (x_ave - min(x_ave)) / (max(x_ave) - min(x_ave))
        s <- (m2 - m1) * weight
        return (s)
    }, TRUE, x, group, ug)
    ## Return the maximum value of the absolute wad score, if data consists
    ## of more than three groups.
    s <- apply(s, 1, function(i) {
             return (i[max(abs(i)) == abs(i)])
    })
    return(s)
}

WAD <- function(data, group, logged = FALSE, floor = 1, sort = FALSE) {
## 
## WAD is the function for identifying DEGs from two-group comaprison data of 
## microarray data. If the data dose not logged, set 'logged = FALSE' for
## logging it internally of WAD function.
##
    data <- as.matrix(data)
    wad <- .wad(x = data, group = group, logged = logged, floor = floor)
    wad <- data.frame(wad = wad, rank = rank(- abs(wad)))
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

