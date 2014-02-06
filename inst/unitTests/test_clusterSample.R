test_clusterSample <- function() {
    data(hypoData)
    tcc <- new("TCC", hypoData, c(1, 1, 1, 2, 2, 2))

    h1 <- clusterSample(hypoData)
    h2 <- clusterSample(unique(hypoData[rowSums(hypoData) > 0, ]))
    h3 <- clusterSample(tcc)

    checkEquals(h1, h2)
    checkEquals(h1, h3)
}


