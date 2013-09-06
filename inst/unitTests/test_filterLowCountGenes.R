test_filterLowCountGenes <- function() {
    data(hypoData)
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", hypoData, group)

    tcc <- filterLowCountGenes(tcc)
    
    filter <- as.logical(rowSums(hypoData) > 0)
    hypoData.filtered <- hypoData[filter, ]

    checkEqualsNumeric(tcc$count, hypoData.filtered)
}

