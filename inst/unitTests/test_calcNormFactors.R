test_calcNormFactors <- function() {
    data(hypoData)
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", hypoData, group)
    tcc <- calcNormFactors(tcc)
}


