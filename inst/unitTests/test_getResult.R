test_getResult <- function() {
    data(hypoData)
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", hypoData, group)
    tcc <- calcNormFactors(tcc)
    tcc <- estimateDE(tcc, test.method = "edger")
    result <- getResult(tcc)
    ma <- plot(tcc)

    checkEqualsNumeric(result$p.value, tcc$stat$p.value)
    checkEqualsNumeric(result$q.value, tcc$stat$q.value)
    checkEqualsNumeric(result$rank, tcc$stat$rank)
    checkEqualsNumeric(result$a.value, ma$a.value)
    checkEqualsNumeric(result$m.value, ma$m.value)
}

