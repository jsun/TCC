## Check all norm.methods in "test_estimateDE".
test_calcNormFactors_increment <- function() {
    data(hypoData)
    tcc <- new("TCC", hypoData, c(1, 1, 1, 2, 2, 2))
    tcc.0 <- calcNormFactors(tcc, iteration = 0)
    tcc.1 <- calcNormFactors(tcc, iteration = 1)
    tcc.0.1 <- calcNormFactors(tcc, increment = TRUE)
    checkEqualsNumeric(tcc.1$norm.factors, tcc.0.1$norm.factors)

    tcc.3 <- calcNormFactors(tcc, iteration = 3)
    tcc.1 <- calcNormFactors(tcc, increment = TRUE)
    tcc.1.1 <- calcNormFactors(tcc.1, increment = TRUE)
    tcc.1.1.1 <- calcNormFactors(tcc.1.1, increment = TRUE)
    checkEqualsNumeric(tcc.3$norm.factors, tcc.1.1.1$norm.factors)

    tcc.1 <- calcNormFactors(tcc, iteration = 1)
    tcc.1.2 <- calcNormFactors(tcc.1, iteration = 2, increment = TRUE)
    checkEqualsNumeric(tcc.3$norm.factors, tcc.1.2$norm.factors)
}


