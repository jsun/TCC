test_plot <- function() {
    data(hypoData)
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", hypoData, group)
    plot(tcc)

    tcc <- calcNormFactors(tcc, iteration = 0)
    plot(tcc)

    tcc <- estimateDE(tcc, test.method = "edger")
    plot(tcc)

    group <- c("A", "A", "A", "B", "B", "B")
    tcc <- new("TCC", hypoData, group)
    plot(tcc)
}

