test_plotFCPseudocolor <- function() {
    tcc <- simulateReadCounts(Ngene = 1000)
    plotFCPseudocolor(tcc)

    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3, 3))
    plotFCPseudocolor(tcc)
}

