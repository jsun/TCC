test_plotFCPseudocolor <- function() {
    tcc <- simulateReadCounts()
    plotFCPseudocolor(tcc)

    tcc <- simulateReadCounts(replicates = c(3, 3, 3))
    plotFCPseudocolor(tcc)
}

