test_simulateReadCounts <- function() {
    tcc <- simulateReadCounts()

    tcc <- simulateReadCounts(replicates = c(2, 3, 4))
  
    tcc <- simulateReadCounts(Ngene = 1000, PDEG = 0.1, 
                              DEG.assign = c(0.6, 0.4))

    tcc <- simulateReadCounts(Ngene = 1000, PDEG = 0.1, 
                              DEG.assign = c(0.6, 0.4), 
                              DEG.foldchange = c(2, 6))
}

