test_estimateDE_crossvalidate <- function() {
    set.seed(201501)
    tcc <- new("TCC")
    test.methods <- tcc$private$available$test.method
    norm.methods <- tcc$private$available$norm.method
    package.name <- rownames(test.methods)
    data.type    <- colnames(test.methods)
    for (d in 1:length(data.type)) {
        x <- switch(data.type[d],
            "TGY" = simulateReadCounts(
                        replicates = c(3, 3),
                        Ngene = 1000, PDEG  = 0.20,
                        DEG.assign = c(0.8, 0.2)),
            "MGY" = simulateReadCounts(
                        replicates = c(3, 3, 3),
                        Ngene = 1000, PDEG  = 0.20,
                        DEG.assign = c(0.6, 0.2, 0.2)),
            "TGNP" = simulateReadCounts(
                        group = data.frame(
                        group  = c("normal", "normal", "normal", "normal",
                                   "tumor", "tumor", "tumor", "tumor"),
                        condition = c("A", "B", "C", "D",
                                      "A", "B", "C", "D")
                        ),
                        DEG.foldchange <- data.frame(
                            group_1 = c(1, 1, 1, 1, 4, 4, 4, 4),
                            group_2 = c(4, 4, 4, 4, 1, 1, 1, 1)
                        ),
                        Ngene = 1000, PDEG = 0.1, DEG.assign = c(0.5, 0.5)
                    ),
            "MF" = simulateReadCounts(
                      group = data.frame(
                      group  = c("WT", "WT", "WT", "WT",
                                 "C1", "C1", "C1", "C1", 
                                 "C2", "C2", "C2", "C2"),
                      condition = c("A", "A", "B", "B", "A", "A", 
                                    "B", "B", "A", "A", "B", "B")
                      ),
                      DEG.foldchange = data.frame(
                        A  = c(4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1),
                        B  = c(1, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 1),
                        C  = c(1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4),
                        AB = c(4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1),
                        BC = c(1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4),
                        CA = c(4, 4, 4, 4, 1, 1, 1, 1, 4, 4, 4, 4),
                        TA = c(4, 4, 1, 1, 4, 4, 1, 1, 4, 4, 1, 1),
                        TB = c(1, 1, 4, 4, 1, 1, 4, 4, 1, 1, 4, 4)
                     ),
                     Ngene = 1000, PDEG = 0.1, 
                     DEG.assign = rep(1 / 8, 8))
        )
        for (p in 1:length(package.name)) {
            if (!test.methods[p, d]) next
            for (n in 1:length(norm.methods)) {
                e <- try(x <- calcNormFactors(x, norm.method = norm.methods[n],
                                     test.method = package.name[p],
                                     iteration = T, samplesize = 10))
                if (class(e) == "try-error") {
                    x <- calcNormFactors(x, norm.method = norm.methods[n],
                                     test.method = package.name[p],
                                     iteration = T, samplesize = 10,
                                     method = "blind", sharingMode = "fit-only")
                }
            }
            e <- try(x <- estimateDE(x, test.method = package.name[p], samplesize = 10))
            if (class(e) == "try-error") {
                x <- estimateDE(x, test.method = package.name[p], samplesize = 10,
                                method = "blind", sharingMode = "fit-only")
            }
            print(calcAUCValue(x))
            checkTrue(calcAUCValue(x) > 0.5)
        }
    }       
}




