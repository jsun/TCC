TCC$methods(.testStrategy = function () {
    fc <- .self$group
    og <- fc[, 1]
    ug <- unique(og)
    ts <- -1    
    if (ncol(fc) > 1) {
        ##  Multi-factors
        ts <- 3
    } else if (ncol(fc) == 1 & length(ug) > 2) {
        ##  Multi-groups & One-factor
        ts <- 2
    } else if (ncol(fc) == 1 & length(ug) == 2) {
        ##  Two-groups & One-factor
        ts <- 1
    }
    return (ts)
})

TCC$methods(.exactTest = function (FDR = NULL, significance.level = NULL,
                                   PDEG = NULL) {
    deg.flg <- rep(0, length = nrow(count))
    if (!is.null(significance.level)) {
        deg.flg <- as.numeric(private$stat$p.value < significance.level)
    } else if (!is.null(FDR)) {
        deg.flg <- as.numeric(private$stat$q.value < FDR)
    } else if (!is.null(PDEG)) {
        deg.flg <- as.numeric(private$stat$rank <= nrow(count) * PDEG)
    } else {
        deg.flg <- private$estimatedDEG #TbT
    }
    return (deg.flg)
})
 

TCC$methods(estimateDE = function (test.method = NULL,
                                   FDR = NULL, 
#                                   paired = FALSE,
                                   PDEG = NULL,
                                   significance.level = NULL,
                                   dispersion = NULL,
                                   fit0 = NULL, fit1 = NULL,
                                   design = NULL,
                                   contrast = NULL, coef = NULL,
                                   comparison = NULL,
                                   samplesize = NULL,
                                   floor.value = 1,
                                   cl = NULL) {
    paired <- FALSE
    if (is.null(test.method)) {
        if (paired)
            test.method = "bayseq"
        else if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
            test.method = "deseq"
        else 
            test.method = "edger"
    }
    pdeg.method <- c("wad", "noiseq", "samseq")
    if (length(grep(test.method, pdeg.method)) > 0) {
        PDEG <- 0.05
    } else if (test.method != "bayseq" && is.null(FDR) && 
               is.null(significance.level)) {
        FDR <- 0.1
    }
    message(paste("TCC::INFO: Identifying DE genes using", test.method, "..."))
    ## calculate statistics values related DE gene.
    private$stat <<- list()
    stat <<- list()
    switch(test.method,
           "edger" = .self$.testByEdger(design = design, 
                                        coef = coef, 
                                        contrast = contrast,  
                                        dispersion = dispersion,
                                        paired = paired),
           "deseq" = .self$.testByDeseq(fit1 = fit1, 
                                        fit0 = fit0,
                                        paired = paired),
           "bayseq" = .self$.testByBayseq(samplesize = samplesize, 
                                          cl = cl, 
                                          comparison = comparison,
                                          paired = paired),
           "noiseq" = .self$.testByNoiseq(paired = paired),
           "ebseq"  = .self$.testByEbseq(samplesize = samplesize,
                                         paired = paired),
           "samseq" = .self$.testBySamseq(samplesize = samplesize,
                                         paired = paired),
           ##"nbpseq"  = .self$.testByNbpseq(),
           "wad" = .self$.testByWad(floor.value = floor.value),
           stop(paste("\nTCC::ERROR: The identifying method of ", 
                      test.method, " doesn't supported.\n"))
    )
    ## identify DE genes with the results of exact test.
    estimatedDEG <<- .self$.exactTest(FDR = FDR, 
                                      significance.level = significance.level,
                                      PDEG = PDEG)
    if (!is.null(private$stat$testStat))
        stat$testStat <<- private$stat$testStat
    if (!is.null(private$stat$prob))
        stat$prob <<- private$stat$prob
    if (!is.null(private$stat$likelihood))
        stat$likelihood <<- private$stat$likelihood
    if (!is.null(private$stat$p.value))
        stat$p.value <<- private$stat$p.value
    if (!is.null(private$stat$q.value))
        stat$q.value <<- private$stat$q.value
    if (!is.null(private$stat$rank)) 
        stat$rank <<- private$stat$rank
    private$estimated <<- TRUE
    message("TCC::INFO: Done.")
})

