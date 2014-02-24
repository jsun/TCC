TCC$methods(.testApproach = function (paired = NULL) {
##
##  This method decides the default test approach for data.
##  The decision information is the given arguments ('paired' and 'group').
##  Returning values:
##     1  :  two-group       & one-factor
##     2  :  multi-group     & one-factor
##     3  :                    multi-factors
##     4  :  two,multi-group & two-factor (paired)
##
    if (is.null(paired)) {
        paired <- FALSE
    }
    sample.design <- .self$group
    sample.condition <- sample.design[, 1]
    unique.condition <- unique(sample.condition)
 
    ## default to 1   
    test.approach <- 1
    if (paired) {
      test.approach <- 4
    } else {
      if (ncol(sample.design) > 1) {
          test.approach <- 3
      } 
      if (ncol(sample.design) == 1) {
          if (length(unique.condition) == 2) {
              test.approach <- 1
          } else {
              test.approach <- 2
          }
      }
    }

    return (test.approach)
})





TCC$methods(.exactTest = function (FDR = NULL, significance.level = NULL,
                                   PDEG = NULL) {
##
## Return the DEGs label with the given threshold.
## The order of priority of threhold is:
##    significance.level > FDR > PDEG > others (TbT)
## 'significance.level' is the threshold for p-value.
## 'FDR' is the threshold for q-value.
## 'PDEG' is the threhold for others. Some methods only outputs the
## probabiligy for DEGs, statistics, or somethings. Those methods
## cannot be able to use p-/q-value, thus we choose top-xx% ranked
## genes as DEGs.
## Because (Kadota original) TbT gives the DEGs label, we used it
## directly.
##
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
                                   paired = NULL,
                                   PDEG = NULL,
                                   significance.level = NULL,
                                   dispersion = NULL,            # edgeR(exactTest)
                                   fit0 = NULL, fit1 = NULL,     # DESeq(GLM)
                                   design = NULL,                # edgeR(GLM)
                                   contrast = NULL, coef = NULL, # edgeR(GLM)
                                   comparison = NULL,            # baySeq('group')
                                   samplesize = NULL,            # baySeq, SAMseq
                                   logged = NULL,                # WAD
                                   floor = NULL,                 # WAD
                                   cl = NULL) {                  # baySeq
##
## Indentify DEGs from data with given method, i.e., 'test.method'.
##
    ## Initialize 'test.method'.
    if (is.null(test.method)) {
        if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) {
            test.method = "deseq"
        } else {
            test.method = "edger"
        }
    }
    ## Initialize threshold for identifying DEGs.
    pdeg.method <- c("wad", "samseq")
    if (length(grep(test.method, pdeg.method)) > 0) {
        PDEG <- 0.05
    } else if (test.method != "bayseq" && is.null(FDR) && 
               is.null(significance.level)) {
        FDR <- 0.1
    }
    message(paste("TCC::INFO: Identifying DE genes using", test.method, "..."))
    ## Calculate statistics.
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
           "samseq" = .self$.testBySamseq(samplesize = samplesize,
                                         paired = paired),
           "wad" = .self$.testByWad(logged = logged,
                                    floor = floor),
           stop(paste("\nTCC::ERROR: The identifying method of ", 
                      test.method, " doesn't supported.\n"))
    )
    ## Identify DEGs with the threshold.
    estimatedDEG <<- .self$.exactTest(FDR = FDR, 
                                      significance.level = significance.level,
                                      PDEG = PDEG)
    ## Set up the statistics to TCC class obejct.
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

