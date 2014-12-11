TCC$methods(.testApproach = function (paired = NULL) {
##
##  Check model structures for deciding the default identification approach.
##  -------------------------------------------------------------------------
##                         i:ncol(group)   j:unique(group[, 1])   approach
##  -------------------------------------------------------------------------
##  two groups                      = 1                   = 2          1
##  two groups (paired)             = 2                   = 2          4
##  multi groups                    = 1                   > 2          2
##  multi factors                  >= 2                  >= 2          3
##  -------------------------------------------------------------------------
##
    test.approach <- -1
    if (is.null(paired)) {
        paired <- FALSE
    }
    i <- ncol(.self$group)
    j <- length(levels(.self$group[, 1]))
    if ((test.approach == -1) && (i == 1 && j == 2)) {
        test.approach <- 1
    }
    if ((test.approach == -1) && ((i == 2 && j == 2) && paired)) {
        test.approach <- 4
    }
    if ((test.approach == -1) && (i == 1 && j > 2)) {
        test.approach <- 2
    }
    if ((test.approach == -1) && (i >= 1 && j >= 2)) {
        test.approach <- 3
    }
    if (test.approach == -1) {
        stop("TCC::ERROR: TCC cannot decide the test approach. Please check your 'group' data.frame. If you cannot resolve it, please conduct the developer [wukong@bi.a.u-tokyo.ac.jp].")
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
                                   full = NULL, reduced= NULL,   # DESeq(GLM)
                                   design = NULL,                # edgeR(GLM)
                                   contrast = NULL, coef = NULL, # edgeR(GLM)
                                   group = NULL,            # baySeq('group')
                                   samplesize = NULL,            # baySeq, SAMseq
                                   logged = NULL,                # WAD
                                   floor = NULL,                 # WAD
                                   cl = NULL, ...) {             # baySeq
##
## Indentify DEGs from data with given method, i.e., 'test.method'.
##
    fcall <- as.list(match.call(expand.dots = TRUE))
    ## Initialize 'test.method'.
    if (is.null(test.method)) {
        if ((ncol(.self$group) == 1) && (min(as.numeric(table(.self$group))) == 1)) {
            test.method = "deseq"
        } else {
            test.method = "edger"
        }
    }
    ## Initialize threshold for identifying DEGs.
    if (test.method == "wad") {
        PDEG <- 0.05
    } else if (is.null(FDR) && is.null(significance.level)) {
        FDR <- 0.1
    }
    message(paste("TCC::INFO: Identifying DE genes using", test.method, "..."))
    ## Calculate statistics.
    private$stat <<- list()
    stat <<- list()
    bayseq.group <- eval(parse(text = fcall$group))
    switch(test.method,
           "edger" = .self$.testByEdger(design   = design, 
                                        coef     = coef, 
                                        contrast = contrast,  
                                        paired   = paired,
                                        ...),
           "deseq" = .self$.testByDeseq(full    = full, 
                                        reduced = reduced,
                                        paired  = paired,
                                        ...),
           "deseq2" = .self$.testByDeseq2(design   = design, 
                                          full     = full,
                                          reduced  = reduced,
                                          contrast = contrast,
                                          paired   = paired,
                                          ...),
           "bayseq" = .self$.testByBayseq(samplesize = samplesize, 
                                          cl         = cl, 
                                          bgroup     = bayseq.group,
                                          paired     = paired,
                                          ...),
           "samseq" = .self$.testBySamseq(samplesize = samplesize,
                                          paired     = paired,
                                          ...),
           "voom" = .self$.testByLimmavoom(design = design, 
                                           coef   = coef, 
                                           paired = paired,
                                           ...),
           "wad" = .self$.testByWad(logged = logged,
                                    floor = floor,
                                    ...),
           "yayoi" = .self$.testByYayoi(...),
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

