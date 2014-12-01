TCC$methods(.normByTmm = function(x){
## 
## TMM normalization
## The 'calcNormFactors' implemented in edgeR is used for calculating.
## 
    suppressMessages(d <- edgeR::DGEList(counts = round(x), 
                                         group = group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    normf <- d$samples$norm.factors
    names(normf) <- colnames(.self$count)
    return(normf)
})





TCC$methods(.normByDeseq = function(x){
##
## Med(DESeq) normalization
## The 'estimateSizeFactors' implemented in DESeq is used for calculating.
## Since 'estimateSizeFactors' outputs size factors, we convert size factor
## to normalization factors before return factors.
##
    if (ncol(group) == 1) {
        suppressMessages(d <- DESeq::newCountDataSet(countData = round(x), 
                                              conditions = group[, 1]))
    } else {
        suppressMessages(d <- DESeq::newCountDataSet(countData = round(x), 
                                              conditions = group))
    }
    suppressMessages(d <- DESeq::estimateSizeFactors(d))
    return(DESeq::sizeFactors(d) / colSums(x))
})




TCC$methods(.normByDeseq2 = function(x){
##
## Med(DESeq2) normalization
## The 'estimateSizeFactors' implemented in DESeq2 is used for calculating.
## Since 'estimateSizeFactors' outputs size factors, we convert size factor
## to normalization factors before return factors.
##
    design <- formula(as.formula(paste("~", paste(colnames(.self$group), collapse = "+"))))
    suppressMessages(d <- DESeq2::DESeqDataSetFromMatrix(
                                      countData = round(x), 
                                      colData = .self$group,
                                      design = design))
    suppressMessages(d <- DESeq2::estimateSizeFactors(d))
    return(DESeq2::sizeFactors(d) / colSums(x))
})





TCC$methods(calcNormFactors = function(norm.method = NULL,
                                       test.method = NULL,
                                       iteration = 1,
                                       FDR = NULL,
                                       floorPDEG = NULL,
                                       increment = FALSE,
                                       ...) {
## 
## The definition of 'calcNormFactors' function in TCC. This function
## performs the DEGES-based normalization. The method generally performs
## three steps, that is X-Y-X. 'X' is the normalization method which is
## given by 'norm.method', and 'Y' is the identification method which is
## given by 'test.method'. The advantage of this function allows user to
## tune X-Y-X pipeline. For example,
##   -------------------------------------------
##   |iteration   |pipeline                    |
##   ------------------------------------------|
##   |FALSE, 0    |X                           |
##   |TRUE,  1    |X-Y-X                       |
##   |3           |X-Y-X-Y-X-Y-X  => X-(Y-X)3  |
##   -------------------------------------------
## In step Y, the DEGs identified by 'test.method' were removed from
## original data before go to next X step. The removing thresholds are
## 'FDR' and 'floorPDEG'.
## The ellipsis (three-dots) is an additional arguments, they will passed
## to the methods for identifying DEGs in Y step.
##
    ## Start to record execution time.
    ex.time <- proc.time()

    ## Initialize arugments.
    ## If 'norm.method', 'test.method', etc were not given, initiailze them
    ## depend on data type ('tcc$count' and 'tcc$group').

    ## Initialize 'start points' of DEGES normlization.
    if ((increment == FALSE) || 
        (increment == TRUE && .self$private$normalized == FALSE)) {
        DEGES$iteration <<- 0
    }
    ## Initialize 'norm.method' an 'test.method'.
    if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1))  {
        if (is.null(norm.method)) norm.method <- "deseq"
        if (is.null(test.method)) test.method <- "deseq"
    } else {
        if (is.null(norm.method)) norm.method <- "tmm"
        if (is.null(test.method)) test.method <- "edger"
    }
    ## Initialize 'FDR' threshold.
    if (is.null(FDR)) {
        FDR <- 0.1
    }
    ## Initialize 'floorPDEG' threshold.
    if (is.null(floorPDEG)) {
            floorPDEG <- 0.05
    }

    ## Initialize 'iteration'.
    if (iteration) {
        if (is.logical(iteration)) {
            iteration <- 1
        }
        if (iteration < 0 && 100 < iteration) {
            stop("TCC::ERROR: The iteration must be given within the limits of from '0' to '100'.")
        }
        ## Print out DEGES pipeline.
        message(paste("TCC::INFO: Calculating normalization factors using DEGES"))
        message(paste("TCC::INFO: (iDEGES pipeline :", norm.method, 
                      "- [", test.method, "-", norm.method, "] X", 
                      iteration + DEGES$iteration, ")"))
        ## Save DEGES pipeline to TCC class object.
        DEGES$pipeline <<- paste(norm.method, "- [", test.method, 
                                 "-", norm.method, "] X", 
                                 iteration + DEGES$iteration)
    } else {
        message(paste("TCC::INFO: Calculating normalization factors using", norm.method, "..."))
        DEGES$pipeline <<- norm.method
    }

    ## The first X step of X-Y-X pipeline.
    if ((increment == FALSE) || 
        (increment == TRUE && private$normalized == FALSE)) {
        norm.factors <<- switch(norm.method,
                                "tmm" = .self$.normByTmm(count),
                                "deseq" = .self$.normByDeseq(count),
                                "deseq2" = .self$.normByDeseq2(count),
                                stop(paste("\nTCC::ERROR: The normalization method of ", 
                                norm.method, " doesn't supported.\n")))
    }
    norm.factors <<- norm.factors / mean(norm.factors) #standardize norm.factors
    DEGES$threshold <<- data.frame(type = "Unused", input = 0, PDEG = 0)

    ## The Y-X step of X-Y-X pipeline.
    if (iteration) {
        ## Repeat Y-X step depends on 'iteration' argument.
        for (i in 1:iteration) {
            ## The Y step of X-Y-X pipeline.
            ## The DEGs will be identified and remvoed in this step.
            ## Identification process.
            private$stat <<- list()
            switch(test.method,
                   "edger" = .self$.testByEdger(...),
                   "deseq" = .self$.testByDeseq(...),
                   "deseq2" = .self$.testByDeseq2(...),
                   "bayseq" = .self$.testByBayseq(...),
                   "voom" = .self$.testByLimmavoom(...),
                   "samseq" = .self$.testBySamseq(...),
                   "wad" = .self$.testByWad(...),
                   #"yayoi" = .self$.testByYayoi(norm = TRUE, ...),
                   stop(paste("\nTCC::ERROR: The identifying method of ", test.method, " doesn't supported.\n"))
            )
            ## Removing process.
            deg.flg <- rep(0, length = nrow(count))            # potential DEG (be used)
            deg.flg.FDR <- .self$.exactTest(FDR = FDR)         # potential DEG (by FDR, TbT)
            deg.flg.floorPDEG <- rep(0, length = nrow(count))  # potential DEG (by floorPDEG, stats, etc)

            if (is.null(.self$private$stat$testStat) &&
                is.null(.self$private$stat$prob)) {
                deg.flg.floorPDEG <- as.numeric(rank(private$stat$p.value, 
                                 ties.method = "min") <= nrow(count) * floorPDEG)
                if (is.null(FDR)) {
                    ## use TbT
                    deg.flg <- deg.flg.FDR
                    DEGES$threshold$type <<- "TbT"
                    DEGES$threshold$input <<- private$tbt$estProps
                    DEGES$threshold$PDEG <<- sum(deg.flg) / length(deg.flg)
                    private$DEGES.PrePDEG <<- deg.flg
                } else {
                    ## use FDR
                    deg.flg <- deg.flg.FDR
                    DEGES$threshold$type <<- "FDR"
                    DEGES$threshold$input <<- FDR
                    DEGES$threshold$PDEG <<- sum(deg.flg) / length(deg.flg)
                    private$DEGES.PrePDEG <<- deg.flg
                }
            } else if (is.null(.self$private$stat$testStat) &&
                       !is.null(.self$private$stat$prob)) {
                deg.flg.floorPDEG <- as.numeric(rank(- abs(private$stat$prob), 
                                 ties.method = "min") <= nrow(count) * floorPDEG)
                private$DEGES.PrePDEG <<- rep(0, length = nrow(count))
            } else {
                deg.flg.floorPDEG <- as.numeric(rank(- abs(private$stat$testStat), 
                                 ties.method = "min") <= nrow(count) * floorPDEG)
                private$DEGES.PrePDEG <<- rep(0, length = nrow(count))
            }
            ## Check the percentage of the potential DEGs.
            ## If the percentage is less than floorPDEG. redefine potential DEGs.
            if (sum(deg.flg != 0) < sum(deg.flg.floorPDEG != 0)) {
                ## use floorPDEG
                deg.flg <- deg.flg.floorPDEG
                DEGES$threshold$type <<- "floorPDEG"
                DEGES$threshold$input <<- floorPDEG
                DEGES$threshold$PDEG <<- sum(deg.flg) / length(deg.flg)
            }
            count.ndeg <- count[(deg.flg == 0), ]
            if (nrow(count.ndeg) == 0) {
                message ("TCC::INFO: No non-DE genes after eliminate DE genes. stop DEGES strategy.")
                break
            }
            ## The last X of X-Y-X pipeline.
            norm.factors <<- switch(norm.method,
                                    "tmm" = .self$.normByTmm(count.ndeg),
                                    "deseq" = .self$.normByDeseq(count.ndeg),
                                    "deseq2" = .self$.normByDeseq2(count.ndeg)
            )
            ## Standarlize normalization factors.
            norm.factors <<- norm.factors * colSums(count.ndeg) / colSums(count)
            norm.factors <<- norm.factors / mean(norm.factors)
            DEGES$iteration <<- DEGES$iteration + 1
        }
        DEGES$potDEG <<- deg.flg
        DEGES$prePotDEG <<- .self$private$DEGES.PrePDEG
    }
    message("TCC::INFO: Done.")
    DEGES$execution.time <<- proc.time() - ex.time
    private$normalized <<- TRUE
})

