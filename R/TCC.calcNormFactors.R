TCC$methods(.normByTmm = function(x){
    suppressMessages(d <- edgeR::DGEList(counts = round(x), 
                                         group = group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    normf <- d$samples$norm.factors
    names(normf) <- colnames(.self$count)
    return(normf)
})

TCC$methods(.normByDeseq = function(x){
    if (ncol(group) == 1) {
        suppressMessages(d <- newCountDataSet(countData = round(x), 
                                              conditions = group[, 1]))
    } else {
        suppressMessages(d <- newCountDataSet(countData = round(x), 
                                              conditions = group))
    }
    suppressMessages(d <- estimateSizeFactors(d))
    return(sizeFactors(d) / colSums(x))
})

#TCC$methods(.normByTwad = function(x, refColumn = NULL, trimWAD, q, AD) {
#    libsize <- colSums(x)
#
#    allzero <- as.logical(rowSums(x) == 0)
#    if (any(allzero))
#      x <- x[!allzero, , drop = FALSE]
#    private$twad.trim <<- gene_id[!allzero]
#
#    ## set reference column
#    y <- sweep(x, 2, 1 / libsize, "*")
#    f75 <- apply(y, 2, function(x) quantile(x, p = 0.75))
#    if (is.null(refColumn)) {
#       refColumn <- which.min(abs(f75 - mean(f75)))
#       if (length(refColumn) == 0 | refColumn < 1 | refColumn > ncol(x))
#           refColumn <- 1
#    }
#    ## norm factors
#    nf <- rep(1, length = ncol(x))
#    for (i in 1:length(nf)) {
#        nf[i] <- .self$.twadcore(obs = x[, i], ref = x[, refColumn],
#                           obs.libsize = libsize[i],
#                           ref.libsize = libsize[refColumn],
#                           trimWAD = trimWAD,
#                           q = q, AD = AD)
#    }
#    nf <- nf / exp(mean(log(nf)))
#    return (nf)
#})

##TCC$methods(.twadcore = function(obs, ref, obs.libsize, ref.libsize, 
##                                 trimWAD, q, AD) {
##    if (all(obs == ref))
##        return (1)
##    ## libsize
##    obs <- as.numeric(obs)
##    ref <- as.numeric(ref)
##    lowcount <- as.logical(obs <= quantile(obs, q) |
##                           ref <= quantile(ref, q))
##    obs <- obs[!lowcount]
##    ref <- ref[!lowcount]
##    private$twad.trim <<- private$twad.trim[!lowcount]
##
##    ## calculate wad
##    wad <- .wad(x = cbind(obs, ref), g = c(1, 2),
##                log.scale = TRUE, AD = AD)
##    rnk <- rank(abs(wad))
##    rnk.sort <- rnk[rev(order(rnk))]
##    min.idx <- min(rnk.sort[1:round(length(obs) * trimWAD)])
##
##    ## calculate normfactors
##    v <- (obs.libsize - obs) / (obs.libsize * obs) +
##         (ref.libsize - ref) / (ref.libsize * ref)
##    v <- v[rnk <= min.idx]
##    obs <- obs[rnk <= min.idx]
##    ref <- ref[rnk <= min.idx]
##    trimmed.geneid <- private$twad.trim[rnk <= min.idx]
##    private$twad.trim <<- rep(0, length = nrow(.self$count))
##    names(private$twad.trim) <<- rownames(.self$count)
##    private$twad.trim[trimmed.geneid] <<- 1
##    nf <- 2^(sum(log2((obs / obs.libsize) / (ref / ref.libsize)) / v,
##             na.rm = TRUE) / (sum(1 / v, na.rm = TRUE)))
##    return (nf)
##})

TCC$methods(calcNormFactors = function(norm.method = NULL,
                                       test.method = NULL,
                                       iteration = 1,
                                       FDR = NULL,
                                       floorPDEG = 0.05,
                                       increment = FALSE,
                                       ...) {
    argus <- list(...)
    if ((increment == FALSE) || 
        (increment == TRUE && private$normalized == FALSE)) {
        DEGES$iteration <<- 0
    }
    ex.time <- proc.time()
    if (is.null(norm.method)) {
        if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
            norm.method = "deseq"
        else 
            norm.method = "edger"
    }
    if (is.null(test.method)) {
        if (!is.null(argus$paired) && argus$paired)
            test.method = "bayseq"
        else if ((ncol(group) == 1) && (min(as.numeric(table(group))) == 1)) 
            test.method = "deseq"
        else 
            test.method = "edger"
    }
    if (norm.method == "edger")
        norm.method <- "tmm" 
    if (test.method != "bayseq" && is.null(FDR))
        FDR <- 0.1
    if (iteration) {
        if (is.logical(iteration))
            iteration <- 1
        if (iteration < 0 && 100 < iteration)
            stop("TCC::ERROR: The iteration must be given within the limits of from '0' to '100'.")
        message(paste("TCC::INFO: Calculating normalization factors using DEGES"))
        message(paste("TCC::INFO: (iDEGES pipeline :", norm.method, 
                      "- [", test.method, "-", norm.method, "] X", 
                      iteration + DEGES$iteration, ")"))
        DEGES$pipeline <<- paste(norm.method, "- [", test.method, 
                                 "-", norm.method, "] X", 
                                 iteration + DEGES$iteration)
    } else {
        message(paste("TCC::INFO: Calculating normalization factors using", norm.method, "..."))
        DEGES$pipeline <<- norm.method
    }
    ## DEGES strategy STEP 1. (First normalization)
    if ((increment == FALSE) || 
        (increment == TRUE && private$normalized == FALSE)) {
        norm.factors <<- switch(norm.method,
                                "tmm" = .self$.normByTmm(count),
                                "deseq" = .self$.normByDeseq(count),
                                ##"twad" = .self$.normByTwad(count, 
                                ##                           trimWAD = trimWAD, 
                                ##                           q = q, AD = AD), 
                                stop(paste("\nTCC::ERROR: The normalization method of ", 
                                norm.method, " doesn't supported.\n")))
    }
    norm.factors <<- norm.factors / mean(norm.factors)
    DEGES$threshold <<- data.frame(type = "Unused", input = 0, PDEG = 0)
    ##  if DEGES not NULL then start DEGES strategy.
    if (iteration) {
        ##  if iteration > 0 then change to iterate DEGES strategy.
        for (i in 1:iteration) {
            ## DEGES strategy  STEP 2. (exact test and remove DEGs.)
            private$stat <<- list()
            switch(test.method,
                   "edger" = .self$.testByEdger(...),#design = design, 
                                                #coef = coef, 
                                                #contrast = contrast, 
                                                #dispersion = dispersion),
                   "deseq" = .self$.testByDeseq(...),#fit1 = fit1,
                                                #fit0 = fit0), 
                   "bayseq" = .self$.testByBayseq(...),#samplesize = samplesize,
                                                  #cl = cl,
                                                  #comparison = comparison),
                   "noiseq" = .self$.testByNoiseq(...),
                   "ebseq" = .self$.testByEbseq(...),
                   "samseq" = .self$.testBySamseq(...),
                   "wad" = .self$.testByWad(...),
                   stop(paste("\nTCC::ERROR: The identifying method of ", test.method, " doesn't supported.\n"))
            )
            ## Remove the DEG from original count data.
            deg.flg <- rep(0, length = nrow(count))
            deg.flg.FDR <- .self$.exactTest(FDR = FDR)
            deg.flg.floorPDEG <- rep(0, length = nrow(count))

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
            ## DEGES strategy STEP 3. (Second normalization)
            norm.factors <<- switch(norm.method,
                                    "tmm" = .self$.normByTmm(count.ndeg),
                                    "deseq" = .self$.normByDeseq(count.ndeg)
                                    ##"twad" = .self$.normByTwad(count.ndeg, 
                                    ##                           trimWAD = trimWAD,
                                    ##                           q = q, AD = AD)
            )
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

