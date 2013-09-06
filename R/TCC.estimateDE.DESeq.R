TCC$methods(.testByDeseq = function(fit1 = NULL, fit0 = NULL, ...) {

.testByDeseq.1 = function() {
    suppressMessages(d <- newCountDataSet(countData = round(.self$count), 
                                          conditions = .self$group[, 1]))
    sizeFactors(d) <- .self$norm.factors * colSums(.self$count)
    if (min(as.numeric(table(.self$group[, 1]))) == 1) {
        ## Single replicates
        e <- try(suppressMessages(d <- estimateDispersions(d, 
                     method = "blind", sharingMode = "fit-only")), 
                 silent = TRUE)
        if (class(e) == "try-error") {
             message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
             message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
             suppressMessages(d <- estimateDispersions(d, 
                      method = "blind", sharingMode = "fit-only",
                      fitType = "local"))
        }
    } else {
        ## Multiple replicates
        e <- try(suppressMessages(d <- estimateDispersions(d)), silent = TRUE)
        if (class(e) == "try-error") {
            message("TCC::WARN: 'estimateDispersions' with method=\"pooled\" in DESeq could not be performed.")
            message("TCC::WARN: 'estimateDispersions' with method=\"blind\" in DESeq was used instead.")
            e <- try(suppressMessages(d <- estimateDispersions(d, 
                         method = "blind", sharingMode = "fit-only")), 
                     silent = TRUE)
            ## try local mode if defaul occurs error
            if (class(e) == "try-error") {
                message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
                message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
                suppressMessages(d <- estimateDispersions(d, 
                         method = "blind", sharingMode = "fit-only", 
                         fitType = "local"))
            }
        }
    }
    ug <- unique(.self$group[, 1])
    suppressMessages(d <- nbinomTest(d, ug[1], ug[2]))
    d$pval[is.na(d$pval)] <- 1
    d$padj[is.na(d$padj)] <- 1
    private$stat$p.value <<- d$pval
    private$stat$q.value <<- d$padj
    private$stat$rank <<- rank(d$pval)
}
.testByDeseq.2 = function(fit1 = NULL, fit0 = NULL) {
    suppressMessages(d <- newCountDataSet(countData = round(count), 
                                              conditions = group[, 1]))
    sizeFactors(d) <- .self$norm.factors * colSums(.self$count)
    if (min(as.numeric(table(.self$group[, 1]))) == 1) {
        ## single replicates
        e <- try(suppressMessages(d <- estimateDispersions(d, 
                      method = "blind", sharingMode = "fit-only")), 
                 silent = TRUE)
        if (class(e) == "try-error") {
             message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
             message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
             suppressMessages(d <- estimateDispersions(d, 
                      method = "blind", sharingMode = "fit-only", 
                      fitType = "local"))
        }
    } else { 
        ## Multiple replicates
        e <- try(suppressMessages(d <- estimateDispersions(d)), silent = TRUE)
        if (class(e) == "try-error") {
            message("TCC::WARN: 'estimateDispersions' with method=\"pooled\" in DESeq could not be performed.")
            message("TCC::WARN: 'estimateDispersions' with method=\"blind\" in DESeq was used instead.")
            e <- try(suppressMessages(d <- estimateDispersions(d, 
                        method = "blind", sharingMode = "fit-only")), 
                     silent = TRUE)
            ## try local mode if detaul occurs error
            if (class(e) == "try-error") {
                message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
                message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
                suppressMessages(d <- estimateDispersions(d, 
                         method = "blind", sharingMode = "fit-only", 
                         fitType = "local"))
            }
        }
    }
    ## GLM for multiple group comparison.
    if (is.null(fit1) && is.null(fit0)) {
        fit1 <- count ~ condition
        fit0 <- count ~ 1
    }
    if (is.null(fit0))
        stop("TCC::ERROR: Need the formula('fit0') to create reduced model regresses for GLM.")
    if (is.null(fit1))
        stop("TCC::ERROR: Need the formula('fit1') to create full model regresses for GLM.")
    capture.output(f0 <- fitNbinomGLMs(d, fit0))
    capture.output(f1 <- fitNbinomGLMs(d, fit1))
    private$stat$p.value <<- nbinomGLMTest(f1, f0)
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}

.testByDeseq.3 = function(fit1 = NULL, fit0 = NULL) {
    suppressMessages(d <- newCountDataSet(countData = round(.self$count), 
                                              conditions = .self$group))
    sizeFactors(d) <- .self$norm.factors * colSums(.self$count)
    ## try default
    e <- try(suppressMessages(d <- estimateDispersions(d)), silent = TRUE)
    ## try blind method
    if (class(e) == "try-error") {
        message("TCC::WARN: 'estimateDispersions' with method=\"pooled\" in DESeq could not be performed.")
        message("TCC::WARN: 'estimateDispersions' with method=\"blind\" in DESeq was used instead.")
        e <- try(suppressMessages(d <- estimateDispersions(d, 
                                  method = "blind", sharingMode = "fit-only")), 
                     silent = TRUE)
        ## try local mode
        if (class(e) == "try-error") {
            message("TCC::WARN: 'estimateDispersions' with sharingMode=\"fit-only\" in DESeq could not be performed.")
            message("TCC::WARN: 'estimateDispersions' with sharingMode=\"local\" in DESeq was used instead.")
            suppressMessages(d <- estimateDispersions(d, 
                             method = "blind", sharingMode = "fit-only", 
                             fitType = "local"))
        }
    }
    if (is.null(fit0))
        stop("TCC::ERROR: Need the formula('fit0') to create reduced model regresses for GLM.")
    if (is.null(fit1))
        stop("TCC::ERROR: Need the formula('fit1') to create full model regresses for GLM.")
    capture.output(f0 <- fitNbinomGLMs(d, fit0))
    capture.output(f1 <- fitNbinomGLMs(d, fit1))
    private$stat$p.value <<- nbinomGLMTest(f1, f0)
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}

ts <- .self$.testStrategy()
if (ts == 1) {
    .testByDeseq.1()
} else if (ts == 2) {
    .testByDeseq.2(fit1 = fit1, fit0 = fit0)
} else if (ts == 3) {
    .testByDeseq.3(fit1 = fit1, fit0 = fit0)
} else {
    stop()
}
})

