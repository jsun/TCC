TCC$methods(.testByEdger = function(design = NULL, coef = NULL, 
                                    contrast = NULL, dispersion = NULL, ...) {

.testByEdger.1 = function(dispersion = NULL) {
    suppressMessages(d <- edgeR::DGEList(counts = round(count), 
                                         group = group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- norm.factors
    if (min(table(group[, 1])) > 1) {
        suppressMessages(d <- edgeR::estimateCommonDisp(d))
        suppressMessages(d <- edgeR::estimateTagwiseDisp(d))
    }
    if (is.null(dispersion)) {
        suppressMessages(d <- edgeR::exactTest(d))
    } else {
        suppressMessages(d <- edgeR::exactTest(d, dispersion = dispersion))
    }
    if (!is.null(d$table$PValue)) {
        private$stat$p.value <<- d$table$PValue
    } else {
        private$stat$p.value <<- d$table$p.value
    }
    private$stat$rank <<- rank(private$stat$p.value)
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
}

.testByEdger.2 = function(design = NULL, coef = NULL,
                                      contrast = NULL){
    if (is.null(design))
         design <- model.matrix(~ as.factor(.self$group[, 1]))
    if (is.null(coef) && is.null(contrast))
         coef <- 2:length(unique(.self$group[, 1]))
    suppressMessages(d <- edgeR::DGEList(counts = round(.self$count), 
                                         group = .self$group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- .self$norm.factors
    suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
    suppressMessages(fit <- edgeR::glmFit(d, design))
    suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef, 
                                          contrast = contrast))
    s <- topTags(lrt, n = nrow(.self$count))
    s <- s$table[rownames(.self$count), ]
    private$stat$p.value <<- s$PValue
    private$stat$rank <<- rank(.self$private$stat$p.value)
    private$stat$q.value <<- s$FDR
}

.testByEdger.3 = function(design = NULL, coef = NULL,
                                      contrast = NULL){
    if (is.null(design))
        stop("TCC::ERROR: Need the design matrix for GLM.")
    suppressMessages(d <- edgeR::DGEList(counts = round(.self$count), 
                                         group = .self$group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- .self$norm.factors
    suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
    suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
    suppressMessages(fit <- edgeR::glmFit(d, design))
    suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef, 
                                          contrast = contrast))
    s <- topTags(lrt, n = nrow(.self$count))
    s <- s$table[rownames(.self$count), ]
    private$stat$p.value <<- s$PValue
    private$stat$rank <<- rank(.self$private$stat$p.value)
    private$stat$q.value <<- s$FDR
}

ts <- .self$.testStrategy()
if (ts == 1) {
   .testByEdger.1(dispersion = dispersion)
} else if (ts == 2) {
   .testByEdger.2(design = design, coef = coef,
                        contrast = contrast)
} else if (ts == 3) {
   .testByEdger.3(design = design, coef = coef,
                        contrast = contrast)
} else {
       stop()
}

})

