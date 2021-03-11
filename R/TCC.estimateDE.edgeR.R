TCC$methods(.testByEdger = function(
        design   = NULL,
        coef     = NULL, 
        contrast = NULL,
        paired   = NULL,
        ...) {



# Two approach can be performed, i.e., exact test and LRT. If the 'design' is
# given, then performe the LRT, or not performe exact test as default.
.testByEdger.1 = function(design, coef, contrast) {
    suppressMessages(d <- edgeR::DGEList(counts = round(.self$count), 
                                         group  = .self$group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- .self$norm.factors

    if (is.null(design)) {
        # ----------------------
        # The following code is the legacy version of eadgeR for two-group test.
        # We changed the deafult two-group test approach from exactTest to GLM based method. 2020.12.
        ##
        ## if (length(unique(.self$group[, 1])) == nrow(.self$group)) {
        ##     suppressMessages(d <- edgeR::estimateGLMCommonDisp(d,
        ##                           method = "deviance", robust = TRUE,
        ##                           subset = NULL))
        ## } else {
        ##     suppressMessages(d <- edgeR::estimateCommonDisp(d))
        ##     suppressMessages(d <- edgeR::estimateTagwiseDisp(d))
        ## }
        ## suppressMessages(d <- edgeR::exactTest(d))
        ## p <- d$table$PValue
        # ----------------------
        
        design <- model.matrix(~ .self$group[, 1])
        suppressMessages(d <- estimateDisp(d, design))
        suppressMessages(fit <- glmQLFit(d, design))
        suppressMessages(out <- glmQLFTest(fit, coef = 2))
        p <- out$table$PValue
    } else {
        if (length(unique(.self$group[, 1])) == nrow(.self$group)) {
            suppressMessages(d <- edgeR::estimateGLMCommonDisp(d,
                                  method = "deviance", robust = TRUE,
                                  subset = NULL))
        } else {
            if (is.null(coef) && is.null(contrast)) coef <- 2
            suppressMessages(d <- estimateDisp(d, design))
            #suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
            #suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
            #suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
        }
        #suppressMessages(fit <- glmFit(d, design))
        #suppressMessages(lrt <- glmLRT(fit, coef = coef, contrast = contrast))
        #p <- lrt$table$PValue
        suppressMessages(fit <- glmQLFit(d, design))
        suppressMessages(out <- glmQLFTest(fit, coef = coef))
        p <- out$table$PValue
    }
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$rank <<- rank(p)
    private$stat$q.value <<- p.adjust(p, method = "BH")
}





.testByEdger.2 = function(design, coef, contrast) {
    ## ANOVA like design if not given
    if (is.null(design)) design <- model.matrix(~ ., data = .self$group)
    if (is.null(coef) && is.null(contrast)) coef <- 2:ncol(design)
    suppressMessages(d <- edgeR::DGEList(counts = round(.self$count), 
                                         group = .self$group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- .self$norm.factors
    if (length(unique(.self$group[, 1])) == nrow(.self$group)) {
        suppressMessages(d <- edgeR::estimateGLMCommonDisp(d,
                              method = "deviance", robust = TRUE,
                              subset = NULL))
    } else {
        suppressMessages(d <- estimateDisp(d, design))
        #suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
        #suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
        #suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
    }
    suppressMessages(fit <- glmQLFit(d, design))
    suppressMessages(out <- glmQLFTest(fit, coef = coef))
    p <- out$table$PValue
    #suppressMessages(fit <- edgeR::glmFit(d, design))
    #suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef, contrast = contrast))
    #p <- lrt$table$PValue
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$rank <<- rank(p)
    private$stat$q.value <<- p.adjust(p, method = "BH")
}




# If the arguments ('design', 'coef', 'contrast') are given, then used them.
# If not given, then perform ANOVA like LRT.
.testByEdger.3 = function(design, coef, contrast) {
    if (is.null(design)) design <- model.matrix(~ ., data = .self$group)
    if (is.null(coef) && is.null(contrast)) coef <- 2:ncol(design)
    suppressMessages(d <- edgeR::DGEList(counts = round(.self$count), 
                                         group = .self$group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- .self$norm.factors
    #suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
    #suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
    #suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
    #suppressMessages(fit <- edgeR::glmFit(d, design))
    #suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef, contrast = contrast))
    #p <- lrt$table$PValue
    suppressMessages(d <- estimateDisp(d, design))
    suppressMessages(fit <- glmQLFit(d, design))
    suppressMessages(out <- glmQLFTest(fit, coef = 2))
    p <- out$table$PValue
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$rank <<- rank(p)
    private$stat$q.value <<- p.adjust(p, method = "BH")
}





.testByEdger.4 <- function(design, coef, contrast) {
    if (is.null(design)) design <- model.matrix(~ ., data = .self$group)
    if (is.null(coef) && is.null(contrast)) coef <- 2
    suppressMessages(d <- DGEList(counts = round(.self$count),
                                  group = .self$group[, 1]))
    suppressMessages(d <- edgeR::calcNormFactors(d))
    d$samples$norm.factors <- .self$norm.factors
    #suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
    #suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
    #suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
    #suppressMessages(fit <- edgeR::glmFit(d, design))
    #suppressMessages(lrt <- edgeR::glmLRT(fit, coef = coef,  contrast = contrast))
    #p <- lrt$table$PValue
    suppressMessages(d <- estimateDisp(d, design))
    suppressMessages(fit <- glmQLFit(d, design))
    suppressMessages(out <- glmQLFTest(fit, coef = 2))
    p <- out$table$PValue
    p[is.na(p)] <- 1
    private$stat$p.value <<-p
    private$stat$rank <<- rank(p)
    private$stat$q.value <<- p.adjust(p, method = "BH")
}





test.approach <- .self$.testApproach(paired = paired)
switch(test.approach,
    "1" = .testByEdger.1(design = design, coef = coef, contrast = contrast),
    "2" = .testByEdger.2(design = design, coef = coef, contrast = contrast),
    "3" = .testByEdger.3(design = design, coef = coef, contrast = contrast),
    "4" = .testByEdger.4(design = design, coef = coef, contrast = contrast),
    stop("TCC::ERROR: TCC does not support such identification strategy.")
)


})

