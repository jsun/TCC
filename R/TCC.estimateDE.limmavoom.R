TCC$methods(.testByLimmavoom = function(
        design   = NULL,
        coef     = NULL, 
        paired   = NULL,
        ...) {



# Two approach can be performed, i.e., exact test and LRT. If the 'design' is
# given, then performe the LRT, or not performe exact test as default.
.testByLimmavoom.1 = function(design, coef) {
    ## ANOVA like design if not given
    if (is.null(design)) design <- model.matrix(~ ., data = .self$group)
    if (is.null(coef)) coef <- 2:ncol(design)
    v <- limma::voom(.self$count, design, plot = FALSE,
              lib.size = colSums(.self$count) * .self$norm.factors)
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit)
    res <- limma::topTable(fit, coef = coef, n = nrow(count), sort.by = "none")
    p <- res$P.Value
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$rank <<- rank(p)
    private$stat$q.value <<- p.adjust(p, method = "BH")
}




.testByLimmavoom.4 <- function(design, coef) {
    if (is.null(design)) design <- model.matrix(~ ., data = .self$group)
    if (is.null(coef)) coef <- 2
    v <- limma::voom(.self$count, design, plot = FALSE,
              lib.size = colSums(.self$count) * .self$norm.factors)
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit)
    res <- limma::topTable(fit, coef = coef, n = nrow(count), sort.by = "none")
    p <- res$P.Value
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$rank <<- rank(p)
    private$stat$q.value <<- p.adjust(p, method = "BH")
}





test.approach <- .self$.testApproach(paired = paired)
switch(test.approach,
    "1" = .testByLimmavoom.1(design = design, coef = coef),
    "2" = .testByLimmavoom.1(design = design, coef = coef),
    "3" = .testByLimmavoom.1(design = design, coef = coef),
    "4" = .testByLimmavoom.4(design = design, coef = coef),
    stop("TCC::ERROR: TCC does not support such identification strategy.")
)


})

