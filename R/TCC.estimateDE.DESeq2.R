TCC$methods(.testByDeseq2 = function(fit0 = NULL, design = NULL, 
                                     DESeq2.test = "LRT", paired = NULL,
                                     ...) {





.testByDeseq2.1 = function(design, DESeq2.test) {
    if (is.null(design)) {
        dgroup <- colnames(.self$group)
        design <- formula(eval(parse(text =
                    paste("~", paste(dgroup, collapse = " + "), sep = " ")
                  )))
    }
    suppressMessages(d <- DESeq2::DESeqDataSetFromMatrix(
                                  countData = round(.self$count),
                                  colData = .self$group,
                                  design = design))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    suppressMessages(DESeq2::sizeFactors(d) <- sz)
    suppressMessages(d <- DESeq2::estimateDispersions(d))
    if (DESeq2.test == "Wald") {
        suppressMessages(d <- DESeq2::nbinomWaldTest(d))
    }
    if (DESeq2.test == "LRT") {
        suppressMessages(d <- DESeq2::nbinomLRT(d, reduced = ~ 1))
    }
    res <- results(d)
    private$stat$p.value <<- res$pvalue
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}





.testByDeseq2.3 = function(design, fit0) {
    if (is.null(design)) {
        dgroup <- colnames(.self$group)
        design <- formula(eval(parse(text =
                    paste("~", paste(dgroup, collapse = " + "), sep = " ")
                  )))
    }
    if (is.null(fit0)) {
        fit0 <- formula(~ 1)
    }
    ## check the reduced model (fit0) fotmat.
    ## if it is DESeq format, change it to DESeq2 format.
    formulatxt <- as.character(fit0)
    if (formulatxt[2] == "count") {
        fit0 <- formula(eval(parse(text =
                  paste("~", formulatxt[3], sep = " ")
                )))
    }
    suppressMessages(d <- DESeq2::DESeqDataSetFromMatrix(
                                  countData = round(.self$count),
                                  colData = .self$group,
                                  design = design))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    suppressMessages(DESeq2::sizeFactors(d) <- sz)
    suppressMessages(d <- DESeq2::estimateDispersions(d))
    suppressMessages(d <- DESeq2::nbinomLRT(d, reduced = ~ fit0))
    res <- results(d)
    private$stat$p.value <<- res$pvalue
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}






##
## main process
##
test.approach <- .self$.testApproach(paired = paired)

switch(test.approach,
    "1" = .testByDeseq2.1(design, DESeq2.test),
    "2" = .testByDeseq2.1(design, DESeq2.test),
    "3" = .testByDeseq2.3(design, fit0),
    stop("TCC::ERROR: TCC does not support such identification strategy.")
)

})
