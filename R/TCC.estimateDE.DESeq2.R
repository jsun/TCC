TCC$methods(.testByDeseq2 = function(design = NULL, full = NULL, reduced = NULL,
                                     contrast = NULL, paired = NULL,
                                     ...) {




# For two-group comparison
# If 'full' or 'reduced' are given, then use 'nbinomLRT' of DESeq2. If 
# 'contrast' is given, then use 'nbinomWaldTest' of DESeq2. If the three arguments
# are 'NULL', then use 'nbinomWaldTest' of DESeq2.
.testByDeseq2.1 = function(design, full, reduced, contrast) {
    suppressMessages(d <- DESeq2::DESeqDataSetFromMatrix(
                                  countData = round(.self$count),
                                  colData = .self$group,
                                  design = design))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    suppressMessages(DESeq2::sizeFactors(d) <- sz)
    suppressMessages(d <- DESeq2::estimateDispersions(d))
    if ((is.null(full) && is.null(reduced)) && is.null(contrast)) {
        contrast <- c(colnames(.self$group), unique(as.character(.self$group[, 1])))
    }
    if ((! is.null(contrast)) && (is.null(full) && is.null(reduced))) {
        suppressMessages(d <- DESeq2::nbinomWaldTest(d))
        suppressMessages(res <- DESeq2::results(d, contrast = contrast))
    } else if ((is.null(contrast)) && ((! is.null(full)) && (! is.null(reduced)))) {
        suppressMessages(d <- DESeq2::nbinomLRT(d, full = full, reduced = reduced))
        suppressMessages(res <- DESeq2::results(d))
    } else {
        stop("TCC::ERROR: TCC requires 'full', 'reduced', or 'contrast' arguments when the 'test.method = \"deseq2\". Please set 'full' and 'reduced' arguments or 'contrast' argument. These arguments are the same as that of DESeq2. Please check the DESeq2 vignettes for studying how to set 'full', 'reduced' or 'contrast'.")
    }
    private$stat$p.value <<- res$pvalue
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}




.testByDeseq2.2 = function(design, full, reduced, contrast) {
    suppressMessages(d <- DESeq2::DESeqDataSetFromMatrix(
                                  countData = round(.self$count),
                                  colData = .self$group,
                                  design = design))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    suppressMessages(DESeq2::sizeFactors(d) <- sz)
    suppressMessages(d <- DESeq2::estimateDispersions(d))
    if ((is.null(full) && is.null(reduced)) && is.null(contrast)) {
        full <- as.formula(paste("~", colnames(.self$group))) 
        reduced <- formula(~ 1)
    }
    if ((! is.null(contrast)) && (is.null(full) && is.null(reduced))) {
        suppressMessages(d <- DESeq2::nbinomWaldTest(d))
        suppressMessages(res <- DESeq2::results(d, contrast = contrast))
    } else if ((is.null(contrast)) && ((! is.null(full)) && (! is.null(reduced)))) {
        suppressMessages(d <- DESeq2::nbinomLRT(d, full = full, reduced = reduced))
        suppressMessages(res <- DESeq2::results(d))
    } else {
        stop("TCC::ERROR: TCC requires 'full', 'reduced', or 'contrast' arguments when the 'test.method = \"deseq2\". Please set 'full' and 'reduced' arguments or 'contrast' argument. These arguments are the same as that of DESeq2. Please check the DESeq2 vignettes for studying how to set 'full', 'reduced' or 'contrast'.")
    }
    private$stat$p.value <<- res$pvalue
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}




.testByDeseq2.3 = function(design, full, reduced, contrast) {
    suppressMessages(d <- DESeq2::DESeqDataSetFromMatrix(
                                  countData = round(.self$count),
                                  colData = .self$group,
                                  design = design))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    suppressMessages(DESeq2::sizeFactors(d) <- sz)
    suppressMessages(d <- DESeq2::estimateDispersions(d))
    if ((is.null(full) && is.null(reduced)) && is.null(contrast)) {
        full <- as.formula(paste("~", paste(colnames(.self$group), collapse = "+")))
        reduced <- formula(~ 1)
    }
    if ((! is.null(contrast)) && (is.null(full) && is.null(reduced))) {
        suppressMessages(d <- DESeq2::nbinomWaldTest(d))
        suppressMessages(res <- DESeq2::results(d, contrast = contrast))
    } else if ((is.null(contrast)) && ((! is.null(full)) && (! is.null(reduced)))) {
        suppressMessages(d <- DESeq2::nbinomLRT(d, full = full, reduced = reduced))
        suppressMessages(res <- DESeq2::results(d))
    } else {
        stop("TCC::ERROR: TCC requires 'full', 'reduced', or 'contrast' arguments when the 'test.method = \"deseq2\". Please set 'full' and 'reduced' arguments or 'contrast' argument. These arguments are the same as that of DESeq2. Please check the DESeq2 vignettes for studying how to set 'full', 'reduced' or 'contrast'.")
    }
    private$stat$p.value <<- res$pvalue
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}






.testByDeseq2.4 = function(design, full, reduced, contrast) {
    suppressMessages(d <- DESeq2::DESeqDataSetFromMatrix(
                                  countData = round(.self$count),
                                  colData = .self$group,
                                  design = design))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    suppressMessages(DESeq2::sizeFactors(d) <- sz)
    suppressMessages(d <- DESeq2::estimateDispersions(d))
    if ((is.null(full) && is.null(reduced)) && is.null(contrast)) {
        full <- as.formula(paste("~", paste(colnames(.self$group), collapse = "+")))
        reduced <- as.formula(paste("~", colnames(.self$group)[2]))
    }
    if ((! is.null(contrast)) && (is.null(full) && is.null(reduced))) {
        suppressMessages(d <- DESeq2::nbinomWaldTest(d))
        suppressMessages(res <- DESeq2::results(d, contrast = contrast))
    } else if ((is.null(contrast)) && ((! is.null(full)) && (! is.null(reduced)))) {
        suppressMessages(d <- DESeq2::nbinomLRT(d, full = full, reduced = reduced))
        suppressMessages(res <- DESeq2::results(d))
    } else {
        stop("TCC::ERROR: TCC requires 'full', 'reduced', or 'contrast' arguments when the 'test.method = \"deseq2\". Please set 'full' and 'reduced' arguments or 'contrast' argument. These arguments are the same as that of DESeq2. Please check the DESeq2 vignettes for studying how to set 'full', 'reduced' or 'contrast'.")
    }
    private$stat$p.value <<- res$pvalue
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$q.value <<- p.adjust(private$stat$p.value, method = "BH")
    private$stat$rank <<- rank(private$stat$p.value)
}






if (is.null(design)) {
    design <- as.formula(paste("~", paste(colnames(.self$group), collapse = "+")))
}
test.approach <- .self$.testApproach(paired = paired)
switch(test.approach,
    "1" = .testByDeseq2.1(design, full, reduced, contrast),
    "2" = .testByDeseq2.2(design, full, reduced, contrast),
    "3" = .testByDeseq2.3(design, full, reduced, contrast),
    "4" = .testByDeseq2.4(design, full, reduced, contrast),
    stop("TCC::ERROR: TCC does not support such identification strategy.")
)

})
