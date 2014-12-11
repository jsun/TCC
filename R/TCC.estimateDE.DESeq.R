TCC$methods(.testByDeseq = function(full = NULL, reduced = NULL,
                                    paired = NULL, ... ) {





.testByDeseq.1 = function(full, reduced, ...) {
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)

    if (is.null(full) && is.null(reduced)) {
        suppressMessages(d <- newCountDataSet(countData = round(.self$count), 
                                              conditions = .self$group[, 1]))
        sizeFactors(d) <- sz
        if ((length(unique(.self$group[, 1])) == nrow(.self$group)) &&
            (length(list(...)) == 0)) {
            suppressMessages(d <- estimateDispersions(d, method = "blind",
                                                  sharingMode = "fit-only"))
        } else {
            suppressMessages(d <- estimateDispersions(d, ...))
        }
        suppressMessages(d <- nbinomTest(d, levels(.self$group[, 1])[1],
                                            levels(.self$group[, 1])[2]))
        p <- d$pval
    } else {
        suppressMessages(d <- newCountDataSet(countData = round(.self$count), 
                                              conditions = .self$group[, 1]))
        sizeFactors(d) <- sz
        if (length(unique(.self$group[, 1])) == nrow(.self$group)) {
            suppressMessages(d <- estimateDispersions(d, method = "blind",
                                                  sharingMode = "fit-only"))
        } else {
            suppressMessages(d <- estimateDispersions(d, ...))
        }
        if (is.na(match("count", as.character(full))))
            full <- as.formula(paste(c("count", as.character(full)), collapse=""))
        if (is.na(match("count", as.character(reduced))))
            reduced <- as.formula(paste(c("count", as.character(reduced)), collapse=""))
        capture.output(f.model <- fitNbinomGLMs(d, full))
        capture.output(r.model <- fitNbinomGLMs(d, reduced))
        p <- nbinomGLMTest(f.model, r.model)
    }
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$q.value <<- p.adjust(p, method = "BH")
    private$stat$rank <<- rank(p)
}





.testByDeseq.2 = function(full, reduced, ...) {
    if (is.null(full)) full <- formula(~ condition) 
    if (is.null(reduced)) reduced <- formula(~ 1)
    suppressMessages(d <- newCountDataSet(countData = round(.self$count), 
                                          conditions = .self$group[, 1]))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    sizeFactors(d) <- sz
    if ((length(unique(.self$group[, 1])) == nrow(.self$group)) &&
        (length(list(...)) == 0)) {
        suppressMessages(d <- estimateDispersions(d, method = "blind",
                                                  sharingMode = "fit-only"))
    } else {
        suppressMessages(d <- estimateDispersions(d, ...))
    }
    if (is.na(match("count", as.character(full))))
        full <- as.formula(paste(c("count", as.character(full)), collapse = ""))
    if (is.na(match("count", as.character(reduced))))
        reduced <- as.formula(paste(c("count", as.character(reduced)), collapse = ""))
    capture.output(f.model <- fitNbinomGLMs(d, full))
    capture.output(r.model <- fitNbinomGLMs(d, reduced))
    p <- nbinomGLMTest(f.model, r.model)
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$q.value <<- p.adjust(p, method = "BH")
    private$stat$rank <<- rank(p)
}





.testByDeseq.3 = function(full, reduced, ...) {
    if (is.null(full)) full <- as.formula(paste("~",
                               paste(colnames(.self$group), collapse = "+"))) 
    if (is.null(reduced)) reduced <- formula(~ 1)
    suppressMessages(d <- newCountDataSet(countData = round(.self$count), 
                                          conditions = .self$group))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    sizeFactors(d) <- sz
    gtags <- apply(.self$group, 1, paste, collapse = "")

    if ((length(unique(gtags)) == length(gtags)) &&
        (length(list(...)) == 0)) {
        suppressMessages(d <- estimateDispersions(d, method = "blind",
                                                  sharingMode = "fit-only"))
    } else {
        suppressMessages(d <- estimateDispersions(d, ...))
    }
    if (is.na(match("count", as.character(full))))
        full <- as.formula(paste(c("count", as.character(full)), collapse=""))
    if (is.na(match("count", as.character(reduced))))
        reduced <- as.formula(paste(c("count", as.character(reduced)), collapse=""))
    capture.output(f.model <- fitNbinomGLMs(d, full))
    capture.output(r.model <- fitNbinomGLMs(d, reduced))
    p <- nbinomGLMTest(f.model, r.model)
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$q.value <<- p.adjust(p, method = "BH")
    private$stat$rank <<- rank(p)
}







.testByDeseq.4 <- function(full, reduced, ...) {
    if (is.null(full)) full <- as.formula(paste("~",
                               paste(colnames(.self$group), collapse = "+"))) 
    if (is.null(reduced)) reduced <- as.formula(paste("~", colnames(.self$group)[2]))
    suppressMessages(d <- newCountDataSet(countData = round(.self$count), 
                                          conditions = .self$group))
    ef.libsizes <- .self$norm.factors * colSums(.self$count)
    sz <- ef.libsizes / mean(ef.libsizes)
    sizeFactors(d) <- sz
    gtags <- apply(.self$group, 1, paste, collapse = "")
    if ((length(unique(gtags)) == length(gtags)) &&
        (length(list(...)) == 0)) {
        suppressMessages(d <- estimateDispersions(d, method = "blind",
                                                  sharingMode = "fit-only"))
    } else {
        suppressMessages(d <- estimateDispersions(d, ...))
    }
    if (is.na(match("count", as.character(full))))
        full <- as.formula(paste(c("count", as.character(full)), collapse = ""))
    if (is.na(match("count", as.character(reduced))))
        reduced <- as.formula(paste(c("count", as.character(reduced)), collapse = ""))
    capture.output(f.model <- fitNbinomGLMs(d, full))
    capture.output(r.model <- fitNbinomGLMs(d, reduced))
    p <- nbinomGLMTest(f.model, r.model)
    p[is.na(p)] <- 1
    private$stat$p.value <<- p
    private$stat$q.value <<- p.adjust(p, method = "BH")
    private$stat$rank <<- rank(p)
}





test.approach <- .self$.testApproach(paired = paired)

switch(test.approach,
    "1" = .testByDeseq.1(full = full, reduced = reduced, ...),
    "2" = .testByDeseq.2(full = full, reduced = reduced, ...),
    "3" = .testByDeseq.3(full = full, reduced = reduced, ...),
    "4" = .testByDeseq.4(full = full, reduced = reduced, ...),
    stop("TCC::ERROR: TCC does not support such identification strategy.")
)


})

