TCC$methods(.testBySamseq = function(...) {
require(samr)

if(! requireNamespace("samr", quietly=TRUE)) {
    stop("TCC::ERROR: TCC needs the 'samr' package for DEG identification if 'test.method' is specifed to 'samseq'. Please install the 'samr' package from CRAN.")
}



.testBySamseq.1 <- function(samplesize = NULL) {
    suppressMessages(capture.output(res <- samr::SAMseq(x = .self$count,
                y = .self$group[, 1],
                resp.type = "Two class unpaired",
                nperms = samplesize)))
    pval <- samr::samr.pvalues.from.perms(res$samr.obj$tt, res$samr.obj$ttstar)
    private$stat$p.value <<- pval
    private$stat$q.value <<- p.adjust(pval, method = "BH")
    private$stat$rank <<- rank(pval)
}






.testBySamseq.2 <- function(samplesize = NULL) {
    suppressMessages(capture.output(res <- samr::SAMseq(x = .self$count,
                y = .self$group[, 1],
                resp.type = "Multiclass",
                nperms = samplesize)))
    pval <- samr::samr.pvalues.from.perms(res$samr.obj$tt, res$samr.obj$ttstar)
    private$stat$p.value <<- pval
    private$stat$q.value <<- p.adjust(pval, method = "BH")
    private$stat$rank <<- rank(pval)
}





.testBySamseq.4 <- function(samplesize = NULL) {
    nreps <- nrow(.self$group) / 2
    g.vec <- c(-(1:nreps), 1:nreps)
    suppressMessages(capture.output(res <- samr::SAMseq(x = .self$count,
                y = g.vec,
                resp.type = "Two class paired",
                nperms = samplesize)))
    pval <- samr::samr.pvalues.from.perms(res$samr.obj$tt, res$samr.obj$ttstar)
    private$stat$p.value <<- pval
    private$stat$q.value <<- p.adjust(pval, method = "BH")
    private$stat$rank <<- rank(pval)
}





##
## main process
##
message("TCC::INFO: TCC does not use the normalization factors when 'test.method = \"samseq\"'.")

add.args <- list(...)
if (is.null(add.args$samplesize)) {
    samplesize <- 100
} else {
    samplesize <- add.args$samplesize
}


test.approach <- .self$.testApproach()

switch(test.approach,
    "1" = .testBySamseq.1(samplesize = samplesize),
    "2" = .testBySamseq.2(samplesize = samplesize),
    "4" = .testBySamseq.4(samplesize = samplesize),
    stop("TCC::ERROR: TCC does not support such identification strategy.")

)

})


