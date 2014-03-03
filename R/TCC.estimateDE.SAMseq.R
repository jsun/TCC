TCC$methods(.testBySamseq = function(...) {





.testBySamseq.1 = function(samplesize = NULL) {
    c <- round(.self$getNormalizedData())
    s <- suppressMessages(samr::SAMseq(x = c, y = .self$group[, 1],
                resp.type = "Two class unpaired",
                nperms = samplesize))
    private$stat$testStat <<- s$samr.obj$tt
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- rep(NA, length = nrow(.self$count))
    private$stat$rank <<- rank(- abs(s$samr.obj$tt))
}





.testBySamseq.1p = function(samplesize = NULL) {
    c <- round(.self$getNormalizedData())
    s <- suppressMessages(samr::SAMseq(x = c, y = .self$group[, 1],
                resp.type = "Two class paired",
                nperms = samplesize))
    private$stat$testStat <<- s$samr.obj$tt
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- rep(NA, length = nrow(.self$count))
    private$stat$rank <<- rank(- abs(s$samr.obj$tt))
}





.testBySamseq.2 = function(samplesize = NULL) {
    c <- round(.self$getNormalizedData())
    s <- suppressMessages(samr::SAMseq(x = c, y = .self$group[, 1],
                resp.type = "Multiclass",
                nperms = samplesize))
    private$stat$testStat <<- s$samr.obj$tt
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- rep(NA, length = nrow(.self$count))
    private$stat$rank <<- rank(- abs(s$samr.obj$tt))
}





##
## main process
##
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
    stop("TCC::ERROR: TCC kernel error.")

)

})


