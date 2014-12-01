TCC$methods(.testByBayseq = function(samplesize = NULL, cl = NULL, paired = NULL, bgroup = NULL, ...) {





.testByBayseq.1 = function(samplesize, cl, ...) {
    capture.output(d <- new("countData", data = round(.self$count),
                            replicates = .self$group[, 1],
                            groups = list(NDE = rep(1, nrow(.self$group)),
                                          DE = .self$group[, 1]),
                            libsizes = .self$norm.factors * colSums(.self$count)))
    suppressMessages(capture.output(d <- baySeq::getPriors.NB(d,
                                    samplesize = samplesize, cl = cl)))
    suppressMessages(capture.output(d <- baySeq::getLikelihoods(d, cl = cl)))
    res <- topCounts(d, group = "DE", number = nrow(.self$count))
    res <- res[order(res$rowID), ]
    private$stat$p.value <<- res$FWER.DE
    private$stat$q.value <<- res$FDR.DE
    private$stat$rank <<- rank(res$FWER.DE)
}




.testByBayseq.3 = function(samplesize, cl, ...) {
    args <- list(...)
    if (is.null(args$bgroup)) {
        g <- colnames(.self$group)[1]
    } else {
        g <- args$bgroup
    }
    grps <- cbind(rep(1, length = nrow(.self$group)), .self$group)
    colnames(grps) <- c("NDE", colnames(.self$group))
    suppressMessages(d <- new("countData", data = round(.self$count),
             replicates = .self$group[, 1],
             groups = grps,
             libsizes = colSums(.self$count) * .self$norm.factors))
    suppressMessages(capture.output(d <- baySeq::getPriors.NB(d, 
                                    samplesize = samplesize, cl = cl)))
    suppressMessages(capture.output(d <- baySeq::getLikelihoods(d, cl = cl)))
    res <- topCounts(d, group = g, number = nrow(.self$count))
    res <- res[order(res$rowID), ]
    private$stat$p.value <<- res[, paste0("FWER.", g)]
    private$stat$q.value <<- res[, paste0("FDR.", g)]
    private$stat$rank <<- rank(res[, paste0("FWER.", g)])
}



.testByBayseq.4 = function(samplesize, cl, ...) {
    n.paires <- nrow(.self$group) / 2
    count_1 <- count[, 1:n.paires]
    count_2 <- count[, (n.paires + 1):(ncol(count))]
    libsize_1 <- colSums(count_1) * .self$norm.factors[1:n.paires]
    libsize_2 <- colSums(count_2) * .self$norm.factors[(n.paires + 1):(ncol(count))]
    grps <- .self$group[1:n.paires, 2]
    patterns <- list(NDE = rep(1, n.paires), DE = grps)
    d <- new("countData", data = list(count_1, count_2),
              replicates = grps, groups = patterns,
              densityFunction = bbDensity,
              libsizes = cbind(libsize_1, libsize_2))
    suppressMessages(capture.output(d <- baySeq::getPriors(d,
                                    samplesize = samplesize, cl = cl)))
    suppressMessages(capture.output(d <- baySeq::getLikelihoods(d, cl = cl)))
    res <- topCounts(d, group = "NDE", number = nrow(.self$count))
    res <- res[order(res$rowID), ]
    private$stat$p.value <<- res$FWER.NDE
    private$stat$q.value <<- res$FDR.NDE
    private$stat$rank <<- rank(res$FWER.NDE)
}





##
## main process
##
if (is.null(samplesize))  samplesize <- 10000

test.approach <- .self$.testApproach(paired = paired)

switch(test.approach, 
    "1" = .testByBayseq.1(samplesize = samplesize, cl = cl, ...),
    "2" = .testByBayseq.1(samplesize = samplesize, cl = cl, ...),
    "3" = .testByBayseq.3(samplesize = samplesize, cl = cl,
                          bgroup = bgroup, ...),
    "4" = .testByBayseq.4(samplesize = samplesize, cl = cl, ...),
    stop("TCC::ERROR: TCC does not support such identification strategy.")
)


})


