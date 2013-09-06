TCC$methods(.testByEbseq = function(...) {

.testByEbseq.1 = function(samplesize = NULL) {
    g <- .self$group[, 1]
    ug <- unique(g)
    suppressMessages(EBout <- EBSeq::EBTest(Data = .self$count,
                     Conditions = as.factor(g),
                     sizeFactors = .self$norm.factors * colSums(.self$count),
                     maxround = samplesize))
    PP <- EBSeq::GetPPMat(EBout)
    df <- matrix(1, ncol = 2, nrow = nrow(.self$count))
    rownames(df) <- rownames(.self$count)
    df[rownames(PP), 1] <- PP[, 1]
    df[rownames(PP), 2] <- PP[, 2]
    df[is.na(df)] <- 0
    private$stat$prob <<- df[, 2]
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- df[, 1]
    private$stat$rank <<- rank(- .self$private$stat$prob)
}

.testByEbseq.2 = function(samplesize = NULL) {
    g <- .self$group[, 1]
    ug <- unique(g)
    gp <- matrix(c(rep(1, length = length(ug)), 1:length(ug)),
                 nrow = 2, byrow = TRUE)
    colnames(gp) <- ug
    rownames(gp) <- c("Pattern1", "Pattern2")
    suppressMessages(MultiOut <- EBSeq::EBMultiTest(.self$count,
                     NgVector = NULL,
                     Conditions = g,
                     AllParti = gp,
                     sizeFactors = .self$norm.factors * colSums(.self$count),
                     maxround = samplesize))
    PP <- EBSeq::GetMultiPP(MultiOut)
    df <- matrix(1, ncol = 2, nrow = nrow(.self$count))
    rownames(df) <- rownames(.self$count)
    df[rownames(PP$PP), 1] <- PP$PP[, 1]
    df[rownames(PP$PP), 2] <- PP$PP[, 2]
    df[is.na(df)] <- 0
    private$stat$prob <<- df[, 2]
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- df[, 1]
    private$stat$rank <<- rank(- .self$private$stat$prob)
}


al <- list(...)
if (is.null(al$samplesize)) {
    samplesize <- 5
} else {
    samplesize <- al$samplesize
}
ts <- .self$.testStrategy()
if (ts == 1) {
   .testByEbseq.1(samplesize = samplesize)
} else if (ts == 2) {
   .testByEbseq.2(samplesize = samplesize)
} else if (ts == 3) {
   stop()
} else {
   stop()
}
})

