TCC$methods(.testByBayseq = function(...) {

.testByBayseq.1p = function(samplesize = NULL, cl = NULL, comparison = NULL) {
    if (is.null(comparison)) {
        comparison <- 1
    } else if(!is.numeric(comparison)) {
        cn <- colnames(.self$group)
        comparison <- (1:length(cn))[comparison == cn]
    }
    ug <- unique(.self$group[, 1])
    cd <- nrow(.self$group) / 2
    el <- colSums(.self$count) * .self$norm.factors
    capture.output(suppressMessages(d <- new("pairedData",
            data = round(.self$count[, 1:cd]),
            pairData = round(.self$count[, (cd + 1):(cd + cd)]),
            replicates = .self$group[1:cd, 2],
            groups = list(NDE = rep(1, length = cd),
                          DE = .self$group[1:cd, 2]),
            libsizes = el[1:cd],
            pairLibsizes = el[(cd + 1):(cd * 2)]
    )))
    capture.output(suppressMessages(d <- getPriors.BB(d, 
             samplesize = samplesize, cl = cl)))
    capture.output(suppressMessages(d <- getLikelihoods.BB(d,
             pET = "BIC", nullProps = 0.5, cl = cl)))
    stat.bayseq <- topCounts(d, group = comparison, number = nrow(.self$count))
    stat.bayseq <- stat.bayseq[rownames(.self$count), ]
    ## private$stat$rank <<- rank(- d@posteriors[, "DE"])
    private$stat$likelihood <<- stat.bayseq$Likelihood
    private$stat$p.value <<- 1 - stat.bayseq$Likelihood
    private$stat$p.value[is.na(private$stat$p.value)] <<- 1
    private$stat$rank <<- rank(private$stat$p.value)
    private$stat$q.value <<- stat.bayseq$FDR
    private$estimatedDEG <<- as.numeric(.self$private$stat$rank < 
                                  (nrow(.self$count) * d@estProps[2]))
    private$tbt$estProps <<- d@estProps[2]
}

.testByBayseq.2 = function(samplesize = NULL, cl = NULL) {
    capture.output(suppressMessages(d <- new("countData",
             data = round(.self$count),
             replicates = .self$group[, 1],
             groups = list(NDE = rep(1, length = nrow(.self$group)),
                           DE = .self$group[, 1]),
             libsizes = colSums(.self$count) * .self$norm.factors)))
    capture.output(suppressMessages(d <- getPriors.NB(d, 
             samplesize = samplesize, estimation = "QL", cl = cl)))
    capture.output(suppressMessages(d <- getLikelihoods.NB(d,
             pET = "BIC", cl = cl)))
    stat.bayseq <- topCounts(d, group = "DE", number = nrow(.self$count))
    stat.bayseq <- stat.bayseq[rownames(.self$count), ]
    private$stat$rank <<- rank(- d@posteriors[, "DE"])
    private$stat$likelihood <<- stat.bayseq$Likelihood
    private$stat$p.value <<- 1 - stat.bayseq$Likelihood
    private$stat$q.value <<- stat.bayseq$FDR
    private$estimatedDEG <<- as.numeric(.self$private$stat$rank < 
                                  (nrow(.self$count) * d@estProps[2]))
    private$tbt$estProps <<- d@estProps[2]
}

.testByBayseq.3 = function(samplesize = NULL, cl = NULL,
                                       comparison = NULL) {
    if (is.null(comparison))
        comparison <- colnames(.self$group)[2]
    gs <- .self$group
    gs <- cbind(rep(1, length = nrow(.self$group)), gs)
    colnames(gs)[1] <- "NDE"
    suppressMessages(d <- new("countData", data = round(.self$count),
             replicates = .self$group[, 1],
             groups = gs,
             libsizes = colSums(.self$count) * .self$norm.factors))
    capture.output(suppressMessages(d <- getPriors.NB(d,
             samplesize = samplesize, estimation = "QL", cl = cl)))
    capture.output(suppressMessages(d <- getLikelihoods.NB(d,
             pET = "BIC", cl = cl)))
    stat.bayseq <- topCounts(d, group = comparison, number = nrow(.self$count))
    stat.bayseq <- stat.bayseq[rownames(.self$count), ]
    private$stat$rank <<- rank(- d@posteriors[, comparison])
    private$stat$likelihood <<- stat.bayseq$Likelihood
    private$stat$p.value <<- 1 - stat.bayseq$Likelihood
    private$stat$q.value <<- stat.bayseq$FDR
    private$estimatedDEG <<- as.numeric(.self$private$stat$rank < 
                                  (nrow(.self$count) * d@estProps[2]))
    private$tbt$estProps <<- d@estProps[2]
}

al <- list(...)
if (is.null(al$samplesize)) samplesize <- 10000
else samplesize <- al$samplesize
if (is.null(al$paired)) al$paired <- FALSE

cl <- al$cl
comparison <- al$comparison
ts <- .self$.testStrategy()
if (al$paired) {
    .testByBayseq.1p(samplesize = samplesize, cl = cl,
                     comparison = comparison)
} else if (ts == 1) {
    .testByBayseq.2(samplesize = samplesize, cl = cl)
} else if (ts == 2) {
    .testByBayseq.2(samplesize = samplesize, cl = cl)
} else if (ts == 3) {
    .testByBayseq.3(samplesize = samplesize, cl = cl,
                    comparison = comparison)
} else {
   stop()
}

})


