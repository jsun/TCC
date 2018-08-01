if (FALSE) {
TCC$methods(.testByYayoi = function(upper.limit = 0.4, deseq = NULL, norm = FALSE, ...) {
    if (is.null(deseq)) {
        deseq <- 1
    }
    if (norm) {
        ## if normalization, treat all libraries as independent will easily remove the outliers from data.
        yayoi_group = 1:length(.self$group[, 1])
    } else {
        yayoi_group = .self$group[, 1]
    }

    r <- YAYOI(count = .self$count,
               group = yayoi_group,
               norm.factors = .self$norm.factors,
               upper.limit = upper.limit,
               deseq = deseq, 
               ...)
    private$stat$rank <<- rank(- r$score)
    private$stat$testStat <<- r$score
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- rep(NA, length = nrow(.self$count))
    private$estimatedDEG <<- rep(0, length = nrow(.self$count))
    private$outmx <<- r$outlier
    private$prob <<- r$prob
})
}
