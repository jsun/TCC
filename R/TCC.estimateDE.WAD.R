TCC$methods(.testByWad = function(logged = NULL, floor = NULL, ...) {
    if (is.null(logged)) {
        logged <- FALSE
    }
    if (is.null(floor)) {
        floor <- 1
    }
    x <- .self$count
    # If normalization factors were set, normalize data with them.
    if (sum(.self$norm.factors == 1) != length(.self$norm.factors)) {
        ef <- colSums(x) * .self$norm.factors
        x <- sweep(x, 2, mean(ef) / ef, "*")
    }
    s <- .wad(x = x,
              group = .self$group[, 1],
              logged = logged,
              floor = floor)
    private$stat$rank <<- rank(- abs(s))
    private$stat$testStat <<- s
    private$stat$p.value <<- rep(NA, length = nrow(.self$count))
    private$stat$q.value <<- rep(NA, length = nrow(.self$count))
    private$estimatedDEG <<- rep(0, length = nrow(count))
})


