##TCC$methods(.testByNbpseq = function() {
##    ts <- .self$.testStrategy()
##    if (ts == 1) {
##       .self$.testByNbpseq.1()
##    } else if (ts == 2) {
##       stop()
##    } else if (ts == 3) {
##       stop()
##    } else {
##       stop()
##    }
##})



##TCC$methods(.testByNbpseq.1 = function() {
##    g <- .self$group[, 1]
##    ug <- unique(g)
##    nbp <- NBPSeq::nbp.test(counts = .self$count,
##                    grp.ids = g,
##                    grp1 = ug[1], grp2 = ug[2],
##                    norm.factors = .self$norm.factors,
##                    print.level = 0)
##    private$stat$p.values <<- nbp$p.values
##    private$stat$q.values <<- nbp$q.values
##    private$stat$rank <<- rank(.self$private$stat$p.value)
##})



