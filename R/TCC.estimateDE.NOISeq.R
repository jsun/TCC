TCC$methods(.testByNoiseq = function(...) {

.testByNoiseq.1 = function() {
    x <- .self$getNormalizedData()
    gl <- data.frame(group = .self$group[, 1])
    nd <- NOISeq::readData(x, gl)
    capture.output(suppressMessages(nr <- NOISeq::noiseq(nd,
                           k = 0.5,
                           norm = "n",
                           replicates = "biological",
                           factor = "group", 
                           conditions = unique(.self$group[, 1]))))
    prob <- nr@results[[1]]$prob
    prob[is.na(prob)] <- 0
    private$stat$prob <<- prob
    private$stat$p.values <<- rep(NA, length = nrow(.self$count))
    private$stat$q.values <<- rep(NA, length = nrow(.self$count))
    private$stat$rank <<- rank(- prob)
}

ts <- .self$.testStrategy()
if (ts == 1) {
   .testByNoiseq.1()
} else if (ts == 2) {
   stop()
} else if (ts == 3) {
   stop()
} else {
   stop()
}

})
