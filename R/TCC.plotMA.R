
TCC$methods(plotMA = function (FDR = NULL,
                       significance.level = NULL,
                       median.lines = FALSE,
                       floor = 0,
                       groups = NULL,
                       col.tag = NULL,
                       normalize = TRUE, ...) {
    arglist <- list(...)
    if (is.null(arglist$xlab))
        arglist$xlab <- expression(A == (log[2] * G2 + log[2] * G1 ) / 2)
    if (is.null(arglist$ylab))
        arglist$ylab <- expression(M == log[2] * G2 - log[2] * G1)
    if (is.null(arglist$cex))
        arglist$cex <- 0.3
    if (is.null(arglist$pch))
        arglist$pch <- 20
    if (is.null(arglist$main))
        arglist$main <- "MA plot"

    ## set up default arguments.
    gro <- .self$group[, 1]
    gru <- unique(as.vector(gro))
    if (is.null(groups)) {
        groups <- c(gru[1], gru[2])
    }
    if (is.null(arglist$col)) {
        if (private$estimated == TRUE) {
            arglist$col <- c(1, rep(6, length = length(gru)))
        } else if (private$simulation == TRUE) {
            arglist$col <- c(1, 4, 2, 4 + 1:(length(gru)))
        } else {
            arglist$col <- rep(1, length = length(gru))
        }
    }
    if (normalize)
      count.normed <- .self$getNormalizedData()  
    else
      count.normed <- .self$count
    mean.i <- rowMeans(as.matrix(count.normed[, gro == groups[1]]))
    mean.j <- rowMeans(as.matrix(count.normed[, gro == groups[2]]))
    norm.i <- mean(norm.factors[gro == groups[1]])
    norm.j <- mean(norm.factors[gro == groups[2]])
    ma.axes <- .self$.getMACoordinates(mean.i, mean.j, floor)
    filter <- as.logical(mean.i > 0 & mean.j > 0)
    a <- ma.axes$a.value
    m <- ma.axes$m.value
  
    if (is.null(arglist$xlim))
        arglist$xlim <- c(min(a), max(a))
    if (is.null(arglist$ylim))
        arglist$ylim <- c(min(m), max(m))
    arglist$x <- c(0, 0)
    arglist$type <- "n"
    do.call(plot, arglist)
    grid(col = "gray", lty = "dotted")
    col.tag.v <- rep(0, length = nrow(count))
    if (private$estimated == FALSE) {
        if (private$simulation == TRUE)
            col.tag.v <- simulation$trueDEG
    } else {
        if ((!is.null(estimatedDEG)) && (length(estimatedDEG != 0))) {
            col.tag.v <- as.numeric(estimatedDEG)
        }
        if (!(is.null(FDR) && is.null(significance.level))) {
            private$stat$q.value <<- stat$q.value
            private$stat$p.value <<- stat$p.value
            col.tag.v <- .self$.exactTest(FDR = FDR, 
                                  significance.level = significance.level)
        }
    }
    if (is.null(col.tag))
        col.tag <- col.tag.v + 1
    if (length(col.tag) != nrow(count))
        stop("\nTCC::ERROR: The length of col.tag has to be equal to the number of genes.\n")
    for (k in unique(col.tag)) {
        points(a[col.tag == k], m[col.tag == k], 
               col = arglist$col[k], pch = arglist$pch, cex = arglist$cex)
    }
    if (median.lines == TRUE) {
        for (k in unique(col.tag)) {
            if (length(setdiff(gru, groups)) != 0 && k == setdiff(gru, groups))
              next
            med <- median(m[(col.tag == k & filter)])
            lines(c(min(a) + 1, max(a)), c(med, med), col = arglist$col[k])
            text(arglist$xlim[2], med + 0.5, sprintf("%.3f", med), col = arglist$col[k],
                 pos = 2, offset = 0)
        }
    }
    invisible(data.frame(a.value = a, m.value = m))
})

TCC$methods(.getMACoordinates = function(g1, g2, floor = 0) {
    m <- rep(0, length = nrow(count))
    a <- rep(0, length = nrow(count))
    g1.min.nonZero <- min(g1[g1 > 0])
    g2.min.nonZero <- min(g2[g2 > 0])
    filter <- as.logical(g1 <= floor | g2 <= floor)
    g1[g1 <= floor] <- g1.min.nonZero
    g2[g2 <= floor] <- g2.min.nonZero
    a <- (log2(g1) + log2(g2)) / 2
    m <- log2(g2) - log2(g1)
    a[filter] <- min(a) - 1
    return(list(m.value = m, a.value = a))
})

