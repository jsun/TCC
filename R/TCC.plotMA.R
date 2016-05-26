
TCC$methods(plotMA = function (FDR = NULL,
                       significance.level = NULL,
                       median.lines = FALSE,
                       floor = 0,
                       group = NULL,
                       col = NULL,
                       col.tag = NULL,
                       normalize = TRUE,
                       showfig = TRUE,
                       ...) {
    fcall <- as.list(match.call(expand.dots = TRUE))
    arglist <- list(...)
    if (is.null(arglist$cex))
        arglist$cex <- 0.3
    if (is.null(arglist$pch))
        arglist$pch <- 20
    if (is.null(arglist$main))
        arglist$main <- "MA plot"

    ## set up default arguments.
    ggroups <- eval(parse(text = fcall$group))
    gro <- .self$group[, 1]
    gru <- unique(as.vector(gro))
    if (is.null(ggroups)) {
        ggroups <- c(gru[1], gru[2])
    }
    if (is.null(arglist$xlab)) {
        g1name <- ggroups[1]
        g2name <- ggroups[2]
        #suppressWarnings(ggname.is.na <- !any(is.na(as.integer(as.character(ggroups)))))
        #if (.self$private$simulation && ggname.is.na) {
        #    g1name <- paste0("G", g1name)
        #    g2name <- paste0("G", g2name)
        #}
        g1name <- paste0("G", g1name)
        g2name <- paste0("G", g2name)
        #arglist$xlab <- paste0("A = (log2(", g2name, ") + log2(", g1name, "))/2")
        arglist$xlab <- bquote(A ~ "=" ~ (log[2](.(g2name)) + log[2](.(g1name))) / 2)
    }
    if (is.null(arglist$ylab)) {
        g1name <- ggroups[1]
        g2name <- ggroups[2]
        
        ## if the group name is give, use group name,
        ## otherwise, use G1, G2, ... as group name.
        #suppressWarnings(ggname.is.na <- !any(is.na(as.integer(as.character(ggroups)))))
        #if (.self$private$simulation && ggname.is.na) {
        #    g1name <- paste0("G", g1name)
        #    g2name <- paste0("G", g2name)
        #}
        
        ## always use G1, G2, ... as group name.
        g1name <- paste0("G", g1name)
        g2name <- paste0("G", g2name)
        
        #arglist$ylab <- paste0("M = log2(", g2name, ") - log2(", g1name, ")")
        arglist$ylab <- bquote(M ~ "=" ~ log[2](.(g2name)) - log[2](.(g1name)))
    }

    if (is.null(col)) {
        if (private$estimated == TRUE) {
            arglist$col <- c("black", rep("magenta", length = length(gru)))
        } else if (private$simulation == TRUE) {
            arglist$col <- 1:length(gru)
        } else {
            arglist$col <- rep("black", length = length(gru))
        }
    } else {
        arglist$col <- col
    }
    if (normalize)
      count.normed <- .self$getNormalizedData()  
    else
      count.normed <- .self$count
    mean.i <- rowMeans(as.matrix(count.normed[, gro == ggroups[1]]))
    mean.j <- rowMeans(as.matrix(count.normed[, gro == ggroups[2]]))
    norm.i <- mean(norm.factors[gro == ggroups[1]])
    norm.j <- mean(norm.factors[gro == ggroups[2]])
    ma.axes <- .self$.getMACoordinates(mean.i, mean.j, floor)
    filter <- as.logical(mean.i > 0 & mean.j > 0)
    a <- ma.axes$a.value
    m <- ma.axes$m.value
    ## plot the M-A plot (default)
    if (showfig) {
        if (is.null(arglist$xlim))
            arglist$xlim <- c(min(a), max(a))
        if (is.null(arglist$ylim))
            arglist$ylim <- c(min(m), max(m))
        arglist$x <- c(0, 0)
        arglist$type <- "n"
        do.call(plot, arglist)
        grid(col = "gray", lty = "dotted")

        if (is.null(col.tag)) {
            col.tag.v <- rep(1, length = nrow(count))
            if (private$estimated == FALSE) {
                if (private$simulation == TRUE)
                    col.tag.v <- simulation$trueDEG + 1
            } else {
                if ((!is.null(.self$estimatedDEG)) && (length(.self$estimatedDEG != 0))) {
                    col.tag.v <- as.numeric(estimatedDEG) + 1
                }
                if (!(is.null(FDR) && is.null(significance.level))) {
                    private$stat$q.value <<- stat$q.value
                    private$stat$p.value <<- stat$p.value
                    col.tag.v <- .self$.exactTest(FDR = FDR, 
                                      significance.level = significance.level) + 1
                }
            }
            col.tag <- col.tag.v
        }
        if (length(col.tag) != nrow(count))
            stop("\nTCC::ERROR: The length of col.tag has to be equal to the number of genes.\n")
        for (k in unique(col.tag)) {
            points(a[col.tag == k], m[col.tag == k], 
               col = arglist$col[k], pch = arglist$pch, cex = arglist$cex)
        }
        if (median.lines == TRUE) {
            if (.self$private$simulation || length(groups) == 2) {
                upG2 <- as.logical(m > 0)
                upG1 <- as.logical(m < 0)
                # nonDEGs
                med <- median(m[(col.tag == 1 & filter)])
                lines(c(min(a) + 1, max(a)), c(med, med), col = arglist$col[1])
                text(arglist$xlim[2], med + 0.5, sprintf("%.3f", med), col = arglist$col[1],
                     pos = 2, offset = 0)
                # DEGs in G1
                med <- median(m[((col.tag == 2 & filter) & upG1)])
                lines(c(min(a) + 1, max(a)), c(med, med), col = arglist$col[2])
                text(arglist$xlim[2], med + 0.5, sprintf("%.3f", med), col = arglist$col[2],
                     pos = 2, offset = 0)
                # DEGs in G2
                med <- median(m[((col.tag == 2 & filter) & upG2)])
                lines(c(min(a) + 1, max(a)), c(med, med), col = arglist$col[2])
                text(arglist$xlim[2], med + 0.5, sprintf("%.3f", med), col = arglist$col[2],
                     pos = 2, offset = 0)
            } else {
                warnings("TCC::WARN:: TCC ploted MA-plot without the 'median.lines = TRUE' option, this option only for two-group simulation data created by 'simulateReadCounts' function.\n")
            }
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

