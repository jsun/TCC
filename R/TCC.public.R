setMethod(
    f = "calcNormFactors",
    signature(tcc = "TCC"),
    definition = function(tcc, norm.method = NULL, test.method = NULL, 
                          iteration = TRUE, FDR = NULL, floorPDEG = NULL,
                          increment = FALSE, ...) {
        obj <- tcc$copy()
        obj$calcNormFactors(norm.method = norm.method, 
                            test.method = test.method, 
                            iteration = iteration, FDR = FDR, 
                            floorPDEG = floorPDEG, increment = increment,
                            ...)
        return(obj)
    }
)

estimateDE <- function(tcc, test.method = NULL, FDR = NULL, paired = NULL,
                       full = NULL, reduced = NULL,
                       design = NULL,  contrast = NULL,
                       coef = NULL, 
                       group = NULL, cl = NULL,
                       samplesize = NULL,
                       logged = NULL, floor = NULL, ...) {
    obj <- tcc$copy()
    obj$estimateDE(test.method = test.method, FDR = FDR, paired = paired,
                   full = full, reduced = reduced, 
                   design = design, contrast = contrast, coef = coef,
                   group = group, samplesize = samplesize,
                   logged = logged, floor = floor, cl = cl, ...)
    return(obj)
}

getResult <- function(tcc, sort = FALSE, ...) {
    if (length(tcc$stat) == 0)
        stop("\nTCC::ERROR: There are no statistics in stat fields of TCC class tcc. Execute TCC.estiamteDE for calculating them.\n")
    ## calculate M-A coordinates
    gru <- unique(tcc$group[, 1])
    m.value <- rep(NA, length = nrow(tcc$count))
    a.value <- rep(NA, length = nrow(tcc$count))
    if ((length(gru) == 2) && (ncol(tcc$group) == 1)) {
      #ma <- tcc$plotMA(showfig = FALSE, ...)
      ma <- plot(tcc, showfig = FALSE, ...)
      m.value <- ma$m.value
      a.value <- ma$a.value
    }
    ## merge statistics to data frame
    if (!is.null(tcc$stat$p.value)) {
        ## show p-values if existed
        df <- data.frame(
            gene_id = rownames(tcc$count),
            a.value = a.value,
            m.value = m.value,
            p.value = tcc$stat$p.value,
            q.value = tcc$stat$q.value,
            rank    = tcc$stat$rank,
            estimatedDEG = tcc$estimatedDEG
        ) 
    } else if (!is.null(tcc$stat$testStat)) {
        ## show probability if existed
        df <- data.frame(
            gene_id  = rownames(tcc$count),
            a.value  = a.value,
            m.value  = m.value,
            testStat = tcc$stat$testStat,
            rank     = tcc$stat$rank,
            estimatedDEG = tcc$estimatedDEG
        ) 
    } else if (!is.null(tcc$stat$prob)) {
        ## show probability if existed
        df <- data.frame(
            gene_id = rownames(tcc$count),
            a.value = a.value,
            m.value = m.value,
            prob    = tcc$stat$prob,
            rank    = tcc$stat$rank,
            estimatedDEG = tcc$estimatedDEG
        ) 
    }
    if (sort)
        df <- df[order(df$rank), ]
    rownames(df) <- NULL
    return (df)
}

filterLowCountGenes <- function(tcc, low.count = 0) {
    obj <- tcc$copy()
    gru <- unique(obj$group[, 1])
    filters <- matrix(0, ncol = length(gru), nrow = nrow(obj$count)) 
    for (i in 1:length(gru)) {
        filters[, i] <- as.numeric(rowSums(
                            as.matrix(obj$count[, (obj$group[, 1] == gru[i])])
                        ) <= low.count)
    }
    left.tag <- as.logical(rowSums(filters) != length(gru))
    obj$count <- obj$count[left.tag, ]
    if (!is.null(obj$simulation$trueDEG) && length(obj$simulation$trueDEG) != 0)
        obj$simulation$trueDEG <- obj$simulation$trueDEG[left.tag]
    if (!is.null(obj$estimatedDEG) && length(obj$estimatedDEG) != 0)
        obj$estimatedDEG <- obj$estimatedDEG[left.tag]
    if (!is.null(obj$stat) && length(obj$stat) != 0) {
        for (i in 1:length(obj$stat)) {
            if (length(obj$stat[[i]]) == length(left.tag))
                obj$stat[[i]] <- obj$stat[[i]][left.tag]
        }
    }
    return (obj)
}

getNormalizedData <- function(tcc) {
    return (tcc$getNormalizedData())
}

