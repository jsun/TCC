setMethod(
    f = "calcNormFactors",
    signature(tcc = "TCC"),
    definition = function(tcc, norm.method = NULL, test.method = NULL, 
                          iteration = TRUE, FDR = NULL, floorPDEG = 0.05,
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

estimateDE <- function(tcc, test.method = NULL, FDR = NULL,
                       dispersion = NULL,
                       fit0 = NULL, fit1 = NULL, design = NULL, contrast=NULL,
                       coef = NULL, comparison = NULL, samplesize = NULL,
                       floor.value = 1, cl = NULL) {
    obj <- tcc$copy()
    obj$estimateDE(test.method=test.method, FDR=FDR,
                   dispersion=dispersion,
                   fit0=fit0, fit1=fit1,
                   design=design, contrast=contrast, coef=coef,
                   comparison=comparison, samplesize=samplesize,
                   floor.value = floor.value, cl=cl)
    return(obj)
}

getResult <- function(tcc, sort = FALSE, floor = 0) {
    if (length(tcc$stat) == 0)
        stop("\nTCC::ERROR: There are no statistics in stat fields of TCC class tcc. Execute TCC.estiamteDE for calculating them.\n")
    ## calculate M-A coordinates
    gru <- unique(tcc$group[, 1])
    m.value <- rep(NA, length = nrow(tcc$count))
    a.value <- rep(NA, length = nrow(tcc$count))
    if ((length(gru) == 2) && (ncol(tcc$group) == 1)) {
        count.normed <- tcc$getNormalizedData()
        mean.exp <- matrix(0, ncol = length(gru), nrow = nrow(tcc$count))
        gru <- unique(as.vector(tcc$group[, 1]))
        mean.i <- rowMeans(as.matrix(count.normed[, tcc$group[, 1] == gru[1]]))
        mean.j <- rowMeans(as.matrix(count.normed[, tcc$group[, 1] == gru[2]]))
        ma.axes <- tcc$.getMACoordinates(mean.i, mean.j, floor)
        m.value <- ma.axes$m.value
        a.value <- ma.axes$a.value
    }
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
    rownames(df) <- NULL
    if (sort)
        df <- df[order(df$rank), ]
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

