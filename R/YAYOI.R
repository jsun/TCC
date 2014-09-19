YAYOI <- function(count = NULL,          # count data
                  group = NULL,          # condition label
                  norm.factors = NULL,   # normalization factors
                  upper.limit = 0.4,      # maximum of seek range
                  sortby = NULL,         # sort by tissues
                  sort = FALSE,          # sort by statistics
                  deseq = 1
) {
    ## the function for finding specific tissues.
    .find_mincluster <- function(x) {
        d <- dist(x)
        m <- as.matrix(d)
        m[upper.tri(m)] <- NA
        diag(m) <- NA
        mincol <- (1:length(x))[x == min(x)][1]
        v1 <- m[mincol, !is.na(m[mincol, ])]
        v2 <- m[!is.na(m[, mincol]), mincol]
        vm <- mean(c(v1, v2))  
        vs <- sd(c(v1, v2))  
        l1 <- l2 <- NA
        if (mincol != 1)
            l1 <- seq(from = 1, to = length(v1), by = 1)[v1 < vm - 2 * vs]
        if (mincol != length(x))
            l2 <- seq(from = length(x), to = length(x) - length(v2) + 1,
                      by = -1)[v1 < vm - 2 * vs]
        v  <- rep(FALSE, length(x))
        v[c(l1, l2, mincol)] <- TRUE
        return (v)
    }

    ## arguments settings.
    n.libs <- ncol(count)
    n.tags <- nrow(count)
    if (is.null(group))
        group <- 1:n.libs

    if (unique(group) < 3) {
        stop("TCC::ERROR: YAYOI requires the data with more than three groups.")
    }


    ## normalize count data.
    if (is.null(norm.factors)) {
        tcc <- new("TCC", count, group)
        tcc <- calcNormFactors(tcc, norm.method = "deseq",
                               test.method = "yayoi",
                               iteration = 3,
                               deseq = deseq)
        norm.factors <- tcc$norm.factors
    }
    ef <- norm.factors * colSums(count)
    ef <- ef / mean(ef)

    ## estimate dispersions by DESeq.
    if (deseq == 1) {
        suppressMessages(d <-
               DESeq::newCountDataSet(countData = count, conditions = group)
             )
        suppressMessages(DESeq::sizeFactors(d) <- ef)
        suppressMessages(d <-
               DESeq::estimateDispersions(d, method = "blind",
                   sharingMode = "fit-only", fitType = "local")
             )
        phi <- d@featureData@data$disp_blind
    }
    if (deseq == 2) {
        group <- data.frame(group = group)
        suppressMessages(d <-
               DESeq2::DESeqDataSetFromMatrix(countData = count,
                     colData = group, design = ~ group)
             )
        suppressMessages(DESeq2::sizeFactors(d) <- ef)
        suppressMessages(d <- DESeq2::estimateDispersions(d))
        phi <- DESeq2::dispersions(d)
    }
    ## estimate means by DESeq approach.
    ##     if there more than one replicates for a tissues,
    ##     treate them as one replicate.
    ##     because we should create expression pattern.
    ugroup <- unique(group)
    nmc <- matrix(0, nrow = nrow(count), ncol = length(ugroup))
    colnames(nmc) <- ugroup
    for (i in 1:length(ugroup)) {
       nmc[, i] <- rowSums(as.matrix(count[, (group == ugroup[i])] / 
                           norm.factors[(group == ugroup[i])])) /
                   sum(group == ugroup[i])
    }
    nmcs <- t(apply(nmc, 1, sort))
    mu <- rowMeans(as.matrix(nmcs[, -c(1:(floor(upper.limit * n.libs / 2) + 1),
          rev(n.libs:(n.libs + 1 - (floor(upper.limit * n.libs / 2) + 1))))]))

    ## sort normalized count data for calcualting probability.
    prob <- pnbinom(nmc, mu = mu, size = 1 / phi)
    updw <- ifelse(prob > .5, 1, -1)          # 1:up, -1:down-regulated
    prob <- ifelse(prob > .5, 1 - prob, prob)
    prob[is.na(prob)] <- .5                   # no specific if error occurs
    yscore <- -log2(apply(prob, 1, min))
    yscore[is.infinite(yscore)] <- max(yscore[!is.infinite(yscore)]) + 0.1
    yscore[is.na(yscore)] <- 0                # no specific if error occurs

    ## create pattern matrix.
    expprof <- matrix(0, ncol = n.libs, nrow = n.tags)
    o <- matrix(0, ncol = length(ugroup), nrow = n.tags)
    for (i in 1:n.tags) {
        v <- .find_mincluster(prob[i, ])
        o[i, v] <- 1
        o[i, ] <- o[i, ] * updw[i, ]
    }
    for (i in 1:length(ugroup)) {
        expprof[, (group == ugroup[i])] <- o[, i]
    }
    rownames(expprof) <- rownames(count)
    colnames(expprof) <- colnames(count)

    if (sort) {
        if (!is.null(sortby)) {
            idx <- order(o[(ugroup == sortby)], yscore, decreasing = TRUE)
        } else {
            idx <- order(yscore, decreasing = TRUE)
        }
        expprof <- expprof[idx, ]
        yscore <- yscore[idx]
    }

    yscore[rowSums(abs(expprof)) == 0] <- min(yscore)
    return (list(outlier = expprof, score = yscore, prob = prob))
}
