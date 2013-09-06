test_calcNormFactors_DEGESedgeR_1 <- function() {
    data(hypoData)
    FDR <- 0.1
    floorPDEG <- 0.05
    count <- hypoData
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", count, group)

    tcc_edgeR <- calcNormFactors(tcc, norm.method = "tmm", iteration = 0)
    tcc_DEGES_edgeR <- calcNormFactors(tcc, norm.method = "tmm", 
                                       test.method = "edger", iteration = 1)
    d <- DGEList(counts = hypoData, group = group)
    d <- calcNormFactors(d)
    nf.1 <- d$samples$norm.factors
    nf.1 <- nf.1 / mean(nf.1)
    d <- estimateCommonDisp(d)
    d <- estimateTagwiseDisp(d)
    r <- exactTest(d)
    q <- p.adjust(r$table$PValue, method = "BH")
    if (sum(q < FDR) > (floorPDEG * nrow(hypoData))) {
        is.DEG <- as.logical(q < FDR)
    } else {
        is.DEG <- as.logical(rank(result$table$PValue, 
                             ties.method = "min") <=
                             nrow(hypoData) * floorPDEG)
    }
    d <- DGEList(counts = hypoData[!is.DEG, ], group = group)
    d <- calcNormFactors(d)
    nf.2 <- d$samples$norm.factors * colSums(hypoData[!is.DEG, ]) / 
                                           colSums(hypoData)
    nf.2 <- nf.2 / mean(nf.2)
    checkEqualsNumeric(nf.1, tcc_edgeR$norm.factors)
    checkEqualsNumeric(nf.2, tcc_DEGES_edgeR$norm.factors)
}

test_calcNormFactors_DEGESDESeq_1 <- function() {
    data(hypoData)
    FDR <- 0.1
    floorPDEG <- 0.05
    count <- hypoData
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", count, group)

    tcc_DESeq <- calcNormFactors(tcc, norm.method = "deseq", iteration = 0)
    tcc_DEGES_DESeq <- calcNormFactors(tcc, norm.method = "deseq", 
                                       test.method = "deseq", iteration = 1,
                                       FDR = FDR, floorPDEG = floorPDEG)
    cds <- newCountDataSet(hypoData, group)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    nf.1 <- sizeFactors(cds) / colSums(hypoData)
    nf.1 <- nf.1 / mean(nf.1)
    r <- nbinomTest(cds, 1, 2)
    r$pval[is.na(r$pval)] <- 1
    r$padj[is.na(r$padj)] <- 1
    if (sum(r$padj < FDR) > (floorPDEG * nrow(hypoData))) {
        is.DEG <- as.logical(r$padj < FDR)
    } else {
        is.DEG <- as.logical(rank(result$table$PValue, 
                             ties.method = "min") <=
                             nrow(hypoData) * floorPDEG)
    }
    cds <- newCountDataSet(hypoData[!is.DEG, ], group)
    cds <- estimateSizeFactors(cds)
    nf.2 <- sizeFactors(cds) / colSums(hypoData)
    nf.2 <- nf.2 / mean(nf.2)
    checkEqualsNumeric(nf.1, tcc_DESeq$norm.factors)
    checkEqualsNumeric(nf.2, tcc_DEGES_DESeq$norm.factors)
}

test_calcNormFactors_DEGESTbT <- function() {
    data(hypoData)
    count <- hypoData
    group <- c(1, 1, 1, 2, 2, 2)
    tcc <- new("TCC", count, group)
    set.seed(1)
    tcc_DEGES_baySeq <- calcNormFactors(tcc, norm.method = "tmm", 
                                        test.method = "bayseq",
                                        iteration = 1, samplesize = 10)
    d <- DGEList(count = hypoData, group = group)
    d <- calcNormFactors(d)
    nf.1 <- d$samples$norm.factors
    nf.1 <- nf.1 / mean(nf.1)
    cD <- new("countData", data = hypoData, replicates = group,
             groups = list(NDE = rep(1, length = length(group)), DE = group),
             libsizes = colSums(hypoData) * nf.1)
    set.seed(1)
    cD <- getPriors.NB(cD, samplesize = 10, estimation = "QL", cl = NULL)
    cD <- getLikelihoods.NB(cD, pET = "BIC", cl = NULL)
    is.DEG <- as.logical(rank(-cD@posteriors[, "DE"]) <
                         (nrow(hypoData) * cD@estProps[2]))
    d <- DGEList(count = hypoData[!is.DEG, ], group = group)
    d <- calcNormFactors(d)
    nf.2 <- d$samples$norm.factors * colSums(hypoData[!is.DEG, ]) /
                    colSums(hypoData)
    nf.2 <- nf.2 / mean(nf.2)

    checkEqualsNumeric(nf.2, tcc_DEGES_baySeq$norm.factors)
}

test_calcNormFactors_DEGESDESeq_classic_single <- function() {
    data(hypoData)
    floorPDEG <- 0.05
    count <- hypoData[, c(1, 4)]
    group <- c(1, 2)
    tcc <- new("TCC", count, group)

    tcc <- calcNormFactors(tcc, norm.method = "tmm", iteration = 0)
    tcc <- calcNormFactors(tcc, norm.method = "deseq", iteration = 0)
    tcc <- calcNormFactors(tcc, norm.method = "deseq", 
                           test.method = "deseq", iteration = 1)
    tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq",
                           iteration = 1, samplesize = 10)
}


test_calcNormFactors_DEGESedgeR_glm <- function() {
    data(hypoData_mg)
    count <- hypoData_mg
    group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    tcc <- new("TCC", count, group)
    
    desgin <- model.matrix(~ 0 + factor(group))
    coef <- 2:3
    tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                           design = desgin, iteration = 1)
}


test_calcNormFactors_DEGESDESeq_glm <- function() {
    data(hypoData_mg)
    count <- hypoData_mg
    group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    tcc <- new("TCC", count, group)
    
    fit1 <- count ~ condition
    fit0 <- count ~ 1
    tcc <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                           fit1 = fit1, fit0 = fit0, iteration = 1)
}

test_calcNormFactors_increment <- function() {
  data(hypoData)
  tcc <- new("TCC", hypoData, c(1, 1, 1, 2, 2, 2))
  tcc.0 <- calcNormFactors(tcc, iteration = 0)
  tcc.1 <- calcNormFactors(tcc, iteration = 1)
  tcc.0.1 <- calcNormFactors(tcc, increment = TRUE)
  checkEqualsNumeric(tcc.1$norm.factors, tcc.0.1$norm.factors)


  tcc.3 <- calcNormFactors(tcc, iteration = 3)
  tcc.1 <- calcNormFactors(tcc, increment = TRUE)
  tcc.1.1 <- calcNormFactors(tcc.1, increment = TRUE)
  tcc.1.1.1 <- calcNormFactors(tcc.1.1, increment = TRUE)
  checkEqualsNumeric(tcc.3$norm.factors, tcc.1.1.1$norm.factors)


  tcc.1 <- calcNormFactors(tcc, iteration = 1)
  tcc.1.2 <- calcNormFactors(tcc.1, iteration = 2, increment = TRUE)
  checkEqualsNumeric(tcc.3$norm.factors, tcc.1.2$norm.factors)
}



