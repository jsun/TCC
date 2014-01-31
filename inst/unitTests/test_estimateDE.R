test_estimateDE_EBSeq_1 <- function() {
    library(EBSeq)
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3))
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    set.seed(1)
    tcc <- estimateDE(tcc, test.method = "ebseq", samplesize = 10)
    auc <- calcAUCValue(tcc)
 
    set.seed(1)
    x <- EBTest(Data = tcc$count,
                Conditions = as.factor(tcc$group[, 1]),
                sizeFactors = tcc$norm.factors * colSums(tcc$count),
                maxround = 10)
    PP <- GetPPMat(x)
    df <- matrix(1, ncol = 2, nrow = nrow(tcc$count))
    rownames(df) <- rownames(tcc$count)
    df[rownames(PP), 1] <- PP[, 1]
    df[rownames(PP), 2] <- PP[, 2]
    df[is.na(df)] <- 0

    checkEqualsNumeric(df[, 2], tcc$stat$prob)
    checkTrue(auc > 0.80)
}

test_estimateDE_EBSeq_2 <- function() {
    library(EBSeq)
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3, 3))
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    set.seed(1)
    tcc <- estimateDE(tcc, test.method = "ebseq", samplesize = 10)
    auc <- calcAUCValue(tcc)

    g <- tcc$group[, 1]
    ug <- unique(g)
    gp <- matrix(c(rep(1, length = length(ug)), 1:length(ug)),
                 nrow = 2, byrow = TRUE)
    colnames(gp) <- ug
    rownames(gp) <- c("Pattern1", "Pattern2")
    set.seed(1)
    x <- EBMultiTest(Data = tcc$count,
                NgVector = NULL,
                Conditions = g,
                AllParti = gp,
                sizeFactors = tcc$norm.factors * colSums(tcc$count),
                maxround = 10)
    PP <- GetMultiPP(x)
    df <- matrix(1, ncol = 2, nrow = nrow(tcc$count))
    rownames(df) <- rownames(tcc$count)
    df[rownames(PP$PP), 1] <- PP$PP[, 1]
    df[rownames(PP$PP), 2] <- PP$PP[, 2]
    df[is.na(df)] <- 0

    checkEqualsNumeric(df[, 2], tcc$stat$prob)
    checkTrue(auc > 0.80)
}


test_estimateDE_SAMseq_1 <- function() {
    library(samr)
    samplesize <- 10
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3))
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    set.seed(1)
    tcc <- estimateDE(tcc, test.method = "samseq", samplesize = samplesize)
    auc <- calcAUCValue(tcc)

    x <- round(getNormalizedData(tcc))
    set.seed(1)
    d <- SAMseq(x = x, y = tcc$group[, 1],
                resp.type = "Two class unpaired",
                nperms = samplesize)

    checkEqualsNumeric(d$samr.obj$tt, tcc$stat$testStat)
    checkTrue(auc > 0.80)
}

##test_estimateDE_SAMseq_1p <- function() {
##    library(samr)
##    samplesize <- 10
##    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3))
##    tcc <- calcNormFactors(tcc, iteration = FALSE)
##    set.seed(1)
##    tcc <- estimateDE(tcc, test.method = "samseq", paired = TRUE,
##                      samplesize = samplesize)
##    auc <- calcAUCValue(tcc)
##
##    x <- round(getNormalizedData(tcc))
##    set.seed(1)
##    d <- SAMseq(x = x, y = tcc$group[, 1],
##                resp.type = "Two class paired",
##                nperms = samplesize)
##
##    checkEqualsNumeric(d$samr.obj$tt, tcc$stat$testStat)
##    checkTrue(auc > 0.80)
##}

test_estimateDE_SAMseq_2 <- function() {
    library(samr)
    samplesize <- 10
    tcc <- simulateReadCounts(replicates = c(3, 3, 3)) 
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    set.seed(1)
    tcc <- estimateDE(tcc, test.method = "samseq", samplesize = samplesize)
    auc <- calcAUCValue(tcc)

    x <- round(getNormalizedData(tcc))
    set.seed(1)
    d <- SAMseq(x = x, y = tcc$group[, 1], 
                resp.type = "Multiclass",
                nperms = samplesize)
    checkEqualsNumeric(d$samr.obj$tt, tcc$stat$testStat)
    checkTrue(auc > 0.80)
}

##test_estimateDE_NOISeq_1 <- function() {
##    library(NOISeq)
##    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3))
##    tcc <- calcNormFactors(tcc, iteration = FALSE)
##    tcc <- estimateDE(tcc, test.method = "noiseq")
##    auc <- calcAUCValue(tcc)
##
##    x <- getNormalizedData(tcc)
##    gl <- data.frame(group = tcc$group[, 1])
##    nd <- NOISeq::readData(x, gl)
##    nr <- noiseq(nd, k = 0.5, norm = "n", replicates = "biological",
##                 factor = "group", conditions = unique(tcc$group[, 1]))
##    prob <- nr@results[[1]]$prob
##    prob[is.na(prob)] <- 0
##    checkEqualsNumeric(prob, tcc$stat$prob)
##    checkTrue(auc > 0.80)
##}

test_estimateDE_baySeq_1 <- function() {
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3))
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    set.seed(1)
    tcc <- estimateDE(tcc, test.method = "bayseq", samplesize = 10)
    auc <- calcAUCValue(tcc)

    group <- c(1, 1, 1, 2, 2, 2)
    el <- colSums(tcc$count) * tcc$norm.factors
    groups <- list(NDE = rep(1, length(group)), DE = group)
    cD <- new("countData", data = tcc$count, replicates = group,
              libsizes = colSums(tcc$count) * tcc$norm.factors, 
              groups = groups)
    set.seed(1)
    cD <- getPriors.NB(cD, samplesize = 10,
                       estimation = "QL", cl = NULL)
    cD <- getLikelihoods.NB(cD, pET = "BIC", cl = NULL)
    tmp <- topCounts(cD, group = "DE", number = nrow(tcc$count))
    tmp <- tmp[rownames(tcc$count), ]
    p <- 1 - tmp$Likelihood

    checkEqualsNumeric(p, tcc$stat$p.value)
    checkTrue(auc > 0.60)
}

##test_estimateDE_baySeq_1p <- function() {
##    tcc <- simulateReadCounts(Ngene = 1000, PDEG = 0.3,
##               group = data.frame(A = c(1, 1, 1, 1, 2, 2, 2, 2),
##                                  B = c(1, 1, 2, 2, 1, 1, 2, 2)),
##               DEG.foldchange = data.frame(F1 = c(4, 4, 4, 4, 1, 1, 1, 1),
##                                           F2 = c(1, 1, 1, 1, 4, 4, 4, 4),
##                                           F3 = c(1, 1, 1/4, 1/4, 1, 1, 4, 4)),
##               DEG.assign = c(0.2, 0.2, 0.6))
##    tcc <- calcNormFactors(tcc, iteration = FALSE)
##    set.seed(1)
##    tcc <- estimateDE(tcc, test.method = "bayseq", paired = TRUE,
##                      samplesize = 10)
##    auc <- calcAUCValue(tcc)
##
##    group <- c(1, 1, 2, 2)
##    el <- colSums(tcc$count) * tcc$norm.factors
##    groups <- list(NDEG = c(1, 1, 1, 1), DE = c(1, 1, 2, 2))
##    cD <- new("pairedData", 
##               data = tcc$count[, 1:4],
##               pairData = tcc$count[, 5:8],
##               replicates = group,
##               groups = groups,
##               libsizes = el[1:4],
##               pairLibsizes = el[5:8])
##    set.seed(1)
##    cD <- getPriors.BB(cD, samplesize = 300,
##                       estimation = "QL", cl = NULL)
##    cD <- getLikelihoods.BB(cD, pET = "BIC", nullProps = 0.5, cl = NULL)
##
##    ## DE between replicate groups    
##    tmp <- topCounts(cD, group = 1, number = nrow(tcc$count))
##    tmp <- tmp[rownames(tcc$count), ]
##    p <- 1 - tmp$Likelihood
##AUC(rocdemo.sca(truth=c(rep(0, 2700),rep(1, 300),rep(0, 7000)), data = -rank(p)))
##AUC(rocdemo.sca(truth=c(rep(1, 2700),rep(0, 300),rep(0, 7000)), data = -rank(p)))
##
##    tmp <- topCounts(cD, group = 2, number = nrow(tcc$count))
##    tmp <- tmp[rownames(tcc$count), ]
##    p <- 1 - tmp$Likelihood
##AUC(rocdemo.sca(truth=c(rep(0, 2700),rep(1, 300),rep(0, 7000)), data = -rank(p)))
##AUC(rocdemo.sca(truth=c(rep(1, 2700),rep(0, 300),rep(0, 7000)), data = -rank(p)))
##
##
##    checkEqualsNumeric(p, tcc$stat$p.value)
##    checkTrue(auc > 0.70)
##}

test_estimateDE_baySeq_2 <- function() {
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3, 3))
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    set.seed(1)
    tcc <- estimateDE(tcc, test.method = "bayseq", samplesize = 10)
    auc <- calcAUCValue(tcc)

    group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    el <- colSums(tcc$count) * tcc$norm.factors
    groups <- list(NDE = rep(1, length(group)), DE = group)
    cD <- new("countData", data = tcc$count, replicates = group,
              libsizes = colSums(tcc$count) * tcc$norm.factors, 
              groups = groups)
    set.seed(1)
    cD <- getPriors.NB(cD, samplesize = 10,
                       estimation = "QL", cl = NULL)
    cD <- getLikelihoods.NB(cD, pET = "BIC", cl = NULL)
    tmp <- topCounts(cD, group = "DE", number = nrow(tcc$count))
    tmp <- tmp[rownames(tcc$count), ]
    p <- 1 - tmp$Likelihood

    checkEqualsNumeric(p, tcc$stat$p.value)
    checkTrue(auc > 0.70)
}

test_estimateDE_baySeq_3 <- function() {
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(4, 4))
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    tcc$group <- data.frame(GROUP = c(1, 1, 1, 1, 2, 2, 2, 2 ),
                            TIME = c(1, 1, 2, 2, 1, 1, 2, 2))
    set.seed(1)
    tcc <- estimateDE(tcc, test.method = "bayseq", samplesize = 10,
                      comparison = "GROUP")
    auc <- calcAUCValue(tcc)

    groups <- tcc$group
    groups <- cbind(rep(1, length = nrow(tcc$group)), groups)
    colnames(groups)[1] <- "NDE"
    el <- colSums(tcc$count) * tcc$norm.factors
    cD <- new("countData", data = tcc$count,
              replicates = c(1, 1, 1, 1, 2, 2, 2, 2),
              libsizes = colSums(tcc$count) * tcc$norm.factors, 
              groups = groups)
    set.seed(1)
    cD <- getPriors.NB(cD, samplesize = 10,
                       estimation = "QL", cl = NULL)
    cD <- getLikelihoods.NB(cD, pET = "BIC", cl = NULL)
    tmp <- topCounts(cD, group = "GROUP", number = nrow(tcc$count))
    tmp <- tmp[rownames(tcc$count), ]
    p <- 1 - tmp$Likelihood

    checkEqualsNumeric(p, tcc$stat$p.value)
    checkTrue(auc > 0.70)
}


test_estimateDE_DESeq_1 <- function() {
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3)) 
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    tcc <- estimateDE(tcc, test.method = "deseq")
    auc <- calcAUCValue(tcc)

    d <- newCountDataSet(tcc$count, tcc$group[, 1])
    sizeFactors(d) <- tcc$norm.factors * colSums(tcc$count)
    d <- estimateDispersions(d)
    r <- nbinomTest(d, 1, 2)
    r$pval[is.na(r$pval)] <- 1

    checkEqualsNumeric(r$pval, tcc$stat$p.value)
    checkTrue(auc > 0.80)
}

test_estimateDE_DESeq_2 <- function() {
    fit1 <- count ~ condition
    fit0 <- count ~ 1

    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3, 3)) 
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    t1 <- estimateDE(tcc, test.method = "deseq")
    t2 <- estimateDE(tcc, test.method = "deseq", fit0 = fit0, fit1 = fit1)
    auc <- calcAUCValue(t1)
    
    d <- newCountDataSet(tcc$count, tcc$group[, 1])
    sizeFactors(d) <- tcc$norm.factors * colSums(tcc$count)
    d <- estimateDispersions(d)
    f0 <- fitNbinomGLMs(d, fit0)
    f1 <- fitNbinomGLMs(d, fit1)
    p <- nbinomGLMTest(f1, f0)
    p[is.na(p)] <- 1
    
    checkEqualsNumeric(p, t1$stat$p.value)
    checkEqualsNumeric(p, t2$stat$p.value)
    checkTrue(auc > 0.80)
}

test_estimateDE_DESeq_3 <- function() {
    group <- data.frame(
               COND = as.factor(c(1, 1, 1, 1, 2, 2, 2, 2)),
               TIME = as.factor(c(1, 1, 2, 2, 1, 1, 2, 2))
             )
    fit1 <- count ~ TIME + COND
    fit0 <- count ~ 1

    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(4, 4)) 
    tcc$group <- group
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    tcc <- estimateDE(tcc, test.method = "deseq", fit0 = fit0, fit1 = fit1)
    auc <- calcAUCValue(tcc)
    
    d <- newCountDataSet(tcc$count, tcc$group)
    sizeFactors(d) <- tcc$norm.factors * colSums(tcc$count)
    d <- estimateDispersions(d)
    f0 <- fitNbinomGLMs(d, fit0)
    f1 <- fitNbinomGLMs(d, fit1)
    p <- nbinomGLMTest(f1, f0)
    p[is.na(p)] <- 1
    
    checkEqualsNumeric(p, tcc$stat$p.value)
    checkTrue(auc > 0.80)
}


test_estimateDE_edgeR_1 <- function() {
    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3)) 
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    tcc <- estimateDE(tcc, test.method = "edger")
    auc <- calcAUCValue(tcc)

    d <- DGEList(counts = tcc$count, group = tcc$group[, 1])
    d$samples$norm.factors <- tcc$norm.factors
    d <- estimateCommonDisp(d)
    d <- estimateTagwiseDisp(d)
    r <- exactTest(d)
    checkEqualsNumeric(r$table$PValue, tcc$stat$p.value)
    checkTrue(auc > 0.80)
}

test_estimateDE_edgeR_2 <- function() {
    group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    coef <- 2:3
    design <- model.matrix(~ as.factor(group))

    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3, 3)) 
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    t1 <- estimateDE(tcc, test.method = "edger")
    t2 <- estimateDE(tcc, test.method = "edger", design = design, coef = coef)
    auc <- calcAUCValue(t1)
 
    d <- DGEList(counts = tcc$count, group = tcc$group[, 1])
    d$samples$norm.factors <- tcc$norm.factors
    d <- estimateGLMCommonDisp(d, design)
    d <- estimateGLMTrendedDisp(d, design)
    d <- estimateGLMTagwiseDisp(d, design)
    fit <- glmFit(d, design)
    lrt <- glmLRT(fit, coef = coef)
    r <- topTags(lrt, n = nrow(tcc$count))
    r <- r$table[rownames(tcc$count), ]

    checkEqualsNumeric(r$PValue, t1$stat$p.value)
    checkEqualsNumeric(r$PValue, t2$stat$p.value)
    checkTrue(auc > 0.80)
}

test_estimateDE_edgeR_3 <- function() {
    group <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    contrast <- c(-1, 0, 1)
    design <- model.matrix(~ 0 + as.factor(group))

    tcc <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3, 3)) 
    tcc <- calcNormFactors(tcc, iteration = FALSE)
    tcc <- estimateDE(tcc, test.method = "edger", design = design, contrast = contrast)
    auc <- calcAUCValue(tcc)
 
    d <- DGEList(counts = tcc$count, group = tcc$group[, 1])
    d$samples$norm.factors <- tcc$norm.factors
    d <- estimateGLMCommonDisp(d, design)
    d <- estimateGLMTrendedDisp(d, design)
    d <- estimateGLMTagwiseDisp(d, design)
    fit <- glmFit(d, design)
    lrt <- glmLRT(fit, contrast = contrast)
    r <- topTags(lrt, n = nrow(tcc$count))
    r <- r$table[rownames(tcc$count), ]

    checkEqualsNumeric(r$PValue, tcc$stat$p.value)
    checkTrue(auc > 0.80)
}

test_estimateDE_crossvalidate <- function() {
    tcc <- new("TCC")
    av <- tcc$private$available$test.method
    ty <- colnames(av)
    pk <- rownames(av)
    for (i in 1:length(ty)) {
        for (j in 1:length(pk)) {
            if (av[j, i]) {
                if (ty[i] == "UnRepTwoGroup")
                    x <- simulateReadCounts(Ngene = 1000, replicates = c(1, 1))
                else if (ty[i] == "TwoGroup" || ty[i] == "PairedTwoGroup")
                    x <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3))
                else 
                    x <- simulateReadCounts(Ngene = 1000, replicates = c(3, 3, 3))
                x <- calcNormFactors(x, norm.method = "tmm",
                                     test.method = pk[j], samplesize = 10)
                x <- estimateDE(x, test.method = pk[j], samplesize = 10)
            }
        }
    }
}




