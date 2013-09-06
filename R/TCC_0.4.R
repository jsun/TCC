MAplot <- function(datalist, FDR_threshold = 0.01){
    data <- datalist$counts
    data.cl <- datalist$group
    norm_f_TbT <- datalist$norm_f_TbT
    x_axis <- datalist$Mval
    y_axis <- datalist$Aval
    plot(x_axis, y_axis, xlab = "A = (log2(B)+log2(A))/2", 
    ylab = "M = log2(B)-log2(A)", pch = 20, cex = .3)
    grid(col = "gray", lty = "dotted")
    points(x_axis[datalist$data$FDR < FDR_threshold], 
           y_axis[datalist$data$FDR < FDR_threshold], 
           col = 2, pch = 20, cex = 0.3)
    baseline_TbT <- log2(mean(norm_f_TbT[data.cl == 2]) / 
                         mean(norm_f_TbT[data.cl == 1]))
    abline(h = baseline_TbT, col = "red", lwd = 1)
}



## generate negative binomial distributed datasets with different frequencies
NBsample <- function(DEG_foldchange = 4, repA = 3, repB = 3, 
                     Ngene = 3000, PDEG = 0.15, PA = 0.2){
    arab <- NULL;rm(arab)  # to avoid note by R CMD check
    data(arab)             # arab dataset from NBPseq
    data.cl <- c(rep(1, 3), rep(2, 3))
    RPM <- sweep(arab, 2, 1000000 / colSums(arab), "*")
    RPM_A <- RPM[,data.cl == 1]
    RPM_B <- RPM[,data.cl == 2]
    RPM_A <- RPM_A[apply(RPM_A, 1, var) > 0,]
    RPM_B <- RPM_B[apply(RPM_B, 1, var) > 0,]
    MEAN <- c(apply(RPM_A, 1, mean), apply(RPM_B, 1, mean))
    VARIANCE <- c(apply(RPM_A, 1, var), apply(RPM_B, 1, var))
    DISPERSION <- (VARIANCE - MEAN) / (MEAN * MEAN)
    mean_disp_tmp <- cbind(MEAN, DISPERSION)
    mean_disp_tmp <- mean_disp_tmp[mean_disp_tmp[,2] > 0,]
    resampling_vector <- sample(1:nrow(mean_disp_tmp), Ngene, replace = TRUE)
    mean_disp <- mean_disp_tmp[resampling_vector,]
    mu <- mean_disp[,1]
    DEG_degree_A <- rep(1, Ngene)
    DEG_degree_A[1:(Ngene*PDEG*PA)] <- DEG_foldchange
    mu_A <- mu * DEG_degree_A
    DEG_degree_B <- rep(1, Ngene)
    DEG_degree_B[(Ngene * PDEG * PA + 1):(Ngene * PDEG)] <- DEG_foldchange
    mu_B <- mu * DEG_degree_B
    DEG_posi_org <- (DEG_degree_A * DEG_degree_B) > 1
    nonDEG_posi_org <- (DEG_degree_A * DEG_degree_B) == 1
    outA <- NULL
    colnamev <-NULL
    for(i in 1:repA){
        outA <- cbind(outA, rnbinom(n = length(mu_A), 
                      mu = mu_A, size = 1 / mean_disp[,2]))
        colnamev <-cbind(colnamev, paste("A", as.character(i), sep = ""))
    }
    outB <- NULL
    for(i in 1:repB){
        outB <- cbind(outB, rnbinom(n = length(mu_B), 
                      mu = mu_B, size = 1 / mean_disp[,2]))
        colnamev <-cbind(colnamev, paste("B", as.character(i), sep = ""))
    }
    out <- cbind(outA, outB)
    colnames(out) <- colnamev
    obj <- rowSums(out) > 0
    RAW <- out[obj,]
    DEG_posi <- DEG_posi_org[obj]
    nonDEG_posi <- nonDEG_posi_org[obj]
    retval <- list(RAW, DEG_posi, nonDEG_posi)
    names(retval) <- c("data", "DEG_posi", "nonDEG_posi")
    return(retval)
}



##  TbT  normalization methods
do_TbT <- function(data, data.cl, sample_num = 10000){
    RAW <- data
    ##  Step 1: first normalization 
    d <- DGEList(counts = data, group = data.cl)
    d <- calcNormFactors(d)
    norm_f_TMM <- d$samples$norm.factors
    names(norm_f_TMM) <- colnames(data)

    ##  Step 2: DEG identification 
    groups <- list(NDE = rep(1, length(data.cl)), DE = data.cl)
    norm_f_RPM = 1000000 / colSums(data)
    RPM <- sweep(data, 2, norm_f_RPM, "*")
    data <- round(RPM)
    once_normalized <- new("countData", data = as.matrix(data), 
                           replicates = data.cl, 
                           libsizes = colSums(data) * norm_f_TMM, 
                           groups = groups)
    once_normalized.NB <- getPriors.NB(once_normalized, 
                                       samplesize = sample_num, 
                                       estimation = "QL", cl = NULL)
    out <- getLikelihoods.NB(once_normalized.NB, pET = "BIC", cl = NULL)
    PDEG <- out@estProps[2]     # proportion of differentially expressed genes
    rank_bayseq <- rank(-out@posteriors[,2])
    NDEG <- (nrow(data) * PDEG) # number of differentially expressed genes

    ##  Step 3: second normalization 
    obj_DEGy <- (rank_bayseq < NDEG)
    obj_DEGn <- (rank_bayseq >= NDEG)
    data <- RAW[obj_DEGn,]
    d <- DGEList(counts = data, group = data.cl)
    d <- calcNormFactors(d)
    norm_f_TbTorg <- d$samples$norm.factors * colSums(data) / colSums(RAW)
    norm_f_TbT <- norm_f_TbTorg / mean(c(mean(norm_f_TbTorg[data.cl == 1]),
                                         mean(norm_f_TbTorg[data.cl == 2])))
    data <- RPM
    meanA <- log2(apply(data[,data.cl == 1], 1, mean))
    meanB <- log2(apply(data[,data.cl == 2], 1, mean))
    Aval <- (meanA + meanB) / 2
    Mval <- meanB - meanA

    ##  calculation of PA value (degree of biased expression)  ###
    RPM_TMM <- sweep(RPM, 2, 1 / norm_f_TMM, "*")
    data <- RPM_TMM
    logratio <- log2(apply(data[,data.cl == 2], 1, mean)) - 
                log2(apply(data[,data.cl == 1], 1, mean))
    PA <- sum(logratio[rank_bayseq < NDEG] < 0) / NDEG

    retval <- list(norm_f_TbT, Aval, Mval, PDEG, PA, obj_DEGn, obj_DEGy, norm_f_TMM, norm_f_TbTorg, data.cl, data)
    names(retval) <- c("norm_f_TbT", "Mval", "Aval", "PDEG", "PA", "nonDEG_posi", "DEG_posi", "norm_f_TMM", "norm_f_TbTorg", "data.cl", "data")
    return(retval)
}


exactTestafterTbT <- function(names, counts, group, sample_num = 10000){
    ##if (!("edgeR" %in% loadedNamespaces()))
    ##  library(edgeR)
    ##edgeR_Version <-sessionInfo()$otherPkgs$edgeR$Version
    ##edgeR_v <- as.integer(strsplit(edgeR_Version, '.', fixed=TRUE)[[1]])
    tbtout <- do_TbT(counts, group, sample_num)
    d <- DGEList(counts = counts, group = group)
    d$samples$norm.factors <- tbtout$norm_f_TbT 
    d <- estimateCommonDisp(d)

    ##if (edgeR_v[[1]] == 2 & edgeR_v[[2]] <= 6) {
    ##if(is.null(span) == FALSE) prop.used <- 
    ## span else if(is.null(prop.used)) prop.used <- 0.5
    d <- estimateTagwiseDisp(d) 
    ##}
    ##if (edgeR_v[[1]] > 2 | edgeR_v[[2]] >= 7)
    ##  if(is.null(prop.used) == FALSE) span <- prop.used
    ##  d <- estimateTagwiseDisp(d, span=span, grid.length=grid.length) 
    ##}
    out <- exactTest(d)
    if(is.vector(out$table$PValue)){        # for current edgeR
        FDR <- p.adjust(out$table$PValue, method = "BH")
    }else if(is.vector(out$table$p.value)){ # for older edgeR
        FDR <- p.adjust(out$table$p.value, method = "BH")
    }else{                                  # something strange
        warning("PValue was not available")
    }
    rank_edgeR <- rank(FDR)
    retval <- cbind(names, out$table, FDR, rank_edgeR) 
    return(list(data = retval, norm_f_TbT = tbtout$norm_f_TbT, 
                Mval = tbtout$Mval, Aval = tbtout$Aval, 
                counts = counts,  group = group))
}

