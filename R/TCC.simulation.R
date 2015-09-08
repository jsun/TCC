makeFCMatrix <- function(Ngene = 10000, PDEG = 0.20, DEG.assign = NULL, replicates = NULL, fc.params = NULL) {
    if (is.null(DEG.assign)) DEG.assign <- c(0.8, 0.2)
    if (is.null(replicates)) replicates <- c(3, 3)

    ## intialise fc parameters
    if (is.null(fc.params)) {
        fc.params <- matrix(0, nrow = length(replicates), ncol = 3)
        colnames(fc.params) <- c("floor", "shape", "scale")
        for (i in 1:nrow(fc.params)) fc.params[i, ] <- c(1.2, 2, 0.5)
    }

    ## create single replicate fold-change matrix
    fc.matrix.single <- matrix(1, ncol = length(replicates), nrow = Ngene)
    fc.index <- c(0, cumsum(round(Ngene * PDEG * DEG.assign)))
    for (i in 2:length(fc.index)) {
        if (fc.index[i] - fc.index[i - 1] == 0) next
        fc.index.range.i <- (fc.index[i - 1] + 1):(fc.index[i])
        fc.matrix.single[fc.index.range.i, i - 1] <- fc.params[i - 1, 1] +
            rgamma(length(fc.index.range.i), shape = fc.params[i - 1, 2], scale = fc.params[i - 1, 3])
    }

    ## replicate
    fc.matrix <- matrix(1, ncol = sum(replicates), nrow = Ngene)
    k <- 0
    for (i in 1:length(replicates)) {
        for (j in 1:replicates[i]) {
            k <- k + 1
            fc.matrix[, k] <- fc.matrix.single[, i]
        }
    }

    message("Fold change matrix has created. Please use the same parameters in 'simulateReadCounts' function to simulate data")
    message(paste0("   Ngene <- ", Ngene))
    message(paste0("   PDEG  <- ", PDEG))
    message(paste0("   DEG.assign <- c(", paste(DEG.assign, collapse = ","), ")"))
    message(paste0("   replicates <- c(", paste(replicates, collapse = ","), ")"))
    message("Example:")
    message(paste0("simulateReadCounts(",
                paste0("Ngene = ", Ngene), ", ",
                paste0("PDEG = ", PDEG), ",",
                paste0("DEG.assign = c(", paste(DEG.assign, collapse = ", "), ")"), ",",
                paste0("replicates = c(", paste(replicates, collapse = ", "), ")"), ",",
                "fc.matrix = fc.matrix)"))

    fc.matrix
}



simulateReadCounts <- function(Ngene = 10000, PDEG = 0.20,
                               DEG.assign = NULL, DEG.foldchange = NULL,
                               replicates = NULL, group = NULL, fc.matrix = NULL) {

    ## initialize the NULL arguments
    if (is.null(group)) {
        ## replicates
        ## [1] 3 3
        if (is.null(replicates))
            replicates <- c(3, 3)

        ## cond.num
        ## [1] 2
        cond.num <- length(replicates)

        ## DEG.assign
        ## [1] 0.9 0.1
        if (is.null(DEG.assign))
            DEG.assign <- c(0.9, rep(0.1 / (cond.num - 1),
                                     length = cond.num - 1))

        ## DEG.foldchange
        ## [1] 4  4
        if (is.null(DEG.foldchange))
            DEG.foldchange <- rep(4, length = cond.num)

        ## group, DEG.fc
        ##   V1  V2
        ## 1  1   1
        ## 2  1   1
        ## 3  1   1
        ## 4  1   1
        ## 5  1   1
        ## 6  1   1
        group <- as.data.frame(matrix(1, nrow = sum(replicates), 
                                         ncol = cond.num))
        DEG.fc <- as.data.frame(matrix(1, nrow = sum(replicates),
                                          ncol = cond.num))

        ## reps
        ## [1] 1 1 1 2 2 2
        reps <- rep(1:cond.num, times = replicates)

        ## group, DEG.fc, DEG.foldchange
        ##   V1 V2
        ## 1  4  1
        ## 2  4  1
        ## 3  4  1
        ## 4  1  4
        ## 5  1  4
        ## 6  1  4
        for (i in 1:cond.num) {
            group[(reps == i), i] <- 2
            DEG.fc[(reps == i), i] <- DEG.foldchange[i]
        }
        DEG.foldchange <- DEG.fc
    }


    ## check the arguments
    if (is.null(group))
        stop("TCC::ERROR: The 'group' argument is required.")
    if (!is.data.frame(group))
        stop("TCC::ERROR: The 'group' argument should be data.frame.")
    if (is.null(fc.matrix)) {
        if (is.null(DEG.assign))
            stop("TCC::ERROR: The 'DEG.assign' argument is required.")
        if (is.null(DEG.foldchange))
            stop("TCC::ERROR: The 'DEG.foldchange' argument is required.")
        if (sum(DEG.assign) > 1)
            stop("TCC::ERROR: The total value of DEG.assign must less than one.")
        if (!is.data.frame(DEG.foldchange))
            stop("TCC::ERROR: The 'DEG.foldchange' argument should be data.frame.")
        if (nrow(group) != nrow(DEG.foldchange))
            stop("TCC::ERROR: The number of rows of 'group' and 'DEG.foldchange' must equal.")
        if (length(DEG.assign) != ncol(DEG.foldchange))
            stop("TCC::ERROR: The length of 'DEG.assign' should equal to the number of columns of 'DEG.foldchange'.")
    }

    ## print the simulation conditions
    message("TCC::INFO: Generating simulation data under NB distribution ...")
    message(paste("TCC::INFO: (genesizes   : ", Ngene, ")"))
    if (is.null(fc.matrix)) {
        if (!is.null(replicates)) {
            message(paste("TCC::INFO: (replicates  : ", paste(replicates, collapse=", "), ")"))
            message(paste("TCC::INFO: (PDEG        : ", paste(PDEG * DEG.assign, collapse=", "), ")"))
        } else {
            message(paste("TCC::INFO: (samples     : ", nrow(group), ")"))
            message(paste("TCC::INFO: (factors     : ", ncol(group), ")"))
            message(paste("TCC::INFO: (PDEG        : ", PDEG, ")"))
        }
    } else {
            message(paste("TCC::INFO: (PDEG        : ", sum(rowSums(fc.matrix) != ncol(fc.matrix)), ")"))
    }



    ## create foldchange matrix for sampling count data with foldchange
    ## fc.matrix
    ##      [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,]    4    4    4    1    1    1
    ## [2,]    4    4    4    1    1    1
    ## [3,]    4    4    4    1    1    1
    ## [4,]    4    4    4    1    1    1
    ## [5,]    4    4    4    1    1    1
    ## [6,]    4    4    4    1    1    1
    if (is.null(fc.matrix)) {
        ## create foldchange matrix for sampling count data with foldchange
        fc.matrix <- matrix(1, nrow = Ngene, ncol = nrow(group))
        fc.index  <- c(0, cumsum(round(Ngene * PDEG * DEG.assign)))
        for (i in 2:length(fc.index)) {
            if (fc.index[i] - fc.index[i - 1] == 0) next
            fc.matrix[(fc.index[i - 1] + 1):(fc.index[i]), ] <-
                fc.matrix[(fc.index[i - 1] + 1):(fc.index[i]), ] *
                matrix(rep(DEG.foldchange[, i - 1],
                           times = fc.index[i] - fc.index[i - 1]),
                       ncol = ncol(fc.matrix), byrow = TRUE)
        }
    }
    trueDEG <- as.numeric(rowSums(abs(fc.matrix)) != ncol(fc.matrix))

    ## prepare the population ('arab') for sampling
    arab <- NULL
    rm(arab)
    data(arab)
    rpm.a <- sweep(arab[, 1:3], 2, 
                   median(colSums(arab[, 1:3])) / colSums(arab[, 1:3]), "*")
    rpm.b <- sweep(arab[, 4:6], 2, 
                   median(colSums(arab[, 4:6])) / colSums(arab[, 4:6]), "*")
    rpm.a <- rpm.a[apply(rpm.a, 1, var) > 0, ]
    rpm.b <- rpm.b[apply(rpm.b, 1, var) > 0, ]
    mean.ab <- c(apply(rpm.a, 1, mean), apply(rpm.b, 1, mean))
    var.ab  <- c(apply(rpm.a, 1, var), apply(rpm.b, 1, var))
    dispersion <- (var.ab - mean.ab) / (mean.ab * mean.ab)
    population <- data.frame(mean = mean.ab, disp = dispersion)
    population <- population[population$disp > 0, ]
    resampling.vector <- sample(1:nrow(population), Ngene, replace = TRUE)
    population <- population[resampling.vector, ]


    ## simulating data
    count <- matrix(0, ncol = ncol(group), nrow = Ngene)
    count <- apply(fc.matrix, 2, function(x, pp = population) {
                      rnbinom(n = Ngene,
                              mu = x * pp$mean,
                              size = 1 / pp$disp)
                   }, population)
    if (!is.null(replicates)) {
        colnames(count) <- paste("G", rep(1:length(replicates),
                                          times = replicates),
                                 "_rep", sequence(replicates), sep = "")
    } else {
        repnm <- apply(group, 1, function(i){paste(i, collapse="")})
        colnm <- repnm
        tb <- table(repnm)
        tbm <- tb + 1
        for (i in 1:length(repnm)) {
          colnm[i] <- paste(repnm[i], paste("rep",
               tbm[repnm[i]] - tb[repnm[i]], sep = ""), sep = "_")
          tb[repnm[i]] <- tb[repnm[i]] - 1
        }
        colnames(count) <- colnm
    }
    rownames(count) <- paste("gene", 1:nrow(count), sep = "_")
    names(trueDEG) <- rownames(population) <- rownames(fc.matrix) <- rownames(count)
    colnames(fc.matrix) <- colnames(count)

    ## TCC constructor
    tcc <- new("TCC", count,
               if(is.null(replicates)) group
               else rep(1:length(replicates), times = replicates))
    tcc$simulation$trueDEG <- trueDEG
    tcc$simulation$DEG.foldchange <- fc.matrix
    tcc$simulation$PDEG <- PDEG * DEG.assign
    tcc$simulation$params <- population
    tcc$private$simulation.rep <- 
               if(is.null(replicates)) group 
               else rep(1:length(replicates), times = replicates)
    tcc$private$simulation <- TRUE
    tcc$private$estimated <- FALSE
    return(tcc)
}

calcAUCValue <- function(tcc, t = 1) {
    if (t < 0 || 1 < t) stop("\nTCC::ERROR: 't' is limited in (0, 1)")

    if (is.null(tcc$simulation$trueDE) || length(tcc$simulation$trueDE) == 0)
        stop("\nTCC::ERROR: No true positive annotations about differential expression genes.\n ")
    if (is.null(tcc$stat$rank) || length(tcc$stat$rank) == 0)
        stop("\nTCC::ERROR: There are no rank informations in TCC tcc. It need run TCC.estimateDE().\n")
      return(pAUC(rocdemo.sca(truth = as.numeric(tcc$simulation$trueDE != 0), 
                             data = - tcc$stat$rank), t0 = t))
}

plotFCPseudocolor <- function(tcc, main = "",
                              xlab = "samples", ylab = "genes") {
    if (is.null(tcc$simulation$trueDEG) || length(tcc$simulation$trueDEG) == 0)
      message("\nTCC::ERROR: There is no annotations about simulation data.\n")
    d <- tcc$simulation$DEG.foldchange
    layout(matrix(data = c(1, 2), nrow = 1, ncol = 2), 
           widths = c(4, 1), heights=c(1, 1))
    maxlevel <- ceiling(max(tcc$simulation$DEG.foldchange))
    minlevel <- ceiling(1 / min(tcc$simulation$DEG.foldchange))
    d[d < 1] <- - 1 / d[d < 1] + 2
    if (min(d) >= 1) {
        colorRamp <- c(
            "#FFFFFFFF",
            cm.colors((maxlevel - 1) * 32)[((maxlevel - 1) 
                                       * 16):((maxlevel - 1) * 32 - 1)]
        )
    } else if (max(d) <= 1) {
        colorRamp <- c(
            cm.colors((minlevel - 1) * 32)[2:((minlevel - 1) * 16)],
            "#FFFFFFFF"
        )
    } else {
        colorRamp <- c(
            cm.colors((minlevel - 1) * 32)[2:((minlevel - 1) * 16)],
            "#FFFFFFFF",
            cm.colors((maxlevel - 1) * 32)[((maxlevel - 1) 
                                       * 16):((maxlevel - 1) * 32 - 1)]
        )
    }

    colorLevels <- seq(2 - minlevel, maxlevel, length = length(colorRamp))
    par(mar = c(3 + ncol(tcc$group) * 0.6, 4.5, 2.5, 2))
    image(1:ncol(d), 1:nrow(d), t(d[rev(1:nrow(d)), ]),
          col = colorRamp, ylab = ylab, xlab = "", main = main, axes = FALSE, 
          zlim = range(2 - minlevel, maxlevel))
    title(xlab = xlab, line = 1 + ncol(tcc$group))
    for (i in 1:ncol(tcc$group)) {
        axis(1, at = 1:nrow(tcc$group), labels = tcc$group[, i],
             cex.axis = 0.8,
             line = i * ifelse(i == 1, 1, 0.6) - ifelse(i == 1, 1, 0.6) ,
             tick = as.logical(i == 1),
             lty = as.numeric(i != 0))
        mtext(colnames(tcc$group)[i], side = 1, at = -0, 
              cex = 0.8, adj = 1, 
              line = i * ifelse(i == 1, 1, 0.6) - ifelse(i == 1, 1, 0.6) + 1)
    }
    ycoor <- unique(c(0, cumsum(nrow(tcc$count) * tcc$simulation$PDEG),
           nrow(tcc$count) - 0.5))
    yaxis <- unique(sprintf("%d", c(0, cumsum(nrow(tcc$count) * tcc$simulation$PDEG),
           nrow(tcc$count))))
    ycoor[ycoor == 0] <- 1
    yaxis[yaxis == "0"] <- "1"
    axis(2, at = nrow(tcc$count) - ycoor, labels = yaxis, 
         cex.axis = 0.7, las = 1)
    box()
    par(mar = c(3 + ncol(tcc$group) * 0.6, 2.5, 2.5, 2))
    image(1, 0:length(colorRamp),
          matrix(colorLevels, ncol = length(colorRamp), nrow = 1),
          col = colorRamp, xlab = "", ylab = "",
          xaxt = "n", yaxt="n")
    lb <- seq(from = - minlevel + 2, to = maxlevel, by = 1)
    lc <- lb
    lc[lc < 1] <- 1 / (2 - lc[lc < 1])
    axis(2,
         at = seq(from = 0, to = length(colorRamp),
                  by = length(colorRamp) / (length(lb) - 1)),
         labels =  c(rev(paste("1/", 1:minlevel, sep = "")[-1]),
                     sprintf("%d", 1:maxlevel)),
         cex.axis = 0.8)
    box()
    layout(1)
}



