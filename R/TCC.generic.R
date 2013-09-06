plot.TCC <- function(x, FDR = NULL, median.lines = FALSE, floor = 0, 
                     groups = NULL, col.tag = NULL, normalize = TRUE, ...) {
    invisible(x$plotMA(FDR = FDR, median.lines = median.lines, floor = floor,
                       groups = groups, col.tag = col.tag,
                       normalize = normalize, ...))
}

subset.TCC <- function(x, subset, ...){
    if(!is.logical(subset)){
        if(is.numeric(subset)){
            new_v = logical(length(x))
            new_v[subset] <- TRUE
            return(subset(x, new_v))
        }
        if(is.character(subset)){
            new_v = logical(length(x))
            names(new_v) <- x$gene_id
            new_v[subset] <- TRUE
            return(subset(x, new_v))
        }
        message("subset called with unsupported type")
        return(F);
    }
    new_tcc <- new("TCC", as.matrix(x$count[subset, ]), 
                   x$group, x$norm.factors, 
                   as.character(x$gene_id[subset]))
    if (x$private$estimated == TRUE) {
        new_tcc$stat$rank <- x$stat$rank[subset]
        new_tcc$stat$p.value <- x$stat$p.value[subset]
        new_tcc$stat$q.value <- x$stat$q.value[subset]
    }
    if (!is.null(x$estimatedDEG) && length(x$estimatedDEG) > 0){
        new_tcc$estimatedDEG <- x$estimatedDEG[subset]
    }
    if (!is.null(x$simulation)){
        if(length(x$simulation$trueDEG)>0)
            new_tcc$simulation$trueDEG <- x$simulation$trueDEG[subset] 
    if(length(x$simulation$fold.change)>0)
        new_tcc$simulation$fold.change <- x$simulation$fold.change[subset] 
        new_tcc$simulation$PDEG <- x$simulation$PDEG
    }
    new_tcc$private <- x$private
    return(new_tcc)
}

show.TCC <- function(object) {
    ## Counts.
    cat("Count:\n")
    print(head(object$count))
    cat("\n")
    ## Conditions and Annotations.
    df <- data.frame(
        norm.factors = object$norm.factors,
        lib.sizes = object$norm.factors * colSums(object$count)
    )
    rownames(df) <- colnames(object$count)
    df <- cbind(object$group, df)
    cat("Sample:\n")
    print(df)
    cat("\n")
    ## Normalized results.
    if (object$private$normalized) {
        cat("DEGES:\n")
        cat(paste("   Pipeline       : ", 
                  object$DEGES$pipeline, 
                  "\n", sep = ""))
        cat(paste("   Execution time : ", 
                   sprintf("%.1f", object$DEGES$execution.time[3]),
                   " sec\n", sep = ""))
        cat(paste("   Threshold type : ", 
                   object$DEGES$threshold$type, 
                   " < ", 
                   sprintf("%.2f", object$DEGES$threshold$input),
                   "\n",
                   "   Potential PDEG : ", 
                   sprintf("%.2f", sum(object$DEGES$potDEG != 0) /
                                   length(object$DEGES$potDEG)), 
                   "\n\n", sep = ""))
    }
    ## Esimated results.
    if (object$private$estimated) {
        df <- getResult(object)
        cat("Results:\n")
        print(head(df))
        cat("\n")
    }
}



setGeneric(
    name = "calcNormFactors",
    def = function(tcc, ...) tcc)
setMethod(
    f = "calcNormFactors",
    signature(tcc = "DGEList"),
    definition = function(tcc, ...) {
        return(edgeR::calcNormFactors(tcc, ...))
    }
)

setMethod(
    f = "names",
    signature(x = "TCC"),
    definition = function(x) {
        return (c("count", "gene_id", "group", "norm.factors", 
                  "DEGES", "stat", "estimatedDEG", "simulation"))
    }
)

setMethod(
    f = "length",
    signature(x = "TCC"),
    definition = function(x) {
        return (nrow(x$count))
    }
)

setMethod(
    f = "[",
    signature(x = "TCC"),
    definition = function(x, i){
        return(subset(x,i))
    }
)

setMethod(
    f = "subset",
    signature(x = "TCC"),
    definition = subset.TCC
)

setMethod(
    f = "show",
    signature(object = "TCC"),
    definition = show.TCC
)

