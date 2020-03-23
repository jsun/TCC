clusterSample <- function (data, dist.method = "spearman",
                           hclust.method = "average", unique.pattern = TRUE) {
##  This function performs hierarchical clustering for samples (tissues or 
##  columns) from expression data.
    if (class(data)[1] == "TCC") {
        data <- data$count
    }
    data <- data[rowSums(data) > 0, ]
    if (unique.pattern) {
        data <- unique(data)
    }
    d <- as.dist(1 - cor(data, method = dist.method))
    h <- hclust(d, method = hclust.method)
    return (h)
}

