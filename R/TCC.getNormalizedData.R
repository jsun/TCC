
TCC$methods(getNormalizedData = function () {
    effective.libsizes <- colSums(count) * norm.factors
    return (sweep(count, 2, 
                  mean(effective.libsizes) / effective.libsizes, "*"))
})

