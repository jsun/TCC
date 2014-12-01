test_new <- function() {
    data(hypoData)
    groupFactor <- factor(c(1, 1, 1, 2, 2, 2))
    groupNum <- c(1, 1, 1, 2, 2, 2)  
    groupStr <- c("G1", "G1", "G1", "G2", "G2", "G2")  

    matrixObj <- as.matrix(hypoData)
    dataframeObj <- as.data.frame(hypoData)

    tccMatrixObj <- new("TCC", matrixObj, groupNum)
    tccDataframeObj <- new("TCC", dataframeObj, groupNum)
    checkEquals(tccMatrixObj, tccDataframeObj)

    tccMatrixObj <- new("TCC", matrixObj, groupStr)
    tccDataframeObj <- new("TCC", dataframeObj, groupStr)
    checkEquals(tccMatrixObj, tccDataframeObj)


    tccMatrixObj <- new("TCC", matrixObj, groupFactor)
    tccDataframeObj <- new("TCC", dataframeObj, groupFactor)
    checkEquals(tccMatrixObj, tccDataframeObj)
}
