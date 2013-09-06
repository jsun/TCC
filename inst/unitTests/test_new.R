test_new <- function() {
    data(hypoData)
    groupNum <- c(1, 1, 1, 2, 2, 2)  
    groupStr <- c("G1", "G1", "G1", "G2", "G2", "G2")  

    matrixObj <- as.matrix(hypoData)
    dataframeObj <- as.matrix(hypoData)

    tccMatrixObj <- new("TCC", matrixObj, groupNum)
    tccDataframeObj <- new("TCC", dataframeObj, groupNum)
    checkEquals(tccMatrixObj, tccDataframeObj)

    tccMatrixObj <- new("TCC", matrixObj, groupStr)
    tccDataframeObj <- new("TCC", dataframeObj, groupStr)
    checkEquals(tccMatrixObj, tccDataframeObj)
}
