test_quantro <- function() {
    library(minfi)
    data(flowSorted)
  
    ## Test 1: Using matrix object and groupFactor
    p <- getBeta(flowSorted, offset = 100)
    pd <- pData(flowSorted)
    qtest <- quantro(object = p, groupFactor = pd$CellType)
    checkEqualsNumeric(quantroStat(qtest), 8.807, tolerance=1.0e-4)

    ## Test 2: Using MethylSet object and groupFactor
    qtest <- quantro(object = flowSorted, groupFactor = pd$CellType)
    checkEqualsNumeric(quantroStat(qtest), 8.807, tolerance=1.0e-4)

    ## Test 3: Using MethylSet object, but with parallelization feature
    if(require(doParallel)){
        registerDoParallel(cores = 4)
        
        pd <- pData(flowSorted)
        qtest <- quantro(flowSorted, pd$CellType, B = 100)
        checkEqualsNumeric(quantroStat(qtest), 8.807, tolerance=1.0e-4)
    }
}

