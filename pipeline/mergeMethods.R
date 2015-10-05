mergeExpression <- function(){
  #fileCheck <- checkConsist()
  fileCheck <- c()
  if (length(fileCheck) == 0){
    myFiles <- list.files()
    curFile <- readLines(myFiles[1])
    myRows <- length(curFile) - 1
    myCols <- length(myFiles) + 1
    newFile <- matrix(nrow = myRows, ncol = myCols)
    geneList <- c("GeneName")
    for (i in 3:length(curFile)){
      curRow <- curFile[i]
      curRow = strsplit(curRow,"\\s+")
      nextGene = curRow[[1]][1]
      geneList <- c(geneList,nextGene)
    }
    newFile[,1] = geneList
    vacControl = 2
    vacTreated = myCols
    for (i in 1:length(myFiles)){
      print(100.0 * (i + 0.0) / (length(myFiles) + 0.0))
      curFile <- readLines(myFiles[i])
      colName <- strsplit(curFile[1],"\\s+")[[1]][3]
      sampleCode = strsplit(colName,"-")[[1]][4]
      sampleCode = substr(sampleCode,1,2)
      sampleCode = as.numeric(sampleCode)
      nextCol = c(colName)
      for (i in 3:length(curFile)){
        curRow <- curFile[i]
        nextCol = c(nextCol,strsplit(curFile[i],"\\s+")[[1]][2])
      }
      if (sampleCode > 9){
        newFile[,vacControl] = nextCol
        vacControl = vacControl + 1
      }
      else{
        newFile[,vacTreated] = nextCol
        vacTreated = vacTreated - 1
      }
    }
    return(newFile)
  }
}

mergeMethylation <- function(){
  #fileCheck <- checkConsist()
  fileCheck <- c()
  if (length(fileCheck) == 0){
    myFiles <- list.files()
    curFile <- readLines(myFiles[1])
    myRows <- length(curFile) - 1
    myCols <- length(myFiles) + 4
    newFile <- matrix(nrow = myRows, ncol = myCols)
    chrList <- c("Chromosome")
    geneList <- c("Gene")
    locList <- c("Genomic_Coordinate")
    refList <- c("Reference")
    for (i in 3:length(curFile)){
      curRow <- curFile[i]
      curRow = strsplit(curRow,"\\s+")
      if (length(curRow[[1]]) == 5){
        nextRef = curRow[[1]][1]
        nextGene = curRow[[1]][3]
        nextChr = curRow[[1]][4]
        nextLoc = curRow[[1]][5]
      }
      else{
        nextRef = curRow[[1]][1]
        nextGene = "NA"
        nextChr = curRow[[1]][3]
        nextLoc = curRow[[1]][4] 
      }
      chrList <- c(chrList,nextChr)
      geneList <- c(geneList,nextGene)
      locList <- c(locList,nextLoc)
      refList <- c(refList,nextRef)
    }
    newFile[,1] = chrList
    newFile[,2] = geneList
    newFile[,3] = refList
    newFile[,4] = locList
    vacControl = 5
    vacTreated = myCols
    for (i in 1:length(myFiles)){
      print(100.0 * (i + 0.0) / (length(myFiles) + 0.0))
      curFile <- readLines(myFiles[i])
      colName <- strsplit(curFile[1],"\\s+")[[1]][3]
      sampleCode = strsplit(colName,"-")[[1]][4]
      sampleCode = substr(sampleCode,1,2)
      sampleCode = as.numeric(sampleCode)
      nextCol = c(colName)
      for (i in 3:length(curFile)){
        curRow <- curFile[i]
        nextCol = c(nextCol,strsplit(curFile[i],"\\s+")[[1]][2])
      }
      if (sampleCode > 9){
        newFile[,vacControl] = nextCol
        vacControl = vacControl + 1
      }
      else{
        newFile[,vacTreated] = nextCol
        vacTreated = vacTreated - 1
      }
    }
    return(newFile)
  }
}


checkConsist <- function(){
  myFiles <- list.files()
  badFiles <- c()
  if (length(myFiles) > 0){
    curFile <- readLines(myFiles[1])
    tempLength <- length(curFile)
    for (i in 2:length(myFiles)){
      curFile <- readLines(myFiles[i])
      if (length(curFile) != tempLength){
        badFiles = c(badFiles, myFiles[i])
      }
    }
  }
  return(badFiles)
}
