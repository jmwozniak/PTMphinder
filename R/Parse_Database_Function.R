parseDB <- function(ref, db_source, filt){
  #Count number of sequences in .fasta database
  numSeqs <- sum(substr(ref, 0, 1) == ">")

  #Paste together protein sequences from .fasta database
  reftab <- rep("",2*numSeqs)
  temp <- NULL
  count <- 1
  for(i in 1:length(ref)){
    if(substr(ref[i], 0, 1) == ">"){
      reftab[count] <- ref[i]
      count <- count + 1
    }
    else{
      temp <- c(temp, ref[i])
      temp <- paste(temp, sep="", collapse="")
    }

    if((substr(ref[i+1], 0, 1) == ">")||(i == length(ref)))
    {
      reftab[count] <- temp
      count <- count + 1
      temp <- NULL
    }
  }

  #Organize reference data table into two columns (description, sequence)
  reftab <- matrix(reftab,ncol=2,byrow=TRUE)

  #Filter out reverse and contaminant sequences if the filter setting is set to TRUE
  if(filt==TRUE){
    reftab2 <- subset(reftab, (tolower(substr(reftab[,1], 2, 12)) != "contaminant"))
    reftab <- subset(reftab2, (tolower(substr(reftab2[,1], 2, 8)) != "reverse"))
  }

  #Extract accession IDs from .fasta descriptions
  accID <- rep("",numSeqs)
  if(db_source=="UP"){

    tempIDs <- gsub("|", "]", reftab[,1], fixed = TRUE)
    left <- gregexpr(']',tempIDs)[[1]][1]
    right <- gregexpr(']',tempIDs)[[1]][2]
    accID <- substr(tempIDs, left+1, right-1)

  } else if(db_source=="RS") {

    tempIDs <- reftab[,1]
    right <- gregexpr(' ',tempIDs)[[1]][1]
    accID <- substr(tempIDs, 2, right-1)

  } else {

    tempIDs <- reftab[,1]
    accID <- substr(tempIDs, 2, nchar(tempIDs))

  }

  #Replace column 1 with the accessionIDs
  reftab[,1] <- accID
  return(reftab)
}

