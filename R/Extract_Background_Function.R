

extractBackground <- function(seqs, char, width){
  library(stringr)
  fullProt <- NULL
  flankSeq <- NULL
  flankList <- rep("",sum(str_count(seqs, char)))
  halfWidth <- floor(width/2)
  len<-length(seqs)
  count <- 0

  for(y in 1:len) {

    fullProt <- seqs[y]
    tempLoc <- gregexpr(pattern = char, fullProt)[[1]]

    if(tempLoc[[1]]!=-1){

      for(i in 1:length(tempLoc)){

        count <- count + 1
        flankSeq <- substr(fullProt, tempLoc[i]-halfWidth, tempLoc[i]+halfWidth)
        flankList[count] <- flankSeq

      }
    }

    fullProt <- NULL
    tempLoc <- NULL
  }

  output <- unique(flankList[which(lapply(flankList, nchar)==width)])

  return(output)
}







