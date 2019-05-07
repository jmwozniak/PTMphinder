#####################################################################################
#                                                                                   #
#         Copyright (C) 2019 Jacob M. Wozniak                                       #
#                                                                                   #
#         This file is part of the PTMphinder R package.                            #
#                                                                                   #
#         PTMphinder is free software: you can redistribute it and/or modify        #
#         it under the terms of the GNU General Public License as published by      #
#         the Free Software Foundation, either version 3 of the License, or         #
#         any later version.                                                        #
#                                                                                   #
#         PTMphinder is distributed in the hope that it will be useful,             #
#         but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#         GNU General Public License for more details.                              #
#                                                                                   #
#         You should have received a copy of the GNU General Public License         #
#         along with PTMphinder. If not, see <https://www.gnu.org/licenses/>.       #
#                                                                                   #
#####################################################################################

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







