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

phindPTMs <- function(data, reftab){
  #Read input data into vectors
  protIDs <- data$Protein_ID
  pepSeqs <- toupper(data$Peptide_Seq)
  ptmLocs <- gsub(";", ".", data$PTM_Loc, fixed=TRUE)
  ptmScore <- data$PTM_Score
  numPTMs <- str_count(data$PTM_Loc,";") + 1
  totLocs <- data$Total_Sites
  numPeptides <- length(ptmLocs)

  ambiguity <- ifelse(numPTMs == totLocs, "Unambiguous", "Ambiguous")

  #Reformat reference table
  refID <- unlist(reftab[,1][reftab[,1] %in% data$Protein_ID])
  refSeq <- unlist(reftab[,2][reftab[,1] %in% data$Protein_ID])
  reftab2 <- c(as.character(refID),as.character(refSeq))
  reftab2 <- matrix(reftab2,ncol=2,byrow=FALSE)

  #Find pepSeqs in the reference database and get new location of phosphosites in full protein
  newLocList <- rep("",numPeptides)
  flankList <- rep("",numPeptides)
  rows <- nrow(reftab2)
  fullProt <- NULL
  temp2 <- NULL
  temp3 <- NULL
  flanking <- NULL
  for(i in 1:numPeptides){
    for(y in 1:rows){
      if(data$Protein_ID[i]==reftab2[y,1])
      {
        fullProt <- c(fullProt, reftab2[y,2])
        tempLoc <- regexpr(pepSeqs[i], reftab2[y,2])[1]
        pepLocs <- strsplit(ptmLocs[i], ".", fixed = TRUE)

        for(x in 1:length(pepLocs[[1]])){
          temp <- as.numeric(substr(pepLocs[[1]][x], 2, nchar(pepLocs[[1]][x])))
          site <- substr(pepLocs[[1]][x], 1, 1)
          newLoc <- temp + tempLoc - 1
          flanking <- substr(fullProt[[i]], newLoc-7, newLoc+7)
          newSiteLoc <- paste(c(site, newLoc), sep="", collapse = "")
          temp2 <- c(temp2, newSiteLoc)
          temp3 <- c(temp3, flanking)
        }
        temp2 <- paste(temp2, sep="", collapse=".")
        newLocList[i] <- temp2
        temp2 <- NULL
        temp3 <- paste(temp3, sep="", collapse=".")
        flankList[i] <- temp3
        temp3 <- NULL
        break
      }
    }
  }

  #Create table with ProteinID, Unique Identifier, PTM Loc(s) on Peptide, PTM Loc(s) on Protein, PTM scores, flanking sequences, ambiguity of PTM, and full protein sequence
  writelist<- c(as.character(data$Identifier), as.character(protIDs), as.character(ptmLocs), as.character(newLocList), as.character(ptmScore), as.character(flankList),  as.character(ambiguity), as.character(fullProt))
  writetab <- matrix(writelist,ncol=8,byrow=FALSE)

  #Name the columns of the table
  colnames(writetab) <- c("Identifier", "Protein_ID", "Pep_Loc",  "Prot_Loc", "PTM_Score", "Flank_Seq", "Ambiguity", "Prot_Seq")

  return(writetab)
}



