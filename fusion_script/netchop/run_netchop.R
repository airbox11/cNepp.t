## to run netchop 

NETCHOP_PATH <- '/icgc/dkfzlsdf/analysis/G200/immuno/tools/netchop-3.1/netchop'
#/desktop-home/huangz/immunoinfo/tools/netchop-3.1/netchop 

Sys.setenv('NETCHOP' = NETCHOP_PATH)

funNetchop <- function(nameGene, peptide){
  ## generate a fasta input file
  tmpFastaFile <- paste0("__temp_fasta_for_netchop_.fa")
  pepFasta <- rbind(paste(">",nameGene,sep=""),peptide)
  write.table(pepFasta, file = tmpFastaFile,quote = F, col.names = F, row.names = F)
  # netchop tmp file
  tmpNetchopResultFile <- paste0("__temp_netchop_result_.out")
  # run netchop 
  netchopCommand <- paste(NETCHOP_PATH,'-tdir .', tmpFastaFile, '>', tmpNetchopResultFile)
  system(netchopCommand)
  # read netchop result file
  resultLines <- readLines(tmpNetchopResultFile)
  # identify the start and end position of each input peptide from defaul netchop output
  # multiple fasta input considered here
  pepStartPos <- which(grepl('^ pos', resultLines))+2
  pepEndPos <- which(grepl('^Number', resultLines))-3
  
  # extract the netchop score information and reform them into a list
  pepList <- list()
  for(i in 1:length(pepStartPos)){
    currentPep <- resultLines[pepStartPos[i]:pepEndPos[i]]
    dfPep <- data.frame(do.call(rbind, (strsplit(currentPep, ' +'))), stringsAsFactors=FALSE)
    dfPep[, c(1)] <- list(NULL) # delete first empty column
    names(dfPep) <- c('pos','AA', 'cleavage', 'score', 'Ident')
    pepList[[i]] <- dfPep
    # transform netchop score to 2 rows.
    # here consider only one fasta input
    netchopScore <- t(cbind(pepList[[i]][,4],pepList[[i]][,3]))
    colnames(netchopScore) <- t(pepList[[i]][,1])
    netchopScore[1,] <- round(as.numeric(netchopScore[1,]),digits=2)
  }

  # remove tmp files
  rmCommand1 <- paste("rm ",tmpFastaFile)
  system(rmCommand1)
  rmCommand2 <- paste("rm ",tmpNetchopResultFile)
  system(rmCommand2)
  # return results
  return(netchopScore)
}

