
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  print('parameters needed for converting to xlsx')
} else{
  workDir <- args[1]
  mhc.input <- args[2]
}

## test ==== ====
test <- TRUE
test <- FALSE
if (test){
  workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MO79'
  mhc.input <- 'MHCI_epitopes_RNAseq_netMHCpan4_1.tab_renameCol'

  workDir <- '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/K9665'
  mhc.input <- 'MHCI_epitopes_TCGA-COAD_netMHCpan4_1.tab_renameCol'
}


main_run <- function() {
  dir1 <- file.path(workDir, '8_chose_neoepitode')
  setwd(dir1)
  ## MHCI table ==============
  
  df.mhcI <- read.table(mhc.input, header=TRUE, stringsAsFactors = FALSE, sep = '\t')
  df.mhcI$LossAllele <- tolower(gsub("HLA-(.)\\*(..):(..)", "hla_\\1_\\2_\\3", df.mhcI$MHC_allele))
  
  
  ## loh table ==============
  file.loh <- '../5_LOHHLA/example-out/example.10.DNA.HLAlossPrediction_CI.xls'
  
  if(!file.exists(file.loh)){
    print(file.path(workDir,'5_LOHHLA/example-out'))
    print('=========== =========== no LOH result exist. ')
    return(1)

  }
  df.loh <- read.table(file.loh, header=TRUE, stringsAsFactors = FALSE, sep = '\t')
  if(ncol(df.loh) == 1){
    df.loh <- t(df.loh)
    colnames(df.loh) <- c('region','HLA_A_type1','HLA_A_type2','HLAtype1Log2MedianCoverage','HLAtype2Log2MedianCoverage','HLAtype1Log2MedianCoverageAtSites','HLAtype2Log2MedianCoverageAtSites','HLA_type1copyNum_withoutBAF','HLA_type1copyNum_withoutBAF_lower','HLA_type1copyNum_withoutBAF_upper','HLA_type1copyNum_withBAF','HLA_type1copyNum_withBAF_lower','HLA_type1copyNum_withBAF_upper','HLA_type2copyNum_withoutBAF','HLA_type2copyNum_withoutBAF_lower','HLA_type2copyNum_withoutBAF_upper','HLA_type2copyNum_withBAF','HLA_type2copyNum_withBAF_lower','HLA_type2copyNum_withBAF_upper','HLA_type1copyNum_withoutBAFBin','HLA_type1copyNum_withoutBAFBin_lower','HLA_type1copyNum_withoutBAFBin_upper','HLA_type1copyNum_withBAFBin','HLA_type1copyNum_withBAFBin_lower','HLA_type1copyNum_withBAFBin_upper','HLA_type2copyNum_withoutBAFBin','HLA_type2copyNum_withoutBAFBin_lower','HLA_type2copyNum_withoutBAFBin_upper','HLA_type2copyNum_withBAFBin','HLA_type2copyNum_withBAFBin_lower','HLA_type2copyNum_withBAFBin_upper','PVal','UnPairedPval','PVal_unique','UnPairedPval_unique','LossAllele','KeptAllele','numMisMatchSitesCov','propSupportiveSites')
    rownames(df.loh) = seq(length=nrow(df.loh))
    write.table(df.loh, file=file.loh, col.names=TRUE, row.names=FALSE,quote = FALSE, sep = '\t')
    df.loh.lite <- as.data.frame(t(df.loh[,c("LossAllele", "PVal_unique")]))
  } else {
    df.loh.lite <- df.loh[,c("LossAllele", "PVal_unique")]
  }
  
  ## merge tables ==============
  df.mhcI.loh <- merge(df.mhcI, df.loh.lite, by.x='LossAllele', by.y='LossAllele', all.x=TRUE)
  df.mhcI.loh[!is.na(df.mhcI.loh$PVal_unique),]$LossAllele <- 'Yes'
  df.mhcI.loh[is.na(df.mhcI.loh$PVal_unique),]$LossAllele <- 'No'
  
  df.mhcI.loh2 <- df.mhcI.loh[,c(2,3,ncol(df.mhcI.loh),1,4:(ncol(df.mhcI.loh)-1))]
  mhc.output <- paste0(mhc.input, '_loh')
  write.table(df.mhcI.loh2, file = mhc.output, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
}
main_run()
