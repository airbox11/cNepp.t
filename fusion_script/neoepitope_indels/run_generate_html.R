
library(knitr)
args = commandArgs(trailingOnly=TRUE)
funBin<-args[1]
visualFile=args[2]
outputFile=args[3]

knit2html(funBin,output = outputFile)
