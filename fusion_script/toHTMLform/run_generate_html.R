
library(knitr)

args<-commandArgs(TRUE)
if (length(args) != 3){
    args <- c("--help")
}

## Help section
if("--help" %in% args) {
cat("
    Required three arguments:
    arg1: Path to function file
    arg2: Path to visualization file
    arg3: Path to your output file
        \n")
    q(save="no")
}


funBin<-args[1]
visualFile=args[2]
outputFile=args[3]

knit2html(funBin,output = outputFile)
