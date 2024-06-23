#!/usr/bin/env Rscript

# A script that predicts proteasome cleavage predictions for proteins to obtain a list of antigen candidates

# usage: $ /path/to/preoteasome_cleavage.R -nt 0.5 input_file.fasta output_file.tsv

# Load needed packages and variables needed for Netchop
library('argparse')
NETCHOP_PATH <- '/icgc/dkfzlsdf/analysis/G200/immuno/tools/netchop-3.1/netchop'
Sys.setenv('NETCHOP' = NETCHOP_PATH)

# Command-line arguments parsing
parser <- ArgumentParser() # Parser object
parser$add_argument('-nt', '--netchop-threshold', type = 'double', default = '0.5', 
                    help = 'The minimum NetChop threshold to keep a predicted cleavage site. ')
parser$add_argument('input_file', type = 'character', default = '', 
                    help = 'The path/name of your .fasta input file, containing one or more proteins.', nargs = 1)
parser$add_argument('output_file', type = 'character', default = '', 
                    help = 'The path/name of the output file (predicted antigenes table)', nargs = 1)

args <- parser$parse_args() # Parse arguments

# Check if the input file exists
inputFile <- args$input_file
if(file.access(inputFile) == -1) {
    stop(sprintf('Input file (%s) for Netchop does not exist.',inputFile))
}

# A function that exploits a data.frame to give results for one protein
getResults <- function(resultsDF, lenAg, lenFlank) {
    # Get cleavage sites position
    cleavageSites <- which(resultsDF[, 2]=='S')
    # Get protein sequence
    peptideSequence <- paste0(resultsDF[, 1],collapse='')
    # Filter only compatible sites (in terms of length)
    cleavageSites <- cleavageSites[(cleavageSites >= lenAg+lenFlank) & (cleavageSites <= nchar(peptideSequence)-lenFlank)]
    # Get antigenes, insert sequence, associated netchop score for each predicted cleavage site.
    insertGenes <- sapply(cleavageSites, function(x) substr(peptideSequence, x-(lenAg+lenFlank-1), x+lenFlank))
    antiGenes <- sapply(insertGenes, function(x) substr(x, lenFlank+1, nchar(x)-lenFlank), USE.NAMES=FALSE)
    ncScores <- as.numeric(resultsDF[cleavageSites, 3])
    # Final data.frame
    return(data.frame('antigene' = antiGenes, 'insertGene' = insertGenes, 'netchopScore' = ncScores))
}

# A function that gets all results (one or more proteins) and save them as output file. :)
fullResults <- function(resultsDFList, lenAg, lenFlank, outputFile){
    write(NULL, file=outputFile)
    count <- 1
    for(i in resultsDFList){
        if(count == 1){
            write.table(getResults(i, lenAg, lenFlank), file = outputFile, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', append = TRUE)
        }
        else {
            write.table(getResults(i, lenAg, lenFlank), file = outputFile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t', append = TRUE)
        }
    count<-count+1
    }
    sprintf('Results saved in %s',outputFile)
}

# Call Netchop
print('Calling Netchop')
#tmpResultFile <- tempfile(pattern = 'netchop_result', fileext = '.out')  # Zhiqin

tmpResultFile <- paste0("_temp_netchop_result_.out")

netchopCommand <- paste(
#    paste0(NETCHOP_PATH,'/bin/netChop'),  # Netchop binary, Zhiqin
    NETCHOP_PATH,
    '-t', args$netchop_threshold,         # Netchop minimum threshold
#    '-tdir', paste0(NETCHOP_PATH,'/tmp'), # Netchop tmp direcotry, Zhiqin
    '-tdir .',	
    args$input_file,                      # Input
    '>', tmpResultFile)
system(netchopCommand) # Execute command

print('Netchop run completed. Parsing results...')
# Parse Netchop results
resultLines <- readLines(tmpResultFile)
# Starting and end position of each peptide of the result table
startPos <- which(grepl('^ pos', resultLines))+2
endPos <- which(grepl('^Number', resultLines))-3

# Get results (list of data.frame, one per protein)
resultList <- list()
for(i in 1:length(startPos)){
    currentProtein <- resultLines[startPos[i]:endPos[i]]
    # To df ! 
    dfProtein <- data.frame(do.call(rbind, (strsplit(currentProtein, ' +'))), stringsAsFactors=FALSE)
    dfProtein[, c(1, 2)] <- list(NULL) # Remove the two first columns (useless)
    names(dfProtein) <- c('AA', 'cleavage', 'score', 'id')
    resultList[[i]] <- dfProtein # Save in result list
}

# Get and save results
fullResults(resultList, 9, 2, args$output_file)
print('Proteasome cleavage prediction done.')
