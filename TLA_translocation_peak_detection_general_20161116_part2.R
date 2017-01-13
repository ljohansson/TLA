#This script will predict peak locations for TLA translocation data using full 10kb bin read count data, including all positions having a coverage of 0.
#To do this first noise is filtered in five steps:
#1. For all bins if the total read count of the two bins before and the two bins after the processed bin is 0 the read count of the bin is set to zero.
#   This will remove isolated peaks, caused by for instance repetitive sequences. This filter is called the TwoBesideFilter.
#2. The bins added up into 100 kb bins starting from the first bin in each chromosome. If there are less then ten bins left at the end of a chromosome,
#   the last bin will be smaller.
#3. In the 100 kb bins only bins having a read count over a set threshold (currently 1 / 30000 of the read count of the MAPQ 3 filtered bam file) are kept. In bins having a lower read count, the read count is set to 0.
#4. For all corrected bins out of the two bins before and the two bins after the processed bin, the total read count of the three bins having the lowest read count
#   is calculated. If the result is 0 the read count of the processed bin is set to 0. This filter is called the LowestThreeTwoBesideFilter.
#5. Only the bins having a read count above a set threshold (currently 1 / 6000 of the read count of the MAPQ 3 filtered bam file) are kept.In bins having a lower read count, the read count is set to 0.
[umcg-ljohansson@calculon TLA_scripts]$ head -20 TLA_translocation_peak_detection_general_20161116.R
#This script will predict peak locations for TLA translocation data using full 10kb bin read count data, including all positions having a coverage of 0.
#To do this first noise is filtered in five steps:
#1. For all bins if the total read count of the two bins before and the two bins after the processed bin is 0 the read count of the bin is set to zero.
#   This will remove isolated peaks, caused by for instance repetitive sequences. This filter is called the TwoBesideFilter.
#2. The bins added up into 100 kb bins starting from the first bin in each chromosome. If there are less then ten bins left at the end of a chromosome,
#   the last bin will be smaller.
#3. In the 100 kb bins only bins having a read count over a set threshold (currently 1 / 30000 of the read count of the MAPQ 3 filtered bam file) are kept. In bins having a lower read count, the read count is set to 0.
#4. For all corrected bins out of the two bins before and the two bins after the processed bin, the total read count of the three bins having the lowest read count
#   is calculated. If the result is 0 the read count of the processed bin is set to 0. This filter is called the LowestThreeTwoBesideFilter.
#5. Only the bins having a read count above a set threshold (currently 1 / 6000 of the read count of the MAPQ 3 filtered bam file) are kept.In bins having a lower read count, the read count is set to 0.
#6. On the resulting data again the TwoBesideFilter (as in step 1) is applied.

#Start file using
#Rscript filename filepath


#Functions

############################FUNCTION MERGE 10 BINS###################################
Merge_10_bins <- function(BinsList){
  message ("Start Merge 10 bins")
  #Put list in dataframe for easy processing
  BinsList_df <- na.fail(as.data.frame(BinsList))
  
  #Put rows into separate variables
  Chrom <- (BinsList_df[,1, drop=FALSE])
  Pos <- (BinsList_df[,2, drop=FALSE])
  RC <- (BinsList_df[,3, drop=FALSE])
  
  #Setting ProcessedBinChrom and BinBeforeChrom for first cycle and set BinCounter for first loop
  ProcessedBinChrom <- 0
  BinBeforeChrom <- 0
  BinCounter <- 1
  TableLine <- 0
  
  #Prepare list for output
  MergedBinsList <- list()
  
  #Setting start and endline
  Startpos <- 1
  Endpos <- nrow(BinsList_df)
  
  #Looping through all rows to merge per 10
  for (Line in Startpos:Endpos){                          #Gets information for all lines 
  
    BinBeforeChrom <- ProcessedBinChrom 
    
    #Transfers non-numerical chromosome numbers to numerical ones X -> 23, Y -> 24, MT and GL -> 25
    if ( suppressWarnings(is.na(as.numeric(as.matrix(Chrom)[Line,]))) ){
      if (as.matrix(Chrom)[Line,] == "X"){
        ProcessedBinChrom <- 23
      } else if (as.matrix(Chrom)[Line,] == "Y"){
        ProcessedBinChrom <- 24
      } else {
        ProcessedBinChrom <- 25
      }
    } else 
      ProcessedBinChrom <- as.numeric(as.matrix(Chrom)[Line,]) #Gets the chromosome of the bin of the current line read in the
    ProcessedBinPos <- as.numeric(as.matrix(Pos)[Line,])     #idem for Chromosomal Position
    ProcessedBinRC <- as.numeric(as.matrix(RC)[Line,])       #idem for Read count
    
    #Catch last bins of chromosomes (Lags a loop behind)
    if (BinBeforeChrom != 0 && BinBeforeChrom != ProcessedBinChrom)  {
      MergedBin <- (c(as.integer(format(MergedBinChrom, scientific = FALSE)), as.integer(format(MergedBinPos, scientific = FALSE)), as.integer(format(MergedBinRC, scientific = FALSE))))
      MergedBinsList[[as.numeric(TableLine)]] <- MergedBin  #Fill position
    }
    
    
    #Set first bin of each chromosome
    if (ProcessedBinChrom != BinBeforeChrom){
      BinCounter <- 1
      message("Processing chromosome ", ProcessedBinChrom)
    }
    
    #Add up bins
    if (BinCounter == 1){
      TableLine <- TableLine + 1
      MergedBinChrom <- ProcessedBinChrom
      MergedBinPos <- ProcessedBinPos
      MergedBinRC <- ProcessedBinRC
    }
    if (BinCounter > 1 ){
      MergedBinRC <- MergedBinRC + ProcessedBinRC
    }  
    
    if (BinCounter == 10){
      MergedBin <- (c(as.integer(format(MergedBinChrom, scientific = FALSE)), as.integer(format(MergedBinPos, scientific = FALSE)), as.integer(format(MergedBinRC, scientific = FALSE))))
      MergedBinsList[[as.numeric(TableLine)]] <- MergedBin  #Fill position
      BinCounter <- 0
    }
    
    #Catch last bins of file 
    if (Line == Endpos)  {
      MergedBin <- (c(as.integer(format(MergedBinChrom, scientific = FALSE)), as.integer(format(MergedBinPos, scientific = FALSE)), as.integer(format(MergedBinRC, scientific = FALSE))))
      MergedBinsList[[as.numeric(TableLine)]] <- MergedBin  #Fill position
    }
    BinCounter <- BinCounter + 1
    
  }
  MergedBinslist_df <- data.frame(t(as.data.frame(MergedBinsList)))
  rownames(MergedBinslist_df) <- 1:nrow(MergedBinslist_df)
  return(MergedBinslist_df) #Note sex chromosomes and other genomic regions are kept in format 23, 24, 25
}





###################################################################################################
##############     MAIN    ###################Main (Make sure analysis file is read with no header)
###################################################################################################
library(tools) #To enable getting basename and extensions
#NOTE THERE IS A DEPENDENCY ON RSAMTOOLS IN THE COUNT ALIGNED READS FUNCTION

args <- commandArgs(TRUE)
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

sample_name <- args[1]
inputpath <- args[2]
tmpfolder <- paste(inputpath,"/tmp/",sep="")
dir.create(tmpfolder)

message ("Analyzing sample ", sample_name)
#AnalysisFile <- read.csv(file = sample_name, sep = "\t", header = FALSE)
#BamFile_name <- paste(substr(sample_name, 1, regexpr('ROI', sample_name ) - 2 ),".bam", sep="") #For logic see https://www.r-bloggers.com/basic-text-string-functions-in-r/

#message ("Step 0 counting aligned read number total bam file")
#Aligned_reads <- CountAlignedReadsBam(BamFile_name)
#message ("The number of aligned reads are ", Aligned_reads)

message ("Step 1 Two Besides Filter has been performed in part 1")
First_correction_TwoBesidesFilter <- read.csv(file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt.txt",sep=""), sep = "\t", header = FALSE)
#First_correction_TwoBesidesFilter  <- TwoBesideFilter(AnalysisFile)
#write.table(x = First_correction_TwoBesidesFilter, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
#rm(AnalysisFile)

message ("Step 2 Merge to 100kb bins")
Mergedto100kbBins <- Merge_10_bins(First_correction_TwoBesidesFilter)
write.table(x = Mergedto100kbBins, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
message ("Done")
