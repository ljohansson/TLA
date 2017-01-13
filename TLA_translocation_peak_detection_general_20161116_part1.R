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

############################FUNCTION TWO BESIDE FILTER###################################
TwoBesideFilter <- function(BinsList){
    message("Filtering bins")
    #Put list in dataframe for easy processing
    BinsList_df <- na.fail(as.data.frame(BinsList))
    
    #Put rows into separate variables
    Chrom <- (BinsList_df[,1, drop=FALSE])
    Pos <- (BinsList_df[,2, drop=FALSE])
    RC <- (BinsList_df[,3, drop=FALSE])
    
    #Prepare list for output
    CorrectedBinsList <- list()
    
    #Setting the chromosome checks for the surrounding bins for the loop. 
    TwoBinBeforeChrom <- 0
    BinBeforeChrom <- 0
    ProcessedBinChrom <- 0
    BinAfterChrom <- 0
    TwoBinAfterChrom <- 0
    
    #Setting the read count checks for the surrounding bins for the loop
    TwoBinBeforeRC <- 0
    BinBeforeRC <- 0
    ProcessedBinRC <- 0
    BinAfterRC <- 0
    TwoBinAfterRC <- 0
    
    #Setting the position for the processed bin. The information of the two bins thereafter is needed to fill the processed bin.
    ProcessedBinAfterPos <- 0
    BinAfterPos <- 0
    TwoBinAfterPos <- 0
    
    #Setting start and endline
    Startpos <- 1
    Endpos <- nrow(BinsList_df)
    
    #Options for troubleshooting
    #Endpos <- 200          #NEAR BEGINNING oF CHROMOSOME 1            
    #Startpos <- 100        #START FORFIRST BINS TO BE FILTERED OR NOT
    #Endpos <- 120          #END FORFIRST BINS TO BE FILTERED OR NOT
    #Startpos <- 24910      #NEAR END oF CHROMOSOME 1
    #Endpos <- 24950        #NEAR BEGINNING oF CHROMOSOME 2     
    #Startpos <- 288100     #NEAR END oF CHROMOSOME 22
    #Endpos <- 288120       #NEAR BEGINNING oF CHROMOSOME X  
    #Startpos <- 303630     #NEAR END oF CHROMOSOME X
    #Endpos <- 303650       #NEAR BEGINNING oF CHROMOSOME Y
    #Startpos <- 309570     #NEAR END oF CHROMOSOME Y
    #Endpos <- 309590       #NEAR BEGINNING oF CHROMOSOME FRAGMENTS MT/GL
    
    #Looping through all rows to assess if the coverage of all surrounding bins is 0. Note that the bin that is processed lies two behind the current line read.
    for (Line in Startpos:(Endpos + 2)){
      if (Line > 4){
      TwoBinBeforeChrom <- BinBeforeChrom                     #Moves the bin before from last iteration to two bins before
      TwoBinBeforeRC <- BinBeforeRC
      }
      
      if (Line > 3){
      BinBeforeChrom <- ProcessedBinChrom                     #Moves the bin  from last iteration to the bin before
      BinBeforeRC <- ProcessedBinRC
      }
      
      if (Line > 2){
      ProcessedBinChrom <- BinAfterChrom                      #Moves the bin after from last iteration to the processed bin
      ProcessedBinRC <- BinAfterRC
      ProcessedBinPos <- BinAfterPos
      }
      
      if (Line > 1){
      BinAfterChrom <- TwoBinAfterChrom                       #Moves the bin two after from last iteration to the bin after
      BinAfterRC <- TwoBinAfterRC
      BinAfterPos <- TwoBinAfterPos
      } 
      #Gets information for all lines 
      if (Line < (Endpos + 1)){
        #Transfers non-numerical chromosome numbers to numerical ones X -> 23, Y -> 24, MT and GL -> 25
        if ( suppressWarnings(is.na(as.numeric(as.matrix(Chrom)[Line,]))) ){
          if (as.matrix(Chrom)[Line,] == "X"){
            TwoBinAfterChrom <- 23
          } else if (as.matrix(Chrom)[Line,] == "Y"){
            TwoBinAfterChrom <- 24
          } else {
            TwoBinAfterChrom <- 25
          }
        } else 
        TwoBinAfterChrom <- as.numeric(as.matrix(Chrom)[Line,]) #Gets the chromosome of the bin of the current line read in the
  
        TwoBinAfterRC <- as.numeric(as.matrix(RC)[Line,])       #idem for Read count
        TwoBinAfterPos <- as.numeric(as.matrix(Pos)[Line,])     #idem for Chromosomal Position
      }
      

      #Calculate read count for surrounding bins of first bin of each chromosome
      if (ProcessedBinChrom != BinBeforeChrom){
        RCAround <- sum(BinAfterRC,TwoBinAfterRC)
        message("Processing chromosome ", ProcessedBinChrom)
      }
      
      #Calculate read count for surrounding bins of second bin of each chromosome
      if (ProcessedBinChrom == BinBeforeChrom && ProcessedBinChrom != TwoBinBeforeChrom ){
        RCAround <- sum(BinBeforeRC,BinAfterRC,TwoBinAfterRC)
      }
      
      #Calculate read count for surrounding bins of bins in the middle of each chromosome
      if (ProcessedBinChrom == TwoBinBeforeChrom && ProcessedBinChrom == TwoBinAfterChrom ){
        RCAround <- sum(TwoBinBeforeRC,BinBeforeRC,BinAfterRC,TwoBinAfterRC)
      }
      
      #Calculate read count for surrounding bins of penultimate bin of each chromosome
      if (ProcessedBinChrom == BinAfterChrom && ProcessedBinChrom != TwoBinAfterChrom ){
        RCAround <- sum(TwoBinBeforeRC,BinBeforeRC,BinAfterRC)
      }
      
      #Calculate read count for surrounding bins of last bin of each chromosome
      if (ProcessedBinChrom != BinAfterChrom){
        RCAround <- sum(TwoBinBeforeRC,BinBeforeRC)
      }
      
      #Correct read count based upon the read count of the surrounding bins
      if (RCAround == 0){
        CorrectedProcessedBinRC <- 0 
      } else
        CorrectedProcessedBinRC <- ProcessedBinRC

      #COllect chromosome, position and corrected read counts and put them in a list per position.      
      if (ProcessedBinChrom > 0) {
        CorrectedBin <- (c(as.integer(format(ProcessedBinChrom, scientific = FALSE)), as.integer(format(ProcessedBinPos, scientific = FALSE)), as.integer(format(CorrectedProcessedBinRC, scientific = FALSE))))
        CorrectedBinsList[[as.numeric(Line-2)]] <- CorrectedBin  #Fill position
      }
      
    }
    
    CorrectedBinslist_df <- data.frame(t(as.data.frame(CorrectedBinsList)))
    rownames(CorrectedBinslist_df) <- 1:nrow(CorrectedBinslist_df)
    return(CorrectedBinslist_df) #Note sex chromosomes and other genomic regions are kept in format 23, 24, 25
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
AnalysisFile <- read.csv(file = sample_name, sep = "\t", header = FALSE)


message ("Step 1 Two Besides Filter")
First_correction_TwoBesidesFilter  <- TwoBesideFilter(AnalysisFile)
write.table(x = First_correction_TwoBesidesFilter, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
message ("Done")
