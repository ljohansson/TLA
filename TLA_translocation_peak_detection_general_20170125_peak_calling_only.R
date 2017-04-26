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


############################FUNCTION MAKE LIST WITH PEAKS FROM START TO STOP###################################
PeakCalling <- function(BinsList){
  message("Calling peaks: covered bins")
  #Put list in dataframe for easy processing
  BinsList_df <- na.fail(as.data.frame(BinsList))
  
  #Put rows into separate variables
  Chrom <- (BinsList_df[,1, drop=FALSE])
  Pos <- (BinsList_df[,2, drop=FALSE])
  RC <- (BinsList_df[,3, drop=FALSE])
  
  #Prepare list for output
  PeaksList <- list()
  
  #Setting the chromosome checks for the surrounding bins for the loop. 
  BinBeforeChrom <- 0
  ProcessedBinChrom <- 0

  #Setting the read count checks for the surrounding bins for the loop
  BinBeforeRC <- 0
  ProcessedBinRC <- 0

  #Setting the position for the processed bin.
  BinBeforePos <- 0
  ProcessedBinPos <- 0

  #Set counter and variables to determine peak width
  Counter <- 0
  PeakRowCounter <- 0

  #Setting start and endline
  Startpos <- 1
  Endpos <- nrow(BinsList_df)
  
  #Setting bin size
  BinSize <- BinsList[[2]][2] - BinsList[[1]][2] + 1

  #Setting peak start and stop position
  PeakStartChrom <- 0
  PeakEndChrom <- 0
  PeakStartPos <- 0 
  PeakEndPos <- 0
  PeakRC <- 0
  
  #Looping through all rows to assess if they have a read count. If yes then the first covered bin is set as start and the last as end of a peak.
  for (Line in Startpos:Endpos){

    if (BinBeforeChrom != ProcessedBinChrom){
      }
    
    #Check if a peak has ended
    if ((BinBeforeRC > 0) && ((ProcessedBinRC == 0) || (BinBeforeChrom != ProcessedBinChrom))){
      PeakEndPos <- ProcessedBinPos 
      PeakEndChrom <- BinBeforeChrom
    } 
    if (PeakEndChrom > 0){
      PeakRowCounter = PeakRowCounter + 1
      PeaksList[[PeakRowCounter]] <- c(PeakStartChrom, PeakStartPos, PeakEndPos, PeakRC) #fill peaks list.
      Counter <- 0          #reset counters after peak to enable start of new peak
      PeakEndChrom <- 0     #reset counters after peak to enable start of new peak
    }
    

    #Check if a peak starts or continues and continue peak or start new one.
    if ( ProcessedBinRC > 0) {
      Counter <- Counter + 1
      if ( (BinBeforeRC == 0) || (BinBeforeChrom != ProcessedBinChrom) ){
        PeakStartChrom <- ProcessedBinChrom
        PeakStartPos <- ProcessedBinPos
        PeakRC <- ProcessedBinRC
      } else PeakRC <- paste(PeakRC, ProcessedBinRC, sep = ";")
    }
    

    #Gets information for all lines 
    if (Line < (Endpos + 1)){
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
      
      BinBeforeChrom <- ProcessedBinChrom
      BinBeforePos <- ProcessedBinPos
      BinBeforeRC <- ProcessedBinRC
      ProcessedBinChrom <- as.numeric(as.matrix(Chrom)[Line,]) #Gets the chromosome of the bin of the current line read in the dataframe
      ProcessedBinRC <- as.numeric(as.matrix(RC)[Line,])       #idem for Read count
      ProcessedBinPos <- as.numeric(as.matrix(Pos)[Line,])     #idem for Chromosomal Position
    }
  }
  #Part 2 of script
  #Merge closeby peaks (10 bin sizes) into a single peak.
  NumberOfCalls <- length(PeaksList)
  message ("Number of unmerged peak calls is ", NumberOfCalls)
  
  #Set variables
  PeaksListMerge <- list()
  LastLine <- 0
  ProcessedLine <- 0 
  Counter <- 0
  Counter_PeakMergedList <- 1
  PeakStartChrom <- 0 
  PeakStartPos <- 0
  PeakEndPos <- 0
  CurrentLinePeakRC <- 0
  PeakRCTotal <- 0

  
  ##MERGING PEAKS##
  #Loop through all lines in peaks list and determine if they belong to the same peak ('split peak')
  for (Line in PeaksList){
    #Set start values for first line.
    Counter <- Counter + 1
    PreviousPeakStartPos <- PeakStartPos
    PreviousPeakStartChrom <- PeakStartChrom
    PreviousPeakEndPos <- PeakEndPos
    PreviousPeakRC <- CurrentLinePeakRC
    CurrentLineStartChrom <- Line[1]
    CurrentLineStartPos <- as.numeric(Line[2])
    CurrentLineEndPos <- as.numeric(Line[3])
    CurrentLinePeakRC <- Line[4]
    message ("Counter is ", Counter)
    #Set first bin of list
    if ( Counter == 1 ){
    PeakStartChrom <- CurrentLineStartChrom
    PeakStartPos <- CurrentLineStartPos
    PeakEndPos <- CurrentLineEndPos
    PeakRCTotal <- CurrentLinePeakRC
    PeakStartChromWrite <- PeakStartChrom 
    PeakStartPosWrite <- PeakStartPos
    PeakEndPosWrite <- PeakEndPos
    PeakRCTotalWrite <- PeakRCTotal
    Counter_Peak <- 0
    }

    
    #Fill existing peak  
    if ( (PreviousPeakEndPos < CurrentLineStartPos ) && (CurrentLineStartChrom == PreviousPeakStartChrom) ){
      if ((CurrentLineStartPos - PreviousPeakEndPos) < ( BinSize * 10 ) ) {
        #Add zero's as Read count for lines in a peak without coverage
        CurrentLineStartPosInter <- CurrentLineStartPos
          while ((CurrentLineStartPosInter - PreviousPeakEndPos) >= BinSize) {
            PeakRCTotal <- paste(PeakRCTotal,0, sep = ";") #Add the zero's
            CurrentLineStartPosInter <- (CurrentLineStartPosInter - BinSize)
          } 
        PeakEndPos <- CurrentLineEndPos
        PeakRCTotal <- paste(PeakRCTotal,CurrentLinePeakRC, sep = ";") #Add peak counts to list
        #Set the values for writing to output
        PeakStartChromWrite <- PeakStartChrom
        PeakStartPosWrite <- PeakStartPos
        PeakEndPosWrite <- (PeakEndPos - BinSize)
        PeakRCTotalWrite <- PeakRCTotal
        }
      }

    
    #Set first bin of new peak
    if ( ( PreviousPeakEndPos > CurrentLineStartPos ) || (CurrentLineStartPos - PreviousPeakEndPos) >= ( BinSize * 10 ) ){
      PeakStartChrom <- CurrentLineStartChrom
      PeakStartPos <- CurrentLineStartPos
      PeakEndPos <- CurrentLineEndPos
      PeakRCTotal <- CurrentLinePeakRC
      Counter_Peak <- Counter_Peak + 1 #Count the number of peaks that will end up in the output
    }
    
    #Check if it is the last bin of a peak and if yes write it to outputlist.
    if ((PeakStartPos != PreviousPeakStartPos) && (Counter_Peak > 0) && (PreviousPeakStartPos != 0) || (Counter == NumberOfCalls)) {
      message ("PeakStartPos is ", PeakStartPos, " and PreviousPeakStartPos is ", PreviousPeakStartPos )
        PeaksListMerge[[Counter_PeakMergedList]] <- c(PeakStartChromWrite, PeakStartPosWrite, PeakEndPosWrite, PeakRCTotalWrite)
        Counter_PeakMergedList <- Counter_PeakMergedList + 1
        message ("peakRCTotalWrite is ", PeakRCTotalWrite)
      
    }
    
    #Reset write for first bin of new peak
    if ( ( PreviousPeakEndPos > CurrentLineStartPos ) || (CurrentLineStartPos - PreviousPeakEndPos) >= ( BinSize * 10 )){
      PeakStartChromWriteOld <- PeakStartChromWrite
      PeakStartPosWriteOld <- PeakStartPosWrite
      PeakRCTotalWriteOld <- PeakRCTotalWrite
      PeakStartChromWrite <- CurrentLineStartChrom
      PeakStartPosWrite <- CurrentLineStartPos
      PeakEndPosWrite <- CurrentLineEndPos
      PeakRCTotalWrite <- CurrentLinePeakRC
      if ((Counter == NumberOfCalls) && !(PeakRCTotalWrite == PeakRCTotalWriteOld) && !(PeakStartChromWrite == PeakStartChromWriteOld) && !(PeakStartPosWrite == PeakStartPosWriteOld)) {  
        PeaksListMerge[[Counter_PeakMergedList]] <- c(PeakStartChromWrite, PeakStartPosWrite, PeakEndPosWrite, PeakRCTotalWrite)
        Counter_PeakMergedList <- Counter_PeakMergedList + 1
        message ("peakRCTotalWrite is ", PeakRCTotalWrite)
      }
    }
  
  }
  
  message ("There are ", Counter_Peak, " peaks found.")
  
  #Create output
  options(scipen=999)     #Disable scientific notation
  PeaksListMerge_df  <- data.frame(t(as.data.frame(PeaksListMerge))) 
  rownames(PeaksListMerge_df ) <- 1:nrow(PeaksListMerge_df ) 
  return(PeaksListMerge_df ) #Note sex chromosomes and other genomic regions are kept in format 23, 24, 25
}


############################FUNCTION FILTER SMALL PEAKS FROM LIST###################################
#Peaks are considered small if the total peak area is less than one percent of the largest peak in the list 
#or if the maximum height of the peak is less than one percent of the highest point of the highest peak

FilterSmallPeaksFromList <- function(PeakFileToFilter){
  message("Filtering noise peaks")
  #Put list in dataframe for easy processing
  ###PeakFileToFilter_df <- na.fail(as.data.frame(PeakFileToFilter))
  
  #Put rows into separate variables
  ChromCol <- (PeakFileToFilter[1])
  StartPosCol <- (PeakFileToFilter[2])
  EndPosCol <- (PeakFileToFilter[3])
  RCCol <- (PeakFileToFilter[4])
  
  #Set variables
  HighestPeakPointsA <- list()
  PeakAreasA <- list()
  HighestPeakPointsB <- list()
  PeakAreasB <- list()
  NumberOfPeaks <- nrow(ChromCol)

  #Check all peaks for higest value and peak area
  for (Line in 1:NumberOfPeaks){
    #Get info from rows
    Chrom <- as.numeric(ChromCol[Line,1])
    StartPos <- as.numeric(StartPosCol[Line,1])
    EndPos <- as.numeric(EndPosCol[Line,1])
    RC <- RCCol[Line,1]

    #Get values for individual peaks
    RCList <- strsplit(as.character(RC),";")                                                #Split ReadCounts to list
    RCListLength <- lengths(RCList)                                                         #Calculate number of rows in list for current peak
    HighestPeakPointsA[[Line]] <- lapply(RCList, function(x) x[which.max(as.numeric(x))])   #Get maximum value from list for current peak
    PeakAreasA[[Line]] <- lapply(RCList, function(x) sum(as.numeric(x)))                    #Get total peak area from list for current peak
  }
  
  #Reformat items in list so I can work with them
  for (Item in 1:NumberOfPeaks){
    HighestPeakPointsB[[Item]] <- unlist(as.numeric(HighestPeakPointsA[[Item]]))
    PeakAreasB[[Item]] <- unlist(as.numeric(PeakAreasA[[Item]]))
  }

  #Get highest peak overall
  HighestPeakIndex <- which.max(HighestPeakPointsB)
  LargestAreaIndex <- which.max(PeakAreasB)
  HighestPeak <- HighestPeakPointsB[[HighestPeakIndex]]
  LargestArea <- PeakAreasB[[LargestAreaIndex]]
  message ("Highest Peak is ",HighestPeak)
  message ("Largest Area is ",LargestArea)
  
  #Write peaks having the highest point higher than 0.3% of the highest peak and having at least 0.3% of the area of the largest peak.
  LargePeakList <- list()
  Counter <- 1
  for (Line in 1:NumberOfPeaks){
    #Get info from rows
    Chrom <- as.numeric(ChromCol[Line,1])
    StartPos <- as.numeric(StartPosCol[Line,1])
    EndPos <- as.numeric(EndPosCol[Line,1])
    RC <- RCCol[Line,1]

    if ((HighestPeakPointsB[[Line]] > ( HighestPeak * 0.003 )) && (PeakAreasB[[Line]] > ( LargestArea * 0.004 ))){  #Maybe set Highest peak to 0.002?
      LargePeakList[[Counter]] <- c(Chrom,StartPos,EndPos,as.character(RC))
      Counter <- Counter + 1
    }
  }
  
  message ("There are ", ( Counter - 1 ), " peaks left.")
  
  #Create output
  options(scipen=999)     #Disable scientific notation
  LargePeakList_df  <- data.frame(t(as.data.frame(LargePeakList))) 
  rownames(LargePeakList_df ) <- 1:nrow(LargePeakList_df ) 
  return(LargePeakList_df ) #Note sex chromosomes and other genomic regions are kept in format 23, 24, 25
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
BamFile_name <- paste(substr(sample_name, 1, regexpr('ROI', sample_name ) - 2 ),".bam", sep="") #For logic see https://www.r-bloggers.com/basic-text-string-functions-in-r/

message ("Step 0-7 already done")

#message ("Step 7 Peak Filter")
PeakFilter_correction  <- read.csv(file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt",sep=""), sep = "\t", header = FALSE)
#PeakFilter_correction  <- PeakFilter(Second_correction_TwoBesidesFilter)
#write.table(x = PeakFilter_correction, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
#message ("Done")


message ("Step 8 Peak calling")
Called_Peaks  <- PeakCalling(PeakFilter_correction)
write.table(x = Called_Peaks, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak_Call.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)

message ("Step 9 Peak calling filter")
Called_Peaks_File <- read.csv(file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak_Call.txt",sep=""), sep = "\t", header = FALSE)
Called_Peaks_Filtered  <- FilterSmallPeaksFromList(Called_Peaks_File)
write.table(x = Called_Peaks_Filtered, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak_Call_filtered.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
