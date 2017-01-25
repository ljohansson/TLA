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




############################FUNCTION THREE OUT OF FOUR TWO BESIDE FILTER###################################
LowestThreeFour_TwoBesideFilter <- function(BinsList){
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

  #Looping through all rows to assess if the coverage of the lowest two of all surrounding bins is 0. Note that the bin that is processed lies two behind the current line read.
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
    
    
    #Calculate read count of two bins with lowest read count for surrounding bins of first bin of each chromosome
    if (ProcessedBinChrom != BinBeforeChrom){
      RCAround <- sum(BinAfterRC,TwoBinAfterRC)
      message("Processing chromosome ", ProcessedBinChrom)
    }
    
    #Calculate read count of three bins with lowest read countfor surrounding bins of second bin of each chromosome
    if (ProcessedBinChrom == BinBeforeChrom && ProcessedBinChrom != TwoBinBeforeChrom ){
      RCAround <- sort(c(BinBeforeRC,BinAfterRC,TwoBinAfterRC))[1] + sort(c(BinBeforeRC,BinAfterRC,TwoBinAfterRC))[2] + sort(c(BinBeforeRC,BinAfterRC,TwoBinAfterRC))[3]
    }
    
    #Calculate read count of two bins with lowest read countfor surrounding bins of bins in the middle of each chromosome
    if (ProcessedBinChrom == TwoBinBeforeChrom && ProcessedBinChrom == TwoBinAfterChrom ){
      RCAround <- sort(c(TwoBinBeforeRC,BinBeforeRC,BinAfterRC,TwoBinAfterRC))[1] + sort(c(TwoBinBeforeRC,BinBeforeRC,BinAfterRC,TwoBinAfterRC))[2] + sort(c(TwoBinBeforeRC,BinBeforeRC,BinAfterRC,TwoBinAfterRC))[3]
    }
    
    #Calculate read count of two bins with lowest read countfor surrounding bins of penultimate bin of each chromosome
    if (ProcessedBinChrom == BinAfterChrom && ProcessedBinChrom != TwoBinAfterChrom ){
      RCAround <- sort(c(TwoBinBeforeRC,BinBeforeRC,BinAfterRC))[1] + sort(c(TwoBinBeforeRC,BinBeforeRC,BinAfterRC))[2] + sort(c(TwoBinBeforeRC,BinBeforeRC,BinAfterRC))[3]
    }
    
    #Calculate read count of two bins with lowest read count for surrounding bins of last bin of each chromosome
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


#######################FUNCTION FILTER IF SMALLER THAN THRESHOLD###########################
FilterLowCoverage <- function(BinsList, Threshold){
  message ("Filtering bins")

  
  #Put list in dataframe for easy processing
  BinsList_df <- na.fail(as.data.frame(BinsList))
  
  #Put rows into separate variables
  Chrom <- (BinsList_df[,1, drop=FALSE])
  Pos <- (BinsList_df[,2, drop=FALSE])
  RC <- (BinsList_df[,3, drop=FALSE])
  
  #Prepare list for output
  FilteredBinsList <- list()
  
  #Setting start and endline
  Startpos <- 1
  Endpos <- nrow(BinsList_df)

  #Looping through all rows to remove read count in bins with read count below threshold
  for (Line in Startpos:Endpos){                          #Gets information for all lines 
    
    ProcessedBinChrom <- as.numeric(as.matrix(Chrom)[Line,]) #Gets the chromosome of the bin of the current line read in the
    ProcessedBinPos <- as.numeric(as.matrix(Pos)[Line,])     #idem for Chromosomal Position
    ProcessedBinRC <- as.numeric(as.matrix(RC)[Line,])       #idem for Read count
    
    #Set read count to 0 for bins having a read count below the threshold
    if (ProcessedBinRC < Threshold){
      FilteredBinRC <- 0
    } else 
      FilteredBinRC <- ProcessedBinRC
    
    #Filling rows for filtered bin
    FilteredBin <- (c(as.integer(format(ProcessedBinChrom, scientific = FALSE)), as.integer(format(ProcessedBinPos, scientific = FALSE)), as.integer(format(FilteredBinRC, scientific = FALSE))))
    FilteredBinsList[[as.numeric(Line)]] <- FilteredBin  #Fill position
  }  
  
  #Creating output
  FilteredBinslist_df <- data.frame(t(as.data.frame(FilteredBinsList)))
  rownames(FilteredBinslist_df) <- 1:nrow(FilteredBinslist_df)
  return(FilteredBinslist_df) #Note sex chromosomes and other genomic regions are kept in format 23, 24, 25
  
}





#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

############################FUNCTION PEAK WITH FOUR OUT OF TEN MINIMUM WIDTH SIX###################################
PeakFilter <- function(BinsList){
  message("Detecting wide peaks: filtering bins")
  #Put list in dataframe for easy processing
  BinsList_df <- na.fail(as.data.frame(BinsList))
  
  #Put rows into separate variables
  Chrom <- (BinsList_df[,1, drop=FALSE])
  Pos <- (BinsList_df[,2, drop=FALSE])
  RC <- (BinsList_df[,3, drop=FALSE])
  
  #Prepare list for output
  CorrectedBinsList <- list()
  
  #Setting the chromosome checks for the surrounding bins for the loop. 
  ProcessedBinChrom <- 0
  OneBinAfterChrom <- 0
  TwoBinAfterChrom <- 0
  ThreeBinAfterChrom <- 0
  FourBinAfterChrom <- 0
  FiveBinAfterChrom <- 0
  SixBinAfterChrom <- 0
  SevenBinAfterChrom <- 0
  EightBinAfterChrom <- 0
  NineBinAfterChrom <- 0

  
  #Setting the read count checks for the surrounding bins for the loop
  ProcessedBinRC <- 0
  OneBinAfterRC <- 0
  TwoBinAfterRC <- 0
  ThreeBinAfterRC <- 0
  FourBinAfterRC <- 0
  FiveBinAfterRC <- 0
  SixBinAfterRC <- 0
  SevenBinAfterRC <- 0
  EightBinAfterRC <- 0
  NineBinAfterRC <- 0


  #Setting the position for the processed bin. The information of the two bins thereafter is needed to fill the processed bin.
  ProcessedBinPos <- 0
  OneBinAfterPos <- 0
  TwoBinAfterPos <- 0
  ThreeBinAfterPos <- 0
  FourBinAfterPos <- 0
  FiveBinAfterPos <- 0
  SixBinAfterPos <- 0
  SevenBinAfterPos <- 0
  EightBinAfterPos <- 0
  NineBinAfterPos <- 0
  
  #Setting bin status 0 is no coverage, 1 is coverage
  ProcessedBinStatus <- 0
  OneBinAfterStatus <- 0
  TwoBinAfterStatus <- 0
  ThreeBinAfterStatus <- 0
  FourBinAfterStatus <- 0
  FiveBinAfterStatus <- 0
  SixBinAfterStatus <- 0
  SevenBinAfterStatus <- 0
  EightBinAfterStatus <- 0
  NineBinAfterStatus <- 0
  
  #Set counter and variables to determine peak width
  Counter <- 0
  startpos <- 1
  stoppos <- 0
  
  #Setting start and endline
  Startpos <- 1
  Endpos <- nrow(BinsList_df)
  
  #Looping through all rows to assess if five of ten bins have a value. If yes then the bins retain their value. If  Note that the bin that is processed lies two behind the current line read.
  for (Line in Startpos:(Endpos + 9)){
    #Reset Counter
    Counter <- 0
    
    if (Line > 9){
      ProcessedBinChrom <- OneBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      ProcessedBinRC <- OneBinAfterRC
      ProcessedBinPos <- OneBinAfterPos
      ProcessedBinStatus <- OneBinAfterStatus
    }
    
    if (Line > 8){
      OneBinAfterChrom <- TwoBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      OneBinAfterRC <- TwoBinAfterRC
      OneBinAfterPos <- TwoBinAfterPos
      OneBinAfterStatus <- TwoBinAfterStatus
    }
    
    if (Line > 7){
      TwoBinAfterChrom <- ThreeBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      TwoBinAfterRC <- ThreeBinAfterRC
      TwoBinAfterPos <- ThreeBinAfterPos
      TwoBinAfterStatus <- ThreeBinAfterStatus
    }
    
    if (Line > 6){
      ThreeBinAfterChrom <- FourBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      ThreeBinAfterRC <- FourBinAfterRC
      ThreeBinAfterPos <- FourBinAfterPos
      ThreeBinAfterStatus <- FourBinAfterStatus
    }
    
    if (Line > 5){
      FourBinAfterChrom <- FiveBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      FourBinAfterRC <- FiveBinAfterRC
      FourBinAfterPos <- FiveBinAfterPos
      FourBinAfterStatus <- FiveBinAfterStatus
    }
    
    if (Line > 4){
      FiveBinAfterChrom <- SixBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      FiveBinAfterRC <- SixBinAfterRC
      FiveBinAfterPos <- SixBinAfterPos
      FiveBinAfterStatus <- SixBinAfterStatus
    }
    
    if (Line > 3){
      SixBinAfterChrom <- SevenBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      SixBinAfterRC <- SevenBinAfterRC
      SixBinAfterPos <- SevenBinAfterPos
      SixBinAfterStatus <- SevenBinAfterStatus
    }
    
    if (Line > 2){
      SevenBinAfterChrom <- EightBinAfterChrom                    #Moves the bin after from last iteration to the processed bin
      SevenBinAfterRC <- EightBinAfterRC
      SevenBinAfterPos <- EightBinAfterPos
      SevenBinAfterStatus <- EightBinAfterStatus
    }
    
    if (Line > 1){
      EightBinAfterChrom <- NineBinAfterChrom                         #Moves the bin after from last iteration to the processed bin
      EightBinAfterRC <- NineBinAfterRC
      EightBinAfterPos <- NineBinAfterPos
      EightBinAfterStatus <- NineBinAfterStatus
    }
    
 
    #Gets information for all lines 
    if (Line < (Endpos + 1)){
      #Transfers non-numerical chromosome numbers to numerical ones X -> 23, Y -> 24, MT and GL -> 25
      if ( suppressWarnings(is.na(as.numeric(as.matrix(Chrom)[Line,]))) ){
        if (as.matrix(Chrom)[Line,] == "X"){
          NineBinAfterChrom <- 23
        } else if (as.matrix(Chrom)[Line,] == "Y"){
          NineBinAfterChrom <- 24
        } else {
          NineBinAfterChrom <- 25
        }
      } else 
      NineBinAfterChrom <- as.numeric(as.matrix(Chrom)[Line,]) #Gets the chromosome of the bin of the current line read in the
      NineBinAfterRC <- as.numeric(as.matrix(RC)[Line,])       #idem for Read count
      NineBinAfterPos <- as.numeric(as.matrix(Pos)[Line,])     #idem for Chromosomal Position
      NineBinAfterStatus <- 0
    }
    
    #Calculate series of covered bins. If at least four bins out of ten have coverage and the last bin is at least five bins further than the read counts should be retained. 
    #Otherwise they should be filtered.
    
    #put read count status to one if there are enough covered bins in the region in a wide enough peak
    if (ProcessedBinRC > 0){
      Counter <- 1
      if ((OneBinAfterRC > 0) && (OneBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((TwoBinAfterRC > 0) && (TwoBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((ThreeBinAfterRC > 0) && (ThreeBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((FourBinAfterRC > 0) && (FourBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((FiveBinAfterRC > 0) && (FiveBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((SixBinAfterRC > 0) && (SixBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((SevenBinAfterRC > 0) && (SevenBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((EightBinAfterRC > 0) && (EightBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      if ((NineBinAfterRC > 0) && (NineBinAfterChrom == ProcessedBinChrom)){
        Counter <- Counter + 1
      }
      
      #Determine peak width and set status to 1 (= peak) if the count of covered bins is at least 4 and the peak width is at least 6.
      if (Counter > 3){
        if (NineBinAfterRC > 0){
          stoppos <- 10
        } else if (EightBinAfterRC > 0){
          stoppos <- 9
        } else if (SevenBinAfterRC > 0){
          stoppos <- 8
        } else if (SixBinAfterRC > 0){
          stoppos <- 7
        } else if (FiveBinAfterRC > 0){
          stoppos <- 6
        } else { 
          stoppos <- 1
        }
        if (stoppos - startpos > 0){
          ProcessedBinStatus <- 1
          OneBinAfterStatus <- 1
          TwoBinAfterStatus <- 1
          ThreeBinAfterStatus <- 1
          FourBinAfterStatus <- 1
          FiveBinAfterStatus <- 1
          SixBinAfterStatus <- 1
          SevenBinAfterStatus <- 1
          EightBinAfterStatus <- 1
          NineBinAfterStatus <- 1
        }
      }
    }
    
    #Reset Binstatus at moving to next chromosome
    if (NineBinAfterChrom != ProcessedBinChrom){
      NineBinAfterStatus <- 0
      if (EightBinAfterChrom != ProcessedBinChrom){
        EightBinAfterStatus <- 0
      }
      if (SevenBinAfterChrom != ProcessedBinChrom){
        SevenBinAfterStatus <- 0
      }
      if (SixBinAfterChrom != ProcessedBinChrom){
        SixBinAfterStatus <- 0
      }
      if (FiveBinAfterChrom != ProcessedBinChrom){
        FiveBinAfterStatus <- 0
      }
      if (FourBinAfterChrom != ProcessedBinChrom){
        FourBinAfterStatus <- 0
      }
      if (ThreeBinAfterChrom != ProcessedBinChrom){
        ThreeBinAfterStatus <- 0
      }
      if (TwoBinAfterChrom != ProcessedBinChrom){
        TwoBinAfterStatus <- 0
      }
      if (OneBinAfterChrom != ProcessedBinChrom){
        OneBinAfterStatus <- 0
        message("Processing chromosome ", OneBinAfterChrom)
      }
    }
    
    #Correct read count based upon the bin status
    if (ProcessedBinStatus == 0){
      CorrectedProcessedBinRC <- 0 
    } else if (ProcessedBinStatus == 1) {
      CorrectedProcessedBinRC <- ProcessedBinRC
    }
    
    #Collect chromosome, position and corrected read counts and put them in a list per position.      
    if (ProcessedBinChrom > 0) {
      CorrectedBin <- (c(as.integer(format(ProcessedBinChrom, scientific = FALSE)), as.integer(format(ProcessedBinPos, scientific = FALSE)), as.integer(format(CorrectedProcessedBinRC, scientific = FALSE))))
      CorrectedBinsList[[as.numeric(Line-9)]] <- CorrectedBin  #Fill position
    }
  
  }  
  CorrectedBinslist_df <- data.frame(t(as.data.frame(CorrectedBinsList))) 
  rownames(CorrectedBinslist_df) <- 1:nrow(CorrectedBinslist_df) 
  return(CorrectedBinslist_df) #Note sex chromosomes and other genomic regions are kept in format 23, 24, 25
}





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






############################FUNCTION COUNT TOTAL ALIGNED READS BAM FILE###################################

CountAlignedReadsBam <- function(BamFile_for_counting){
###Dependency Rsamtools
message ("Loading Rsamtools")
library("Rsamtools")
  
#load bam and count reads
message ("Loading bam sample")

#Counting aligned reads in bamfile
AlignedReadsBamFile <- Rsamtools::countBam(file = BamFile_for_counting, index = paste(BamFile_for_counting,".bai", sep=""), param=ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isUnmappedQuery=FALSE) ) )$records

message ("Detaching Rsamtools")
detach("package:Rsamtools", unload = T)

return(AlignedReadsBamFile)
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

message ("Step 0 counting aligned read number total bam file")
Aligned_reads <- CountAlignedReadsBam(BamFile_name)
message ("The number of aligned reads are ", Aligned_reads)

message ("Step 1 Two Besides Filter")
First_correction_TwoBesidesFilter  <- TwoBesideFilter(AnalysisFile)
write.table(x = First_correction_TwoBesidesFilter, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
rm(AnalysisFile)
gc()

message ("Step 2 Merge to 100kb bins")
Mergedto100kbBins <- Merge_10_bins(First_correction_TwoBesidesFilter)
write.table(x = Mergedto100kbBins, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
rm(First_correction_TwoBesidesFilter)
gc()

message ("Step 3 Filter bins having coverage below 1 / 30000 the number of reads in the total bamfile")
FilterThreshold <- Aligned_reads / 30000
Filtered1000 <- FilterLowCoverage(Mergedto100kbBins,FilterThreshold)
write.table(x = Filtered1000, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
rm(Mergedto100kbBins)
gc()

message ("Step 4 Lowest two out of four Two Besides Filter")
Corrected_LowestThreeFour_TwoBesidesFilter <- LowestThreeFour_TwoBesideFilter(Filtered1000)
write.table(x = Corrected_LowestThreeFour_TwoBesidesFilter, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
rm(Filtered1000)
gc()

message ("Step 5 Filter bins having coverage below 1 / 6000 the number of reads in the total bamfile")
FilterThreshold <- Aligned_reads / 6000
Filtered5000 <- FilterLowCoverage(Corrected_LowestThreeFour_TwoBesidesFilter,FilterThreshold)
write.table(x = Filtered5000, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
rm(Corrected_LowestThreeFour_TwoBesidesFilter)
gc()

message ("Step 6 Two Besides Filter")
Second_correction_TwoBesidesFilter  <- TwoBesideFilter(Filtered5000)
write.table(x = Second_correction_TwoBesidesFilter, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
rm(Filtered5000)
gc()

message ("Step 7 Peak Filter")
PeakFilter_correction  <- PeakFilter(Second_correction_TwoBesidesFilter)
write.table(x = PeakFilter_correction, file = paste(tmpfolder,file_path_sans_ext(basename(sample_name)),"_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
message ("Done")



##TEST PEAK CALLING###
sample_name <- "C:/TLA/FilterQ10_20_30/35_MP1_161216_NB501043_0093_AHTJHMBGXY_L1234_TCCGGAGA-CTTCGCCT_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_4_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt"
message ("Analyzing sample ", sample_name)
AnalysisFile <- read.csv(file = sample_name, sep = "\t", header = FALSE)

message ("Peak calling")
Called_Peaks  <- PeakCalling(AnalysisFile)
write.table(x = Called_Peaks, file = paste("C:/TLA/FilterQ10_20_30/",file_path_sans_ext(basename(sample_name)),"_Call.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)


##TEST PEAKLIST FILTER###
sample_name <- "C:/TLA/FilterQ10_20_30/35_MP1_161216_NB501043_0093_AHTJHMBGXY_L1234_TCCGGAGA-CTTCGCCT_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_4_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak_Call.txt"
message ("Analyzing sample ", sample_name)
AnalysisFilePeakList <- read.table(file = sample_name, sep = "\t", header = FALSE)

message ("Peak calling")
Called_Peaks_Filtered  <- FilterSmallPeaksFromList(AnalysisFilePeakList)
write.table(x = Called_Peaks_Filtered, file = paste("C:/TLA/FilterQ10_20_30/",file_path_sans_ext(basename(sample_name)),"_filtered.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)


##TEST PEAK CALLING2##
sample_name <- "C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/4_MP1_160725_NB501043_0045_AHWJKCBGXX_L1234_ATTACTCG-ACGTCCTG_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_17_Peak.txt"
#sample_name <- "C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/5_MP1_160725_NB501043_0045_AHWJKCBGXX_L1234_TCCGGAGA-AGGCTATA_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_17_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt"
message ("Analyzing sample ", sample_name)
AnalysisFile <- read.csv(file = sample_name, sep = "\t", header = FALSE)

message ("Peak calling")
Called_Peaks  <- PeakCalling(AnalysisFile)
write.table(x = Called_Peaks, file = paste("C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/",file_path_sans_ext(basename(sample_name)),"_Call.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)

##TEST PEAKLIST FILTER2###
sample_name <- "C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/4_MP1_160725_NB501043_0045_AHWJKCBGXX_L1234_ATTACTCG-ACGTCCTG_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_17_Peak_Call.txt"
message ("Analyzing sample ", sample_name)
AnalysisFilePeakList <- read.table(file = sample_name, sep = "\t", header = FALSE)

message ("Peak calling")
Called_Peaks_Filtered  <- FilterSmallPeaksFromList(AnalysisFilePeakList)
write.table(x = Called_Peaks_Filtered, file = paste("C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/",file_path_sans_ext(basename(sample_name)),"_filtered.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)


##TEST PEAK CALLING3##
#sample_name <- "C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/2_MP1_160725_NB501043_0045_AHWJKCBGXX_L1234_ATTACTCG-AGGATAGG_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_17_Peak.txt"
#sample_name <- "C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/5_MP1_160725_NB501043_0045_AHWJKCBGXX_L1234_TCCGGAGA-AGGCTATA_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_17_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt"
sample_name <- "C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_16/17_MP1_160725_NB501043_0045_AHWJKCBGXX_L1234_ATTCAGAA-AGGCTATA_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_16_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt"
message ("Analyzing sample ", sample_name)
AnalysisFile <- read.csv(file = sample_name, sep = "\t", header = FALSE)

message ("Peak calling")
Called_Peaks  <- PeakCalling(AnalysisFile)
write.table(x = Called_Peaks, file = paste("C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/",file_path_sans_ext(basename(sample_name)),"_Call.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)

##TEST PEAKLIST FILTER3###
sample_name <- "C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/2_MP1_160725_NB501043_0045_AHWJKCBGXX_L1234_ATTACTCG-AGGATAGG_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_17_Peak_Call.txt"
message ("Analyzing sample ", sample_name)
AnalysisFilePeakList <- read.table(file = sample_name, sep = "\t", header = FALSE)

message ("Peak calling")
Called_Peaks_Filtered  <- FilterSmallPeaksFromList(AnalysisFilePeakList)
write.table(x = Called_Peaks_Filtered, file = paste("C:/TLA/analysis_reanalysis_Q30_Dilution2/Peak_files/ROI_17/",file_path_sans_ext(basename(sample_name)),"_filtered.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)




##TEST PEAK CALLING###
sample_name <- "C:/TLA/161216/bin_counts/ROI_4/51_MP1_161216_NB501043_0093_AHTJHMBGXY_L1234_GAATTCGT-CTTCGCCT_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_4_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt"
sample_name2 <- "C:/TLA/161216/bin_counts/ROI_14/43_MP1_161216_NB501043_0093_AHTJHMBGXY_L1234_GAGATTCC-CTTCGCCT_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_14_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak.txt"
message ("Analyzing sample ", sample_name)
AnalysisFile <- read.table(file = sample_name, sep = "\t", header = FALSE)

Called_Peaks  <- PeakCalling(AnalysisFile)
write.table(x = Called_Peaks, file = paste("C:/TLA/161216/Peak_files/",file_path_sans_ext(basename(sample_name)),"_Call.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)

message ("Analyzing sample ", sample_name2)
AnalysisFile2 <- read.table(file = sample_name2, sep = "\t", header = FALSE)

Called_Peaks2  <- PeakCalling(AnalysisFile2)
write.table(x = Called_Peaks2, file = paste("C:/TLA/161216/Peak_files/",file_path_sans_ext(basename(sample_name2)),"_Call.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)


  
##TEST PEAKLIST FILTER###
sample_name <- "C:/TLA/161216/Peak_files//ROI_14/43_MP1_161216_NB501043_0093_AHTJHMBGXY_L1234_GAGATTCC-CTTCGCCT_nodup_trimmed_digested_MP1and2_uniquely_aligned_ROI_14_coverage_2BesFilt_100kb_Filt1000_Low3_4TwoBes_Filt5000_2BesFil_Peak_Call.txt"
message ("Analyzing sample ", sample_name)
AnalysisFilePeakList <- read.table(file = sample_name, sep = "\t", header = FALSE)

message ("Peak calling")
Called_Peaks_Filtered  <- FilterSmallPeaksFromList(AnalysisFilePeakList)
write.table(x = Called_Peaks_Filtered, file = paste("C:/TLA/161216/Peak_files/ROI_5/",file_path_sans_ext(basename(sample_name)),"_filtered.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)

 