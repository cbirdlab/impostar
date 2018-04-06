#' Logistic Regression Test
#' 
#' @description This function performs a binomial logistic regression test for a genetic cline among three populations, implementing the appropriate resampling strategies to model the total sampling error associated with either Sanger sequencing (population sampling error) or next-generation sequencing (population + sequencer sampling error). 
#' @usage runLogRegTest(Input, OutputBase=NULL, mainDirectory=getwd(), NumBS=1000, ModelSeqSampErrSwitch=T, nCores=NULL)
#' @param Input A dataframe where each row is a SNP with a unique snpID, and columns represent the total number of alleles (2n, TotAlleles), read depth (DP), reference alleles (RefAl), alternate alleles (AltAl), reference reads (RD), and alternate reads (AD) per pool (as indicated by the number in the column names). Use simulate_data to generate an example and view required format. Required.
#' @param OutputBase  A character string indicating the base of the output filename, to which the number of bootstraps and if ModelSeqSampError was incorporated (ResampReads) will be appended. If NULL (default) then no output file is written.
#' @param mainDirectory A character string indicating the working directory to use for input and optional output files. Default is getwd().
#' @param NumBS An integer >=2 indicating the number of bootstrap iterations to perform. Default is 1000.
#' @param ModelSeqSampErrSwitch A logical indicating whether the resampling strategy for next-generation sequencing data (TRUE, Default) or Sanger sequencing data (FALSE) should be implemented.
#' @param nCores An integer indicating the number of cores to use in parallel processing, or NULL (default) to use the maximum available cores minus 1. Note that parallel processing not supported by Windows.
#' @return The runLogRegTest function returns a dataframe and optional csv file written to the working directory. Each row is a SNP, where ObsSlope refers to the observed rate of change in the log odds ratio of the logistic regression line, and SlopeP is used to evalute its significance.
#' @author  Rebecca M. Hamner, Jason D. Selwyn, Evan Krell, Scott A. King, Christopher E. Bird
#' @references Hamner, R.M., J.D. Selwyn, E. Krell, S.A. King, and C.E. Bird. In review. Modeling next-generation sequencer sampling error in pooled population samples dramatically reduces false positives in genetic structure tests.
#' @seealso \code{\link{simulate_data}}, \code{\link{runAMOVA}}
#' @keywords logistic regression, genetic cline, next-generation sequencing, SNPs, pool
#' @import dplyr parallel 
#' @export
#' @examples
#' # create example data file
#' Input <- simulate_data(rep(50, 3), rep(100, 3), rep(0.5, 3), 1000, file_name=T)
#' # run logistic regression test with default options
#' LogRegResults <- runLogRegTest(simdata)
#' # run logistic regression test specifying more options
#' LogRegResults <- runLogRegTest(Input, OutputBase=NULL, mainDirectory=getwd(), NumBS=1000, ModelSeqSampErrSwitch=T, nCores=NULL)


runLogRegTest <- function(Input, OutputBase=NULL, mainDirectory=getwd(), NumBS=1000, ModelSeqSampErrSwitch=T, nCores=NULL){
  # Input:	a dataframe where each row is a SNP with a unique snpID, and columns represent the total number of alleles (2n, TotAlleles), read depth (DP), reference alleles (RefAl), alternate alleles (AltAl), reference reads (RD), and alternate reads (AD) per pool (as indicated by the number in the column names). Use simulate_data to generate an example and view required format.
  # OutputBase:	 a character string indicating the base of the output filename, to which the number of bootstraps and if ModelSeqSampError was incorporated (ResampReads) will be appended. If NULL (default) then no output file is written.
  # mainDirectory: a character string indicating the working directory to use for input and optional output files. Default is getwd().
  # NumBS:	an integer >=2 indicating the number of bootstrap iterations to perform.
  # ModelSeqSampErrSwitch: 	a logical indicating whether the resampling strategy for next-generation sequencing data (TRUE, Default) or Sanger sequencing data (FALSE) should be implemented.
  # nCores:	an integer indicating the number of cores to use in parallel processing, or NULL (default) to use the maximum available cores minus 1. Note that parallel processing not supported by Windows.
  
  ################################
  # READ IN INPUT DATA & SETTINGS
  
  # SET number of iterations to run
  nBS <- NumBS  
  
  ModelSeqSampErr<-ModelSeqSampErrSwitch
  
  # Turn on (TRUE) or off (FALSE) boostrap to model sequencer sampling error
  if(ModelSeqSampErr==TRUE) {ResampYN <- "ResampReads_"} else {ResampYN <- "noResampReads_"}
  
  # Create and set output directory & subdirectory
  mainDir <- mainDirectory
  if (file.exists(mainDirectory)){
    setwd(mainDirectory)
  } else {
    dir.create(mainDirectory)
    setwd(mainDirectory)}
  
  # Read in input file 
  if(class(Input)=='character'){
    SNPset <- read.csv(Input, sep=",", header = TRUE, row.names=1, fill = TRUE)
  } else if (class(Input)=='data.frame'){
    SNPset <- Input
  } else {
    stop("Input needs to be either a data frame or character string of csv")
  }
  
  # View head to ensure proper read-in
  print(head(SNPset))
  
  # Create a variable for the number of pools (populations or groups being compared) by counting the number of times a column name begins with RD and any number
  pools <- length(grep("^RD[0-9]", colnames(SNPset), value=TRUE))  
  
  # Create variable for number of SNPs
  nSNPs <- nrow(SNPset)    
  
  ######################################################################################
  # CREATE FUNCTION Allelotyper
  # Converts reads to allele counts to allele frequencies
  
  Allelotyper <- function(x){    
    
    # CONVERT READS TO ALLELE COUNTS
    # make function to round >=0.5 up and <0.5 down
    rnd <- function(y) trunc(y+sign(y)*0.5)
    
    # SCALE the REFERENCE READS to ALLELES, for each pool
    x$RefAlleles1 <- rnd((x$RD1 / x$DP1) * x$TotAlleles1)
    x$RefAlleles2 <- rnd((x$RD2 / x$DP2) * x$TotAlleles2)
    x$RefAlleles3 <- rnd((x$RD3 / x$DP3) * x$TotAlleles3)
    
    # CALCULATE ALTERNATE ALLELES (TotAlleles - RefAlleles), for each pool      
    x$AltAlleles1 <- x$TotAlleles1 - x$RefAlleles1
    x$AltAlleles2 <- x$TotAlleles2 - x$RefAlleles2
    x$AltAlleles3 <- x$TotAlleles3 - x$RefAlleles3
    
    
    # CONVERT ALLELE COUNTS TO REF ALLELE FREQUENCIES
    # Create columns with Allele Frequencies for each pool
    x$RefFreq1 <- x$RefAlleles1/x$TotAlleles1
    x$RefFreq2 <- x$RefAlleles2/x$TotAlleles2
    x$RefFreq3 <- x$RefAlleles3/x$TotAlleles3
    
    return(x)      
    
  } # end of FUNCTION: Allelotyper
  ######################################################################################
  
  
  ######################################################################################
  # Create function LogisticRegScratch - use when modeling NGS sampling error
  # Calculates Logistic regression of RefFreqs (y) by shore height (x) for each SNP
  
  suppressPackageStartupMessages(library(dplyr))
  
  LogisticRegScratch <- function(Y, x){  
    
    # Reformat data
    RefSetToMelt <- x[Y,c("snpID", "RefAlleles1", "RefAlleles2", "RefAlleles3")]
    AltSetToMelt <- x[Y,c("snpID", "AltAlleles1", "AltAlleles2", "AltAlleles3")]
    
    AllelesReformated1 <- as.data.frame(rep(RefSetToMelt$snpID, x$TotAlleles1[1]))
    colnames(AllelesReformated1)[1] <- "snpID"
    AllelesReformated1$Alleles <-c(rep(1, RefSetToMelt$RefAlleles1), rep(0, AltSetToMelt$AltAlleles1))
    AllelesReformated1$Pool <- 1
    
    AllelesReformated2 <- as.data.frame(rep(RefSetToMelt$snpID,  x$TotAlleles2[1]))
    colnames(AllelesReformated2)[1] <- "snpID"
    AllelesReformated2$Alleles <-c(rep(1, RefSetToMelt$RefAlleles2), rep(0, AltSetToMelt$AltAlleles2))
    AllelesReformated2$Pool <- 2
    
    AllelesReformated3 <- as.data.frame(rep(RefSetToMelt$snpID,  x$TotAlleles3[1]))
    colnames(AllelesReformated3)[1] <- "snpID"
    AllelesReformated3$Alleles <-c(rep(1, RefSetToMelt$RefAlleles3), rep(0, AltSetToMelt$AltAlleles3))
    AllelesReformated3$Pool <- 3
    
    AllelesReformated <- rbind(AllelesReformated1, AllelesReformated2, AllelesReformated3)
    
    ########################################
    # Define function LogisticRegScratchGuts  #adapted from http://pingax.com/logistic-regression-r-step-step-implementation-part-2/
    
    LogisticRegScratchGuts <- function(X, y){
      #Predictor variable
      X <- as.matrix(X)
      # Add ones to X
      X <- cbind(rep(1,nrow(X)),X)
      # Response variable
      Y <- as.matrix(y)
      
      # sigmoid function
      sigmoid <- function(z)
      {
        g <- 1/(1+exp(-z))
        return(g)
      } # end sigmoid function
      
      # Cost Function
      cost <- function(theta)
      {
        m <- nrow(X)
        g <- sigmoid(X%*%theta)
        J <- (1/m)*sum((-Y*log(g)) - ((1-Y)*log(1-g)))
        return(J)
      } # end cost function
      
      # Intial theta
      initial_theta <- rep(0,ncol(X))
      
      # Derive theta using gradient descent using optim function
      theta_optim <- optim(par = initial_theta, fn = cost)
      
      # pull out theta ("slope"; rate of change of log odds ratio)
      theta <- theta_optim$par
      return(theta[2])
      
    } # End of FUNCTION: LogisticRegScratchGuts
    ###########################################
    
    # Run function LogisticRegScratchGuts
    TempList <- list()
    TempList[Y] <- LogisticRegScratchGuts(AllelesReformated$Pool,AllelesReformated$Alleles)
    TempUnlisted <- unlist(TempList)
    
    return(TempUnlisted)
  } #end function LogisticRegScratch
  ######################################################################################
  
  
  ######################################################################################
  # Create function LogisticRegScratch_noResampReads - use when NOT modeling NGS sampling error
  # Calculates Logistic regression of RefFreqs (y) by shore height (x) for each SNP
  
  LogisticRegScratch_noResampReads <- function(Y, x){  
    
    # Reformat data        
    RefSetToMelt <- x[Y,c("snpID", "shufRefAlleles1", "shufRefAlleles2", "shufRefAlleles3")]
    AltSetToMelt <- x[Y,c("snpID", "shufAltAlleles1", "shufAltAlleles2", "shufAltAlleles3")]
    
    AllelesReformated1 <- as.data.frame(rep(RefSetToMelt$snpID, x$TotAlleles1[1]))
    colnames(AllelesReformated1)[1] <- "snpID"
    AllelesReformated1$Alleles <-c(rep(1, RefSetToMelt$shufRefAlleles1), rep(0, AltSetToMelt$shufAltAlleles1))
    AllelesReformated1$Pool <- 1
    
    AllelesReformated2 <- as.data.frame(rep(RefSetToMelt$snpID,  x$TotAlleles2[1]))
    colnames(AllelesReformated2)[1] <- "snpID"
    AllelesReformated2$Alleles <-c(rep(1, RefSetToMelt$shufRefAlleles2), rep(0, AltSetToMelt$shufAltAlleles2))
    AllelesReformated2$Pool <- 2
    
    AllelesReformated3 <- as.data.frame(rep(RefSetToMelt$snpID,  x$TotAlleles3[1]))
    colnames(AllelesReformated3)[1] <- "snpID"
    AllelesReformated3$Alleles <-c(rep(1, RefSetToMelt$shufRefAlleles3), rep(0, AltSetToMelt$shufAltAlleles3))
    AllelesReformated3$Pool <- 3
    
    AllelesReformated <- rbind(AllelesReformated1, AllelesReformated2, AllelesReformated3)
    
    
    #########################################
    # Define function LogisticRegScratchGuts  #adapted from http://pingax.com/logistic-regression-r-step-step-implementation-part-2/
    
    LogisticRegScratchGuts <- function(X, y){
      #Predictor variable
      X <- as.matrix(X)
      # Add ones to X
      X <- cbind(rep(1,nrow(X)),X)
      # Response variable
      Y <- as.matrix(y)
      
      # sigmoid function
      sigmoid <- function(z)
      {
        g <- 1/(1+exp(-z))
        return(g)
      } # end sigmoid function
      
      # Cost Function
      cost <- function(theta)
      {
        m <- nrow(X)
        g <- sigmoid(X%*%theta)
        J <- (1/m)*sum((-Y*log(g)) - ((1-Y)*log(1-g)))
        return(J)
      } # end cost function
      
      # Intial theta
      initial_theta <- rep(0,ncol(X))
      
      # Derive theta using gradient descent using optim function
      theta_optim <- optim(par = initial_theta, fn = cost)
      
      # pull out theta ("slope"; rate of change of log odds ratio)
      theta <- theta_optim$par
      return(theta[2])
    } # End of FUNCTION: LogisticRegScratchGuts
    ###########################################
    
    # Run function LogisticRegScratchGuts
    TempList <- list()
    TempList[Y] <- LogisticRegScratchGuts(AllelesReformated$Pool,AllelesReformated$Alleles)
    TempUnlisted <- unlist(TempList)
    
    return(TempUnlisted)
  } #end function LogisticRegScratch_noResampReads
  ######################################################################################
  
  
  ######################################################################################
  # Run Functions to CALCULATE OBSERVED ALLELE FREQUENCIES & RATE OF CHANGE IN LOG ODDS (refered to here as 'slope')
  
  # Run Allelotyper
  SNPsetFreq <- Allelotyper(SNPset)
  print(head(SNPsetFreq))  
  
  # Run Function LogisticRegLogisticRegScratch
  TempResults <-  lapply(1:nrow(SNPsetFreq), LogisticRegScratch, x=SNPsetFreq)
  SNPsetFreqSlope <- SNPsetFreq
  SNPsetFreqSlope$ObsSlope <- unlist(TempResults)
  SNPsetFreqSlope$AbsSlope <- abs(SNPsetFreqSlope$ObsSlope)
  head(SNPsetFreqSlope)
  
  # Rename columns to avoid overwriting later
  names(SNPsetFreqSlope) <- gsub("AbsSlope", "ObsAbsSlope", names(SNPsetFreqSlope))
  names(SNPsetFreqSlope) <- gsub("Ref", "ObsRef", names(SNPsetFreqSlope))
  names(SNPsetFreqSlope) <- gsub("Alt", "ObsAlt", names(SNPsetFreqSlope))
  print(head(SNPsetFreqSlope))
  
  # Create subset to send through BootStrapper (to avoid maxing out memory)                
  colsToKeep <- c("snpID", "^ObsRefAlleles", "^TotAlleles", "^DP", "^RD", "ObsAbsSlope", "Obsuared") 
  SNPsetFreqSlopeToBS <- SNPsetFreqSlope[,grep(paste(colsToKeep, collapse = "|"), colnames(SNPsetFreqSlope))]
  ######################################################################################
  
  
  ######################################################################################
  # START OF RESAMPLED DATA
  ######################################################################################
  # CREATE FUNCTION AlleleShuffler (rhyper)
  # samples without replacement to shuffle alleles among groups & converts to allele frequency
  
  AlleleShuffler <- function(x){  
    
    # Create sum of RefAlleles for all pools 
    x[,"TotObsRefAlleles"] <- x$ObsRefAlleles1 + x$ObsRefAlleles2 + x$ObsRefAlleles3
    
    # Create TotAlleles (ref and alt over all pools)
    x$TotAlleles <- x[1,"TotAlleles1"] + x[1,"TotAlleles2"] + x[1,"TotAlleles3"]
    
    # Create total Alt Alleles
    x$TotObsAltAlleles <- x$TotAlleles - x$TotObsRefAlleles
    
    # shuffle ref alleles across groups, i.e. the total number of ref alleles obs across all pools, and sample size each pool remains constant, but the distribution of ref alleles among the pools is randomized    
    # rhyper(nn, m, n, k)    # nn*k=m+n OR nn is number of elements in vector k and sum(k)=m+n
    # rhyper(#iterations sampled from same 'urn', #TotRef, #TotAlt, total alleles in a given group)
    
    # Draw References alleles for pool1
    x$shufRefAlleles1 <- t(t(apply(x[,c('TotObsRefAlleles', 'TotObsAltAlleles', 'TotAlleles1')], 1, function(i) rhyper(1, i['TotObsRefAlleles'], i['TotObsAltAlleles'], i['TotAlleles1']))))
    
    # Calculate Alternate alleles for pool1
    x$shufAltAlleles1 <- x$TotAlleles1 - x$shufRefAlleles1
    
    # Calculate number of Refs and Alts left in 'urn' after pool1 drawn
    x$RefLeft <- x$TotObsRefAlleles - x$shufRefAlleles1
    x$AltLeft <- x$TotObsAltAlleles - x$shufAltAlleles1
    
    # Draw Reference alleles for pool2, from those remaining in 'urn' after pool1 drawn
    x$shufRefAlleles2 <- t(t(apply(x[,c('RefLeft', 'AltLeft', 'TotAlleles2')], 1, function(i) rhyper(1, i['RefLeft'], i['AltLeft'], i['TotAlleles2']))))
    
    # Calculate Alternate alleles for pool2
    x$shufAltAlleles2 <- x$TotAlleles2 - x$shufRefAlleles2
    
    # Calculate Reference alleles for pool3
    x$shufRefAlleles3 <- x$TotObsRefAlleles - x$shufRefAlleles1 - x$shufRefAlleles2
    
    # Calculate Alternate alleles for pool3 (only necessary to verify totals below; could remove from final)
    x$shufAltAlleles3 <- x$TotObsAltAlleles - x$shufAltAlleles1 -  x$shufAltAlleles2
    
    # convert shuffledRefAls to shuffledRefAlFreq (shuffledRefAls/TotAlleles1)  
    x[, "RefFreq1"] <- x$shufRefAlleles1/x$TotAlleles1
    x[, "RefFreq2"] <- x$shufRefAlleles2/x$TotAlleles2
    x[, "RefFreq3"] <- x$shufRefAlleles3/x$TotAlleles3
    
    return(x)
    
  } # End of FUNCTION: AlleleShuffler
  ######################################################################################
  
  
  ######################################################################################
  # CREATE FUNCTION ResampReads
  # samples READS from tubes with replacement to model sequencer sampling error
  
  ResampReads <- function(x){    
    
    # Rename column names to avoid duplicates and confusion later
    names(x)[grep("^RD", colnames(x))] <- gsub("RD", "ObsRD", names(x)[grep("^RD", colnames(x))])
    names(x)[grep("^AD", colnames(x))] <- gsub("AD", "ObsAD", names(x)[grep("^AD", colnames(x))])
    
    # Create simulated Reference READS for all SNPs in each pool by choosing from a binomial distribution using the modified proportion of ref alleles observed as the probabilty of selecting that allele
    x$RD1 <- rbinom(nrow(x), x[,"DP1"], (x$shufRefAlleles1/x$TotAlleles1))
    x$RD2 <- rbinom(nrow(x), x[,"DP2"], (x$shufRefAlleles2/x$TotAlleles2))
    x$RD3 <- rbinom(nrow(x), x[,"DP3"], (x$shufRefAlleles3/x$TotAlleles3))
    
    x$AD1 <- x$DP1 - x$RD1
    x$AD2 <- x$DP2 - x$RD2
    x$AD3 <- x$DP3 - x$RD3
    
    return(x)
    
  } # End of FUNCTION: ResampReads
  ######################################################################################
  
  
  ######################################################################################
  # CREATE FUNCTION Pcounter
  # Implemented within BootStrapper
  # counts number of times observed slope is >= simulated slope
  # Note, here 'slope' refers to rate of change in log odds ratio
  
  Pcounter <- function(x){     
    
    # if simulated slope is greater than or equal to observed slope, add 1 to PcountSlope, if not keep it as is
    x$PcountSlope <- with(x, ifelse(AbsSlope >= ObsAbsSlope, PcountSlope + 1, PcountSlope))
    
    # if simulated  is greater than or equal to observed , add 1 to PcountSlope, if not keep it as is
    # x$Pcount <- with(x, ifelse(uared >= Obsuared, Pcount + 1, Pcount))
    
    return(x)
  } # End of FUNCTION: Pcounter
  ######################################################################################
  
  
  ######################################################################################
  # CREATE FUNCTION BootStrapper_Logistic to repeat pipepline of functions nBS times
  # models both population and sequencer sampling error
  # Functions: 
  #   AlleleShuffler 
  #   ResampReads
  #   Allelotyper
  #   LogisticRegScratch
  #   Pcounter   
  
  BootStrapper_Logistic <- function(x) {
    # if PcountSlope column exists (i.e., the first cycle was completed), replace them with the equavlent columns output at the end of the last cycle; 
    # if they do not exist (i.e., this is the first cycle), create them populated with 0s
    if(!"PcountSlope" %in% colnames(x)){x$PcountSlope <- 0} else{x$PcountSlope <- simData$PcountSlope}
    
    # function pipeline
    simData <- AlleleShuffler(x)
    simData <- ResampReads(simData)
    simData <- Allelotyper(simData)
    TempResults <-  lapply(1:nrow(simData), LogisticRegScratch, x=simData)
    simData$AbsSlope <- abs(unlist(TempResults))
    simData <- Pcounter(simData)
    return(simData)
  } # END of FUNCTION: BootStrapper_Logistic
  ######################################################################################
  
  
  
  ######################################################################################
  # CREATE FUNCTION BootStrapper_Logistic_noResampReads to repeat pipepline of functions nBS times
  # models population sampling error, but not sequencer sampling error
  # Functions: 
  #   AlleleShuffler 
  #   LogisticRegScratch
  #   Pcounter   
  
  BootStrapper_Logistic_noResampReads <- function(x) {
    
    # if PcountSlope and Pcount columns exisit (i.e., the first cycle was completed), replace them with the equavlent columns output at the end of the last cycle; 
    # if they do not exist (i.e., this is the first cycle), create them populated with 0s
    if(!"PcountSlope" %in% colnames(x)){x$PcountSlope <- 0} else{x$PcountSlope <- simData$PcountSlope}
    
    # function pipeline
    simData <- AlleleShuffler(x)
    TempResults <-  lapply(1:nrow(simData), LogisticRegScratch_noResampReads, x=simData)
    simData$AbsSlope <- abs(unlist(TempResults))
    simData <- Pcounter(simData)
    return(simData)
  } # END of FUNCTION BootStrapper_Logistic_noResampReads
  ######################################################################################
  
  
  ######################################################################################
  # Run LogRegAnalysis: BootStrapper, nBS times (1:nBS) and calculate P values
  
  library(parallel)
  
  # Identify total number of cores possible - parallel processing not supported on Windows    
  if(is.null(nCores)){
    if(Sys.info()['sysname']=='Windows'){
      nCores <- 1
    } else {
      nCores <- detectCores()-1
    }
  }
  
  # Select appropriate pipeline to run  
  if(ModelSeqSampErr==TRUE){
    print("ModelSeqSampError TRUE; Running BootStrapper_Logistic")
    parallelBSresults <- mclapply(1:nBS, function(i) BootStrapper_Logistic(SNPsetFreqSlopeToBS), mc.cores = nCores)
  } else {
    print("ModelSeqSampError FALSE; Running BootStrapper_Logistic_noResampReads")
    parallelBSresults <- mclapply(1:nBS, function(i) BootStrapper_Logistic_noResampReads(SNPsetFreqSlopeToBS), mc.cores = nCores)
  }
  
  # Gather results from BootStrapper and CALCULATE P VALUE
  # extract PcountSlope columns across all nBS runs
  AllPcountSlopes <- as.data.frame(sapply(parallelBSresults, function(x) x[["PcountSlope"]]))
  
  # create column of nBS count + 1 for obs  (one list created for each nBS)   
  SNPsetFreqSlope$nBSwithObs <- nrow(summary(parallelBSresults)) + 1
  
  # sum PcountSlope columns across all nBS runs plus 1 for ObsSlope
  SNPsetFreqSlope$SumPcountSlopes <- rowSums(AllPcountSlopes[,grep("nBSwithObs", colnames(AllPcountSlopes), invert = TRUE)]) + 1
  
  # Create columns with pvalues 
  SNPsetFreqSlope$SlopeP <- SNPsetFreqSlope$SumPcountSlopes/SNPsetFreqSlope$nBSwithObs 
  
  # Optional: SAVE results table as csv (does not include the results of each iteration, only Pvalues resulting from them)
  if(!is.null(OutputBase)){
    filename_SlopeCSV <- paste(OutputBase,"_", ResampYN, nBS, "bs_", "SNPsetFreqSlope.csv", sep = "")
    write.csv(SNPsetFreqSlope, file = filename_SlopeCSV)
    print(paste("FIN!", "___", OutputBase,"_", ResampYN, nBS, "bs_", sep = ""))
  } else {
    print(paste("FIN!", "___", ResampYN, nBS, "bs_", sep = ""))
  }
  # END of LogRegTest analysis
  ####################################################################################
  
  return(SNPsetFreqSlope)
}
