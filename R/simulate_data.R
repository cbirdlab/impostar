#' Generate Simulated Data 

#' 

#' @description This function is used to generate a simulated SNP dataset for any number of populations with specified numbers of individuals per population, reads per SNP per population, reference allele frequencies per population, and SNPs to be simulated. Intended for use with the runLogRegTest or runAMOVA functions.

#' @usage simulate_data(nIndiv, nReads, RefProb, nSNPs, file_name=F)

#' @param nIndiv a vector of integers specifying the number of individuals in each pool. The vector length indicates the number of pools and MUST be of equal length to nReads and RefProb. Required.

#' @param nReads a vector of integers specifying the number of reads per SNP in each pool. The vector length indicates the number of pools and MUST be of equal length to nIndiv and RefProb. Required.

#' @param RefProb a vector of numerics specifying the true reference allele frequency for all SNPs in each pool. The vector length indicates the number of pools and MUST be of equal length to nIndiv and nReads. Required.

#' @param nSNPs an integer indicating the number of SNPs to simulate. Required. 

#' @param file_name a logical or string indicating if the resulting simulated dataset should be written to a 'csv' file in the working directory. If TRUE, a file will be written and named based on the number of individuals, reads, and allele frequencies. If a string, the string will be appended with '.csv' and written to the working directory. If FALSE (default), no file will be written.

#' @return The simulate_data function returns a dataframe and optional csv file. Each row is a SNP and columns represent the total number of alleles (2n, TotAlleles), read depth (DP), reference alleles (RefAl), alternate alleles (AltAl), reference reads (RD), and alternate reads (AD) per population (as indicated by the number in the column names).

#' @author  Rebecca M. Hamner, Jason D. Selwyn, Evan Krell, Scott A. King, Christopher E. Bird

#' @references Hamner, R.M., J.D. Selwyn, E. Krell, S.A. King, and C.E. Bird. In review. Modeling next-generation sequencer sampling error in pooled population samples dramatically reduces false positives in genetic structure tests.

#' @seealso \code{\link{runLogRegTest}}, \code{\link{runAMOVA}}

#' @keywords simulated SNPs 
 
#' @examples 
#' # if arguments are read into variables
#' simdata <- simulate_data(nIndiv, nReads, RefProb, nSNPs, file_name=F)
#' # if three pools with equal individuals, reads and reference allele frequencies and with csv written
#' simdata <- simulate_data(rep(50, 3), rep(100, 3), rep(0.5, 3), 1000, file_name=T)
#' # if arguments supplied directly to function for four pools with specified csv file name.
#' simdata <- simulate_data(c(20,45,50,55), c(100,115,200,150), c(0.1,0.15,0.5,0.9), 1000, file_name="mySimulatedData")

#' @export 
simulate_data <- function(nIndiv, nReads, RefProb, nSNPs, file_name=F) {
  # nIndiv: a vector of integers specifying the number of individuals in each pool. The vector length indicates the number of pools and MUST be of equal length to nReads and RefProb.
  # nReads: a vector of integers specifying the number of reads per SNP in each pool. The vector length indicates the number of pools and MUST be of equal length to nIndiv and RefProb.
  # RefProb: a vector of numerics specifying the true reference allele frequency for all SNPs in each pool. The vector length indicates the number of pools and MUST be of equal length to nIndiv and nReads.
  # nSNPs: an integer indicating the number of SNPs to simulate
  # file_name: a logical or string indicating if the resulting simulated dataset should be written to a  'csv ' file in the working directory. If TRUE, a file will be written and named based on the number of individuals, reads, and allele frequencies. If a string, the string will be appended with '.csv' and written to the working directory. If FALSE (default), no file will be written.
  
  # Create columns required for subsequent analysis pipeline
  CHROM <- rep(capture.output(cat("simSNPprob", RefProb, sep = "_")), nSNPs)
  POS <- seq(1,nSNPs, 1)
  snpID <- paste(CHROM, POS, sep="_")
  
  # Create dataframe to hold sim SNPs  
  simSNPs <- as.data.frame(snpID)
  
  # Add total number of alleles per pool to dataframe 
  simSNPs[paste('TotAlleles',1:length(nIndiv),sep='')]<-2*matrix(rep(nIndiv,nSNPs),ncol=length(nIndiv),byrow = T)
  
  # Set total number of reads to simulate per pool
  simSNPs[paste('DP',1:length(nIndiv),sep='')]<-matrix(rep(nReads,nSNPs),ncol=length(nReads),byrow = T)
  
  # Randomly generate ref ALLELES for each pool with specified probability of ref/alt 
  # represents researcher sampling from population
  for(pop in 1:length(nIndiv)){
    simSNPs[paste('RefAl',pop,sep='')] <- rbinom(nSNPs, simSNPs[[paste('TotAlleles',pop,sep='')]], RefProb[pop])
  }
  
  # Calculate alt ALLELES
  simSNPs[paste('AltAl',1:length(nIndiv),sep='')] <- simSNPs[paste('TotAlleles',1:length(nIndiv),sep='')] - simSNPs[paste('RefAl',1:length(nIndiv),sep='')]
  
  
  # SAMPLE READS from ALLELES
  # randomly generate reference reads for each pool with probability of ref = proportion of ref obs in previous step
  for(pop in 1:length(nIndiv)){
    simSNPs[paste('RD',pop,sep='')] <- rbinom(nSNPs, simSNPs[[paste('DP',pop,sep='')]], simSNPs[[paste('RefAl',pop,sep='')]]/simSNPs[[paste('TotAlleles',pop,sep='')]])
  }
  
  # Calculate alt READS
  simSNPs$AD1 <- simSNPs$DP1 - simSNPs$RD1
  simSNPs$AD2 <- simSNPs$DP2 - simSNPs$RD2
  simSNPs$AD3 <- simSNPs$DP3 - simSNPs$RD3
  
  simSNPs[paste('AD',1:length(nIndiv),sep='')] <- simSNPs[paste('DP',1:length(nIndiv),sep='')] - simSNPs[paste('RD',1:length(nIndiv),sep='')]
  
  # Optional: write output to file
  if(!(file_name==F)){
    if(file_name==T){
      catnIndiv <- capture.output(cat(nIndiv, sep="_"))
      catnReads <- capture.output(cat(nReads, sep="."))
      catRefProb <- capture.output(cat(RefProb, sep="_"))
      filename <- paste0("simSNP10000n",catnIndiv, "reads", catnReads, "prob", catRefProb, ".csv")
      write.csv(simSNPs, file = filename)#, row.names = F)
    } else {
      write.csv(simSNPs, file = paste(file_name,'csv',sep='.'))#, row.names = F)
    }
  }
  
  simSNPs
}
