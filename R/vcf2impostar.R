#' Convert VCF to impostar format
#'
#' This function reformats a VCF file into the format required by impostar functions runLogRegTest and runAMOVA
#' @usage vcf2impostar(vcfFile="filename", designFile=NULL, ploidy=2, savecsv=FALSE)
#' @param vcfFile a character string indicating the vcf file you want to convert to impostar format. Required.
#' @param designFile a character string indicating the design file, which contains the number of indiviudals in each pool. Only required for pools. Default is NULL and should be left this way for individuals. The file can be .csv or .txt (tab-delimited)
#' @param ploidy a numeric indicating the ploidy of the focal organism. Default is 2. Note for pooled datasets: ploidy should be set to the ploidy of a single organism, the number of alleles included in a pool should be indicated in the designFile.
#' @param savecsv a logical indicating if you want to save the output as a csv file. If FALSE (default), the output will only be stored as an R object.
#' @return The vcf2impostar function returns a dataframe and optional csv file written to the working directory. The dataframe contains a row for each SNP locus and columns containing a unique name for each SNP based on its sequence location (snpID), reference allelic state (REF), alternate allelic state (ALT), reference reads (RD), alternate reads (AD), total read depth (DP), and total number of alleles (2n, TotAlleles) per individual or pool (as indicated by sequential numbers in the column names).
#' @author Rebecca M. Hamner
#' @references Hamner, R.M., J.D. Selwyn, E. Krell, S.A. King, and C.E. Bird. In review. Modeling next-generation sequencer sampling error in pooled population samples dramatically reduces false positives in genetic structure tests.
#' @seealso \code{\link{runLogRegTest}}, \code{\link{runAMOVA}}
#' @examples
#' # load example vcf file
#' myvcfFile <- impostar::: exVCF
#' # generate an example designFile
#' designFile <- data.frame(n=rep(20,6), Sample=c(seq(1:6)), Region=c(1,1,2,2,3,3))
#' # run vcf2impostar format converter
#' data.impostar <- vcf2impostar(vcfFile = vcfFile, designFile = designFile, ploidy = 2, savecsv = FALSE)
#' # to feed it directly into runAMOVA
#' testAMOVA <- runAMOVA(designFile = designFile, dataFile = data.impostar, NresamplesToStop = 5, maxPermutations = 10)
#' @import VariantAnnotation
#' @export
vcf2impostar <- function(vcfFile, designFile=NULL, ploidy=2, savecsv=FALSE) {
    # vcfFile: a character string indicating the vcf file you want to convert to impostar format. Required.
    # designFile: a character string indicating the design file, which contains the number of individuals in each pool. Only required for pools. Default is NULL and should be left this way for individuals.
    # ploidy: a numeric indicating the ploidy of the focal organism. Default is 2. Note for pooled datasets: ploidy should be set to the ploidy of a single organism, the number of alleles included in a pool should be indicated in the designFile.
    # savecsv: a logical indicating if you want to save the output as a csv file. If FALSE (default), the output will only be stored as an R object.

    # Read in vcf
    library(VariantAnnotation)

    if (class(vcfFile) == "character") {
      vcf <- readVcf(vcfFile)
    } else {
      vcf <- vcfFile
    }

    # if individuals and no design file, create necessary info
    if (is.null(designFile)){
      designFile <-rep(2, dim(vcf)[2])
    }

    # determine if design file is file name or R object, and read in as necessary
    if (class(designFile) != "data.frame") {
      fileExtention <- unlist(strsplit(designFile,"[.]"))[2]
        if (fileExtention == "csv"){
          designFile <- read.csv(designFile)
        } else {
          designFile <- read.table(designFile, header=TRUE)
        }
    }

    # extract necessary info from vcf
    snpID <- matrix(row.names(vcf), dimnames=list(seq(1:dim(vcf)[1]),"snpID"))
    # snpID <- matrix(seqnames(vcf), dimnames=list(row.names(vcf),"snpID"))
    REF <- as.vector(rowRanges(vcf)$REF)
    ALT <- as.vector(unlist(rowRanges(vcf)$ALT))

    # Extract Ref and Alt reads
    prefixRef <- "RD"
    prefixAlt <- "AD"
    prefixTotAlleles <- "TotAlleles"
    seqNums <- seq(1:dim(vcf)[2])

    dimnamesRef <- list()
    dimnamesRef[[1]] <- row.names(vcf) # rows
    dimnamesRef[[2]] <- paste0(prefixRef, seqNums) # columns

    dimnamesAlt <- list()
    dimnamesAlt[[1]] <- row.names(vcf) # rows
    dimnamesAlt[[2]] <- paste0(prefixAlt, seqNums) # columns

    RDall <- matrix(as.numeric(unlist(geno(vcf)$RO)), nrow=nrow(geno(vcf)$RO), dimnames = dimnamesRef)
    ADall <- matrix(as.numeric(unlist(geno(vcf)$AO)), nrow=nrow(geno(vcf)$AO), dimnames = dimnamesAlt)

    # calculate total depth(DP)
    DPall <- RDall + ADall
    colnames(DPall) <- gsub("RD", "DP", colnames(DPall))

    # extract total allles per pool or individual
    TotAlleles <- designFile[,1]*ploidy

    TotAllelesPerPool <- list()
    for (i in 1:dim(vcf)[2]){
      TotAllelesPerPool[[i]] <- rep(TotAlleles[i], nrow(vcf))
    }

    dimnamesTotAlleles <- list()
    dimnamesTotAlleles[[1]] <- row.names(vcf) # rows
    dimnamesTotAlleles[[2]] <- paste0(prefixTotAlleles, seqNums) # columns

    TotAllelesAll <- matrix(as.numeric(unlist(TotAllelesPerPool)), nrow = dim(vcf)[1], dimnames = dimnamesTotAlleles)

    # consolidate into a dataframe
    dataset <- as.data.frame(cbind(snpID, REF, ALT, RDall, ADall, DPall, TotAllelesAll))

    # update appropriate classes of columns
    for (i in 1:ncol(dataset)) {
      dataset[,i] <- as.character(dataset[,i])
    }

    for (i in 4:ncol(dataset)) {
      dataset[,i] <- as.numeric(dataset[,i])
    }

    # OPTION to write csv file
    if(savecsv==TRUE){
      fileBase <- unlist(strsplit(vcfFile,"[.]"))[1]
      mycsvFile <- paste0(fileBase, "_impostar.csv")
      write.csv(dataset, file = mycsvFile)
    }

    return(dataset)
}
