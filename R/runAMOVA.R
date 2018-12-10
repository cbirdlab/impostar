#' Run AMOVA
#'
#' This function performs an analysis of molecular variance (AMOVA) test for genetic structure, implementing the appropriate resampling strategies to model the total sampling error associated with either Sanger sequencing (population sampling error) or next-generation sequencing (population + sequencer sampling error). The current implementation works with pooled data, no missing data and one level of population structure.
#' @usage runAMOVA(designFile, dataFile, outputFile=NULL, NresamplesToStop=1000, maxPermutations=10000, multi.core = TRUE, do.bootstrap = FALSE, save.distributions = FALSE)
#' @param designFile 	a data.frame or character string indicating the path to the experimental design file. Required.
#' @param dataFile a data.frame or character string indicating the path to the data file. Required.
#' @param outputFile a character string indicating the name of the optional .rds output file written to working directory. If NULL (default), output is only retained as an R object. Note that writing large .rds files will increase time to completion.
#' @param NresamplesToStop an integer indicating the number of iterations to complete, after which resampling will stop if resampled F >= observed F. Default = 1000.
#' @param maxPermutations an integer indicating the maximum number of iterations to complete. Default = 10000
#' @param do.bootstrap a logical that determines whether sequencer sampling error (uncertainty in genotype calls) is modelled by bootstrapping the reads in each permutation in the permutation test.  not to be confused with the "rexact" shuffle method
#' @param multi.core a logical or integer indicating the number of cores to use in parallel processing. If FALSE will run on one core. If TRUE (default) will detect OS and use the number of available cores minus 1. If integer will run on specified number of cores. Note that parallel processing not supported by Windows.
#' @param NGSdata a logical indicating whether the resampling strategy for next-generation sequencing data (TRUE, Default) or Sanger sequencing data (FALSE) should be implemented.
#' @param save.distributions a logical indicating if the distributions of the F values during bootstrapping be returned. Note that this might require a very large amount of space.
#' @param ploidy an integer indicating the ploidy
#' @param permutationMethod a string that determines the method of shuffling.  "exact" is Fisher's exact permutation method and should be used in most cases.  "freely" will freely shuffle the smallest units of observation. "bird" is an experimental combination of exact and freely, "rexact" is an experimental bootstrapping procedure that should not be used 
#' @param preshuffle a logical that uses same shuffling in permutation tests for each SNP, such that permutation results can be combined across loci
#' @param multi.node a logical that must be set to FALSE (in devel)
#' @return The runAMOVA function returns a list and optional .rds file written to the working directory. Each list element contains the AMOVA table for one SNP, in the same order as the input file.
#' @author  Scott A. King, Christopher E. Bird, Rebecca M. Hamner, Jason D. Selwyn, Evan Krell
#' @references Hamner, R.M., J.D. Selwyn, E. Krell, S.A. King, and C.E. Bird. In review. Modeling next-generation sequencer sampling error in pooled population samples dramatically reduces false positives in genetic structure tests.
#' @seealso \code{\link{simulate_data}}, \code{\link{runLogRegTest}}
#' @examples
#' # create design file
#' design <- data.frame(n=rep(20,3), Sample=c(1,2,3))
#' # simulate data file
#' simdata <- simulate_data(rep(50, 3), rep(100, 3), rep(0.5, 3), 5, file_name=T)
#' # run AMOVA
#' AMOVAresults <- runAMOVA(designFile=design, dataFile=simdata, outputFile="AMOVAoutput", NresamplesToStop=10, maxPermutations=100, multi.core = T, do.bootstrap = T)
#' @import Rcpp VariantAnnotation parallel
#' @export
##
## SAK Amova - Version 1.05 - June 22 - 2017
runAMOVA <- function(
  designFile,                       # required user input
  dataFile,                         # required user input
  NresamplesToStop=1000,            # optional user input, this stops the resampling if resampled F >= observed F at the user specified value
  ploidy=2,                         # ploidy
  maxPermutations=100000,           # maximum number of permutation to make
  permutationMethod = "exact",      # "exact" "rexact", "freely", "bird", "alleles"
  multi.core = TRUE,                # use multiple core (but shared memory)
  preshuffle = FALSE,               # use a preshuffle, so all SNPs have the same shuffle
  do.bootstrap = FALSE,             # should the reads be bootstrapped during shuffling?
  save.distributions = FALSE,      # should the distributions of the fvals during shuffling be saved and returned?
  multi.node = FALSE                # use multiple nodes (message passing)
)
{
sfw.version <- " 1.04 - 3 June 2017"
permcntsToReportOn <- 10^(1:20)
  # force evaluation so won't have issues when sending to other nodes in a cluser!
  force(NresamplesToStop)
  force(maxPermutations)
  force(permutationMethod)
  force(ploidy)
  force(multi.core)
  force(multi.node)
  force(preshuffle)
  force(do.bootstrap)
  force(save.distributions)
  do.bootstrapReadsOld <- FALSE;  
  do.bootstrapReads <- FALSE;  
  do.bootstrapAlleles <- FALSE;  # First bootstrap method of bootstrapping alleles.
  do.bootstrapReads <- do.bootstrap;   # New method
  force(do.bootstrapAlleles)
  force(do.bootstrapReads)
  library(Rcpp)  

  cat( "NresamplesToStop:", NresamplesToStop, " maxPermutations:", maxPermutations, " permutationMethod:", permutationMethod, " multi.core:", multi.core, " preshuffle:", preshuffle,
  " do.bootstrap:", do.bootstrap, 
  " do.bootstrapReads:", do.bootstrapReads,
  " do.bootstrapAlleles:", do.bootstrapAlleles,
  " save.distributions:", save.distributions,
  " multi.node:", multi.node)
  

  # List of variables and their meanings
  #
  #  these variables don't change for each SNP
  #
  # source1.alleles.indices  - indices for the alleles, grouped by source1
  #   srcX.in.srcY[[1]][[2]]
  # source2.alleles.indices  - indices for the alleles, grouped by source2
  #   srcX.in.srcY[[1]][[3]]
  
  # source2.source1.indices          - indices for source1, grouped by source2 
  #   srcX.in.srcY[[2]][[3]]
  
  #srcX.in.srcY <- list(list(srcX.in.srcY[[1]][[2]], srcX.in.srcY[[1]][[2]]), list(NULL, srcX.in.srcY[[1]][[3]]))
  
  
  # source1.numalleles        - number of alleles in each source1
  # num.alleles.in.src[[2]]
  # source2.numalleleXs        - number of alleles in each source2
  # num.alleles.in.src[[3]]
  
  # source1.of.alleles       - index of the source1 for each allele
  #   srcX.in.srcY[[2]][[1]]
  # source2.of.alleles       - index of the source2 for each allele
  #   srcX.in.srcY[[3]][[1]]
  # source2.of.sourceX1       - index of the source2 for each source1
  #   srcX.in.srcY[[3]][[2]]
  
  
  # source1.names            - name of the source1s   - probably just an int id
  #  source.names[[2]]
  # source2.names            - name of the source2s   - probably just an int id
  #  source.names[[3]]
  # n.source1                - number of source1s
  #  num.source[2]
  # n.source2                - number of source2s
  #  num.source[3]
  # n.total.alleles          - number of alleles
  #  num.source[1]
  
  
  #  Should be per SNP, read into a matrix.
  #
  # alternate.reads.N        - alternate reads in the source1                Not used now
  # reference.reads.N        - reference reads in the source1                Not used now
  # total.reads              - Total number of reads (alternate+reference)   Not used now
  
  #
  # These variables change per SNP - and during shuffling
  #
  # alleles.ref.N            - number of reference alleles in each source1  
  # num.reference.alleles[[2]]
  # alleles.alt.N            - number of alternate alleles in each source1
  # num.alternate.alleles[[2]]
  # alleles.state            - state of each allele,  0 - ref,  1 - alternate
  
  # alts.source2             - number of alternate alleles in each source2
  # num.alternate.alleles[[3]]
  # 
  # source1.alt.frequencies  - frequencies of alternate alleles in each source1 (alts/(alts+refs))
  # alt.frequencies[[2]]
  # source2.alt.frequencies  - frequencies of alternate alleles in each source2 (alts/(alts+refs))
  # alt.frequencies[[3]]
  
  
  
  # If not defined in this scope, can't be seen.
  
  SSW <- vector();           # Sum of Squares - Within
  
  cat("SAK Amova - version ", sfw.version, "\n")
  
  
  if(class(dataFile)=='data.frame'){
    experimentalData<-dataFile
    fileExtension<-'none'
  } else {
    print(paste("processing", designFile, dataFile, "with permutation method", permutationMethod))
    fileExtension <- unlist(strsplit(dataFile,"[.]"))[2]
    
    if (fileExtension == "vcf"){
      library(VariantAnnotation)  # To read vcf files  
      print("Reading VCF file\n")
      vcf <- readVcf(dataFile, "foo")          #Not sure what the "foo" name is for. what is readVcf?
    } else {
      #print("Reading CSV Selwyn format\n")
      experimentalData <- read.csv(dataFile)
    }
  }
  
  if(class(designFile)=='data.frame'){
    expDesign<-designFile
  } else {
    expDesign <- read.csv(designFile) 
  }
  
  num.sources.of.variation <-dim(expDesign)[2]                      #calculate from experimental design file
  cat("There are ", num.sources.of.variation-1, " sources of variation: ", names(expDesign)[2:num.sources.of.variation], "\n")

## SAK took this check out for debugging, I suspect this is true though.
## if (num.sources.of.variation == 2) permutationMethod <- "freely"  # effectly only option
  
  num.alleles.in.src <- sapply(1:(num.sources.of.variation+1), function(x) list())  
  
  if (fileExtension == "vcf"){
    reference.reads.N <- matrix(as.numeric(unlist(geno(vcf)$RO)), nrow=nrow(geno(vcf)$RO))   
    alternate.reads.N <- matrix(as.numeric(unlist(geno(vcf)$AO)), nrow=nrow(geno(vcf)$AO))
    num.alleles.in.src[[2]] <- read.csv("King/individuals.per.pool.csv")[,3]*ploidy
  } else {
    if (length(grep("RD[0-9]+",colnames(experimentalData))) > 0) {
      reference.reads.N <-matrix(as.numeric(unlist(experimentalData[,grep("RD[0-9]+",colnames(experimentalData))])), nrow=dim(experimentalData)[1])
      alternate.reads.N <-matrix(as.numeric(unlist(experimentalData[,grep("AD[0-9]+",colnames(experimentalData))])), nrow=dim(experimentalData)[1])
    } else if (length(grep("RO[0-9]+",colnames(experimentalData))) > 0) {
      reference.reads.N <-matrix(as.numeric(unlist(experimentalData[,grep("RO[0-9]+",colnames(experimentalData))])), nrow=dim(experimentalData)[1])
      alternate.reads.N <-matrix(as.numeric(unlist(experimentalData[,grep("AO[0-9]+",colnames(experimentalData))])), nrow=dim(experimentalData)[1])
    } else {
      print("Can't find the allele read counts.  Expecting column header pairs of either RO#/AO# or RD#/AD#")
    }
    num.alleles.in.src[[2]] <- expDesign[,1]*ploidy
    numSNPs = dim(experimentalData)[1]
  }
  num.alleles.in.src[[1]] <- 1;   # These values are not actully correct but make other calculations easier.
  total.reads <- reference.reads.N + alternate.reads.N
  num.source <- vector()
  num.source[2] <- length(num.alleles.in.src[[2]])
  all.source1s <- 1:num.source[2]   # this vector used a lot so store it
  
  
  # index 1 is for alleles
  # index 2 is for src1/individuals with comes from the N column of the expDesign.
  
  # The expdesign file basically tells you where each src1 is stored, or the srcx index of each src 1:
  # To get the alleles, you would just repeat each of these by the ploidy,  but is that neeeded?
  num.source[1] = sum(num.alleles.in.src[[2]])  # Total number of individuals
  
  # The following ranges of indices are used many times, so instead of calculating them each time, we will store them.
  one.to.num.sources.of.variation <- 1:num.sources.of.variation      
  if (num.sources.of.variation < 3) {
      variationSourcesList <- NULL
  } else {
      variationSourcesList <- 2:(num.sources.of.variation-1)
  }

  if (num.sources.of.variation < 2) {
      variationSourcesPlusAllelesList <- NULL
  } else {
      variationSourcesPlusAllelesList <- 1:(num.sources.of.variation-1)
  }


  # The srcX.in.srcY list of lists  is used multiple ways.   When X is < Y is telling membership, when X is > Y it is telling you the list of members.
  # And since R starts at 1, the alleles have to be in here,  so when looking at the experiment, its source 1, (whether that be individuals or pools)  will be stored
  # with an index of 2.
  #
  srcX.in.srcY <- sapply(one.to.num.sources.of.variation, function(x) list()) 
  srcX.in.srcY[[2]][[1]] <- rep(1:length(num.alleles.in.src[[2]]), num.alleles.in.src[[2]])      # src1 in alleles
  if (num.sources.of.variation > 2)
    for (x in variationSourcesList) {
# RMH start addition - don't think this is necessary
      #browser()
       expDesign[, c(x + 1)] <- factor(expDesign[, c(x + 1)], levels = as.character(unique(expDesign[, c(x + 1)])))
# RMH end addition      
      srcX.in.srcY[[x+1]][[x]] <- as.numeric(unlist(unique(expDesign[,c(x, x+1)])[2]))
      num.source[[x]]=length(srcX.in.srcY[[x+1]][[x]])
      if (x > 2)
	for (k in 1:(x-2)) {
	  srcX.in.srcY[[x+1]][[x-k]] <- srcX.in.srcY[[x+1]][[x]][srcX.in.srcY[[x]][[x-k]]]
	} 
    }
  # num.source  list stores how many things are members of each source of variation.  Remember that index 1 is for alleles.
  num.source[[num.sources.of.variation]] <- length(unique(expDesign[,num.sources.of.variation]))
  
  if (num.sources.of.variation > 2)
    for (x in 3:(num.sources.of.variation)) {
      for (y in 2:(x-1)) {
	srcX.in.srcY[[y]][[x]] <- lapply(1:num.source[x], function(a) which(srcX.in.srcY[[x]][[y]] == a))
      }
      # to get the number of alleles at a level, just sum up the alleles for that group at the earlier level.
      num.alleles.in.src[[x]]<-vapply(srcX.in.srcY[[x-1]][[x]], function(a) sum(num.alleles.in.src[[x-1]][a]),0)  
      srcX.in.srcY[[x]][[1]] <- rep(1:length(num.alleles.in.src[[x]]), num.alleles.in.src[[x]])         # 
    }
  num.alleles.in.src[[num.sources.of.variation+1]]<- sum(num.alleles.in.src[[num.sources.of.variation]])  
  
  if (num.sources.of.variation > 1) 
    for (x in 2:num.sources.of.variation) {
      srcX.in.srcY[[1]][[x]] <- lapply(1:length(num.alleles.in.src[[x]]), function(a) which(srcX.in.srcY[[x]][[1]]==a))
    }
  srcX.in.srcY[[1]][[num.sources.of.variation+1]] <- list()
  srcX.in.srcY[[1]][[num.sources.of.variation+1]][[1]] <- 1:num.alleles.in.src[[num.sources.of.variation+1]][1]

  
  resetSG <- sapply(one.to.num.sources.of.variation, function(x) vector())   ## 
  resetN <- sapply(one.to.num.sources.of.variation, function(x) list(c(1)))  # list of vectors, and set the first element to 1 for all of them.    ## 
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#  i in 2:3
#  j in 2:2  or 2:3

#          --- 
#          \             num.alleles.in.src[[j]]^2
#sg[i][j] = >  --------------------------------------------------
#          /   num.alleles.in.src[[i+1]][srcX.in.srcY[[i+1]][[j]]]
#          ---
#          --- 
#          \             num.alleles.in.src[[2]]^2
#sg[2][2] = >  --------------------------------------------------
#          /   num.alleles.in.src[[3]][srcX.in.srcY[[3]][[2]]]                 800
#          ---
#          --- 
#          \             num.alleles.in.src[[2]]^2
#sg[3][2] = >  --------------------------------------------------
#          /   num.alleles.in.src[[4]][srcX.in.srcY[[4]][[2]]]                 800
#          ---
#          --- 
#          \             num.alleles.in.src[[3]]^2
#sg[3][3] = >  --------------------------------------------------
#          /   num.alleles.in.src[[4]][srcX.in.srcY[[4]][[3]]]                 16
#          ---


# Note if all of the num.alleles.in.src[[x]] are the same, only need this once!!!
# for symmetric experiments, this is the case, but maybe not with pools
  same.sizes <- all(sapply(num.alleles.in.src, function(x) min(x) == max(x))) 

# regardless  srcX.in.srcY  does not change.  
  num.source[num.sources.of.variation+1] <- 1
  DF <<- sapply(one.to.num.sources.of.variation, function (i) num.source[i] - num.source[i+1])                           
  DF.total <<- num.source[1] - 1                                   # I'm invariant

  calcSG <- function(srcNs) {
    sg <- resetSG
    for (i in variationSourcesList) {  # skip src1 - alleles/error
      for (j in 2:i) {  # skip src1 - alleles/error
        sg[[i]][j] <- sum( (srcNs[[j]]^2)/ srcNs[[i+1]][srcX.in.srcY[[i+1]][[j]]])
      }
    }
    for (j in 2:num.sources.of.variation) {  # skip src1 - alleles/error
      sg[[num.sources.of.variation]][j] <- sum((srcNs[[j]]^2) / num.source[1])
    }
    return(sg)
  }
  sg.noshuffle <- calcSG(num.alleles.in.src)
  tsg <- sg.noshuffle    ## init
  ##browser()

  calcN <- function(sg) {
    n <- resetN
    #for (i in 2:num.sources.of.variation) {
      #for (j in 2:i) {
        #if (i == j) n[[i]][j] <- (num.source[1] - sg[[i]][j])/DF[i]        ##  7%    #3
        #else        n[[i]][j] <- (sg[[i-1]][j] - sg[[i]][j])/DF[i]
      #}
    #}
    for (i in 2:num.sources.of.variation)
      n[[i]][i] <- (num.source[1] - sg[[i]][i])/DF[i]

    if (num.sources.of.variation > 2)  # SAK fix me
      for (i in 3:num.sources.of.variation) {
	for (j in 2:(i-1)) {                             
	  n[[i]][j] <- (sg[[i-1]][j] - sg[[i]][j])/DF[i]  
	}
      }
    return(n)
  }
  n.noshuffle <- calcN(sg.noshuffle)
  tn <- n.noshuffle    ## init
  if (same.sizes) {
      print("Experiment is symettric so optimizing")
      calcSG <- function(srcNs) {return(sg.noshuffle)}
      calcN <- function(sg) {return(n.noshuffle)}
  }

#---------------------------------------------------------------------------------------------------------
# amova 
#
#
# SSW - Sum of Squares  - within    == sum (Alts*Refs/(Alts+Refs))    So as these approch each other (50/50) it is maximized
# DF  - Degrees of freedom          ==
# ssa - Sum of Square   - among
# MS  - Mean Square = ssa/DF
# 
# Values are dependant upon
# SSW which is calculated outside since it doesn't change very often
# srcNs - number of alleles in each grouping 
# sg  - well, this is calculated outside and doesn't appear to be needed anymore, n depends on it, but cacluated outside again
# n   - cacluated outside and it is used to calcuate variance (sigma squared)
# DF, calculated outside and won't change

  #amova<-function(srcNs, sg, n, fval.only = FALSE)
  amova<-function(srcNs, n, fval.only = FALSE)
  {
    ssa <- vector()
    MS  <- vector()
    f <- vector()
    ssa <- SSW[one.to.num.sources.of.variation+1] - SSW[one.to.num.sources.of.variation] # all        # 3% #7
#cat(paste("ssa is"), ssa,"\n")
    MS <- ssa / DF                            # src2indicies     # 4% #6
#cat(paste("MS is"), MS,"\n")
    
    sig.sq <- vector()           # .12s
    sig.sq[1] <- MS[1]           # .08s

    for (i in 2:num.sources.of.variation) {
      #tsum <- MS[i]
      tsum <- 0

      for (j in 2:i) {
        tsum <- tsum+(n[[i]][j-1]*sig.sq[j-1])
      }
      sig.sq[i] <- (MS[i] - tsum) / n[[i]][i]
    }
    sig.sq[num.sources.of.variation+1] <- sum(sig.sq[one.to.num.sources.of.variation])
    #browser()    
    
    #numerator <- sig.sq[1]

    # for (i in 2:num.sources.of.variation) {
    #   numerator <- numerator + sig.sq[i]
    #   f[i] <- sig.sq[i]/numerator
    # }
    # f[1] <- (numerator - sig.sq[1])/numerator
    #f[1] <- (numerator - sig.sq[1])/numerator 
    
    numerator<-0
    denominator <- sig.sq[1]
    for (i in 2:num.sources.of.variation) {
      numerator <- numerator + sig.sq[i]
      denominator <- denominator + sig.sq[i]
      if (denominator != 0) {
        f[i] <- sig.sq[i]/denominator
      } else {
        f[i] <- 0
      }
    }
    
    if(num.sources.of.variation >= 3){
      if (denominator != 0) {
        f[1] <- numerator/denominator 
      } else {
        f[1] <- 0
      }
    }
    
    
    
    if (fval.only) {
        return(f)
    } else {
    result = matrix(NA, num.sources.of.variation+1, 7)
    #result[1, ] <- c(DF[3], ssa[3],    MS[3],     sig.sq[3],  f[3],  NA, NA)
    #result[2, ] <- c(DF[2], ssa[2],    MS[2],     sig.sq[2],  f[2],  NA, NA)
    #result[3, ] <- c(DF[1], ssa[1],    MS[1],       sig.sq[1],    f[1], NA, NA)
    #result[4, ] <- c(DF.total, SSW[num.sources.of.variation+1], MS.total, sig.sq[num.sources.of.variation+1], NA, NA, NA);
    #   colnames(result) <- c("df","SS","MS","sigma squared","f","p","#permutations")
    #browser()
    for (i in num.sources.of.variation:1) {
      result[num.sources.of.variation-i+1, ] <- c(DF[i], ssa[i],    MS[i],     sig.sq[i],  f[i],  NA, NA)
      #  4-4+1 = 1 --s4     4-3+1 = 2 --s3   4-2+1 = 3 --s2    4-1+1 = 4 --s1  (error)
    }
    result[num.sources.of.variation+1, ] <- c(DF.total, SSW[num.sources.of.variation+1], MS.total, sig.sq[num.sources.of.variation+1], NA, NA, NA);
    
    #   rownames(result) <- c("source2", "source1", "error", "total")
    

    return(result)
    }

    ## Instead of return the matrix,  should just return the array..
    #return(c((DF, ssa,    MS,     sig.sq,  f,  NA, NA)   
  } 
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #AMOVA FUNCTION END
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#---------------------------------------------------------------------------------------------------------
  

##  Here are the different ways to do shuffling.  We set a variable to the correct one to call, so the rest of the code 
#   doesn't care which one is being called.  They all return either the number of alts in each entitity being
#   shuffled, or the indices of the entities being shuffled.

  # We shuffle the alleles via alleles.state
  # Then add each pair together, and we store the odd and even indices for a speed increase
  # This version optimized for individuals, that is 2 src 0 per src 1.
  shuffleSource1FreelyIndividuals <- function() {
    m<-sample(alleles.state)
    return(m[oddi] + m[eveni])   # speed hack for individuals!
  }

  shuffleSource1FreelyIndividualsPre <- function(p) {
    m <- alleles.state[preshuf[[1]][[p]]]
    return(m[oddi] + m[eveni])   # speed hack for individuals!
  }
  
  # shuffle the alleles by resample the alleles.state variable
  # then sum up alternate alleles by counting 1s in each src 1 grouping.
  # This version is when we have an unknown number of src 0s per src 1
  shuffleSource1FreelyAll <- function() {
    m<-sample(alleles.state)
    return(vapply(srcX.in.srcY[[1]][[2]], function(x) sum(m[x]),0))
  }

  shuffleSource1FreelyAllPre <- function(p) {
    m<-alleles.state[preshuf[[1]][[p]]]
    return(vapply(srcX.in.srcY[[1]][[2]], function(x) sum(m[x]),0))
  }
  
  #
  # If we know the number of source 0's that we are going to add for each src 1, and they are all the same
  # and in the case of individuals it is two, that means every pair get added together so lets store the indicies
  # for these pairs so we don't have to determine more than once.
  data.is.individuals <- FALSE
  if (max(num.alleles.in.src[[2]]) == 2 && min(num.alleles.in.src[[2]]) == 2) { # we have individuals
    print("Data is individuals so optimizing for that case.")
    data.is.individuals <- TRUE
    oddi <- seq(1,num.alleles.in.src[[num.sources.of.variation+1]], 2)
    eveni <- seq(2,num.alleles.in.src[[num.sources.of.variation+1]], 2)
#    shuffleSource1Freely <- shuffleSource1FreelyIndividuals
#    if (preshuffle) {
#      print("i1 Using shuffleSource1FreelyIndividualsPre")
#      shuffleSource1Freely <- function() {shuffleSource1FreelyIndividualsPre(curPermutationi) }
#    }
   }
#    else {
#    shuffleSource1Freely <- shuffleSource1FreelyAll
#  }
 


  #shuffleSource1.6b <- function()
  # this was version 1.6b, one attempt to be as fast as possible.  It is still unclear which is the best
  # version, althought the optimization for individuals is hard to beat.
  # The next one is about the same speed, but simpler, so we'll use it.
  shuffleSource1old <- function() {
    c <- 0;
    refs<-vector(length=num.source[2])
    for(i in 1:num.source[3]) {
      m<-sample(srcX.in.srcY[[1]][[3]][[i]], replace=FALSE)          
      for (s1 in srcX.in.srcY[[2]][[3]][[i]]) {
        refs[s1] <- sum(alleles.state[m[srcX.in.srcY[[1]][[2]][[s1]]-c]])
      }
      c <- c + length(m)
    }
    return(refs)
  }


# cppFunction('
#   NumericVector cshuf2(List l, NumericVector state, NumericVector Shuf) {
#     int n = l.size();
#     NumericVector results(n);
#     
#     for(int i=0; i<n ; ++i) {
#       NumericVector y(l[i]);
#       int m=y.size();
#       results[i] = 0;
#       for (int j=0; j<m; j++) {
#         results[i] += Shuf[y[j]-1];
#       }
#     }
#     return results;
#   }
# ')

# Shuffle source1 basic, but sum using c++
#  shuffled.alleles <- unlist(sapply(1:num.source[3], function(i) sample(alleles.state[srcX.in.srcY[[1]][[3]][[i]]],replace=FALSE)))
shuffleSource1withc <- function(p) {
  return(cshuf2(srcX.in.srcY[[1]][[2]], alleles.state, sample(alleles.state)))
}

# Shuffle source1 basic
  shuffleSource1 <- function(p) {
    shuffled.alleles <- unlist(sapply(1:num.source[3], function(i) sample(alleles.state[srcX.in.srcY[[1]][[3]][[i]]],replace=FALSE)))
    return(vapply(srcX.in.srcY[[1]][[2]], function(x) sum(shuffled.alleles[x]),0))
  }

# cppFunction('
# NumericVector cshuf(List l, NumericVector state, NumericVector pre) {
#   int n = l.size();
#   NumericVector results(n);
#   
#   for(int i=0; i<n ; ++i) {
#     NumericVector y(l[i]);
#     int m=y.size();
#     results[i] = 0;
#     for (int j=0; j<m; j++) {
#       results[i] += state[pre[y[j]-1]-1];
#     }
#   }
#   return results;
# }
#             ')

  shuffleSource1Pre <- function(p) {                      
    return(cshuf(srcX.in.srcY[[1]][[2]], alleles.state, preshuf[[1]][[p]]))
  }

  shuffleSource1Preold <- function(p) {
    #shuffled.alleles <- preshuf[[1]][[p]]
    shuffled.alleles <- alleles.state[preshuf[[1]][[p]]]
    return(vapply(srcX.in.srcY[[1]][[2]], function(x) sum(shuffled.alleles[x]),0))
  }

  # Same as version 1.6b, but not just for src 1.
  #
  # This will shuffle the alleles for src i, by constraining withing src i+1, and return the alternate counts
  # for src 1.
  # but it isn't doing that yet!

  preshuf <- vector(mode="list", length= (num.sources.of.variation-1))
  for (i in variationSourcesPlusAllelesList) {
     preshuf[[i]] <- vector(mode="list", length=maxPermutations)
  }

  preshuffleSourcei.alleles <- function(i, N) {
    for (j in 1:N) {
      preshuf[[i]][[j]] <<- unlist(lapply(srcX.in.srcY[[1]][[i+2]], function(x) sample(x)))
    }
  }

  # will only call for source 1
  preshuffleSourcei.freely <- function(i, N) {
    for (j in 1:N) {
      #preshuf[[i]][[j]] <<- sample(alleles.state)
      preshuf[[i]][[j]] <<- sample(1:num.source[i], replace=FALSE)
    }
  }

  preshuffleSourcei <- function(i, N) {

    for (j in 1:N) {
       #preshuf[[i]][[j]] <<- unlist(lapply(srcX.in.srcY[[1]][[i+2]], function(x) sample(alleles.state[x])))
       #preshuf[[i]][[j]] <<- unlist(lapply(srcX.in.srcY[[1]][[i+2]], function(x) sample(x)))
       if (i == (num.sources.of.variation-1))
	 preshuf[[i]][[j]] <<- sample(1:num.source[i], replace=FALSE)
       else
	 preshuf[[i]][[j]] <<- unlist(lapply(srcX.in.srcY[[i]][[i+2]], function(x) sample(x)))
    }
           #return(shuffle <- unlist(lapply(srcX.in.srcY[[level-1]][[level+1]], sample)))
  }


  getSourceiShuffleNoPre <- function(src) {
    if (src == num.sources.of.variation-1)                        
      return(sample(1:num.source[src], replace=FALSE))
    else
      return(shuffle <- unlist(lapply(srcX.in.srcY[[src]][[src+2]], sample)))
  }

  getSourceiShufflePre <- function(i,p) {
    return(preshuf[[i]][[p]])
  }

  if (preshuffle) {
    print("g1 Using getSourceiShufflePre")
    getSourceiShuffle <- function(i) {getSourceiShufflePre(i, curPermutationi) }
  } else {
    getSourceiShuffle <- getSourceiShuffleNoPre
  }


  shuffleSourceiAllelesIndividualsPre <- function(i,p) {
    #shuf <- preshuf[[i]][[p]]
    shuf <- alleles.state[preshuf[[i]][[p]]]
    return(shuf[oddi]+shuf[eveni])
  }

  shuffleSourceiAllelesIndividuals <- function(i) {
    shuf <- unlist(lapply(srcX.in.srcY[[1]][[i+2]], function(x) sample(alleles.state[x])))
    return(shuf[oddi]+shuf[eveni])
  }

  shuffleSourceiAllelesPre <- function(i,p) {                      
    #shuf <- preshuf[[i]][[p]]
    shuf <- alleles.state[preshuf[[i]][[p]]]
    #return(vapply(1:num.source[[2]], function(x) sum(shuf[srcX.in.srcY[[1]][[2]][[x]]]),0)) # 
    return(vapply(srcX.in.srcY[[1]][[2]], function(x) sum(shuf[x]),0))
  }

  shuffleSourceiAlleles <- function(i) {                      
    shuf <- unlist(lapply(srcX.in.srcY[[1]][[i+2]], function(x) sample(alleles.state[x])))  # 
    #return(vapply(fred, function(x) sum(shuf[x]),0)) # 
    return(vapply(srcX.in.srcY[[1]][[2]], function(x) sum(shuf[x]),0)) # 
    #return(vapply(1:num.source[[2]], function(x) sum(shuf[srcX.in.srcY[[1]][[2]][[x]]]),0)) # 
    #return(sapply(1:num.source[[2]], function(x) sum(shuf[srcX.in.srcY[[1]][[2]][[x]]])))
  }

  # when we have individuals, this is very fast.
  shuffleSource1Individuals <- function() {   
    shuf <- unlist(sapply(srcX.in.srcY[[1]][[3]], function(x) sample(alleles.state[x])))              # 
    return(shuf[oddi]+shuf[eveni])
  }

  shuffleSource1IndividualsPre <- function(p) {   
    #shuf <- preshuf[[1]][[p]]
    shuf <- alleles.state[preshuf[[1]][[p]]]
    return(shuf[oddi]+shuf[eveni])
  }

  num.reference.alleles <- sapply(one.to.num.sources.of.variation, function(x) list())  
  num.alternate.alleles <- sapply(one.to.num.sources.of.variation, function(x) list()) 
  alt.frequencies <- sapply(one.to.num.sources.of.variation, function(x) list())  

    shuffleSource1wreplacement <- function() {
      rbinom(num.source[2], num.alleles.in.src[[2]], alt.frequencies[[2]])
    } 
    

## Lots of different ways to sum up the number of alt/ref allels at each level.  Each one is efficient for different situations.
#  This section of code will try to use the most effecient method.  This can have a huge impact on running time.
#

##browser()
    if (permutationMethod == "rexact") {
      shuffleSource1Function <- shuffleSource1wreplacement
    } else if (permutationMethod == "freely") {
      if (max(num.alleles.in.src[[2]]) == 2 && min(num.alleles.in.src[[2]]) == 2) { # we have individuals
         if (preshuffle) {
	   print("Using f1 shuffleSource1FreelyIndividualsPre")
	   shuffleSource1Function <- function() {shuffleSource1FreelyIndividualsPre(curPermutationi) }
         } else {
	   print("Using f2 shuffleSource1FreelyIndividuals")
           shuffleSource1Function <- shuffleSource1FreelyIndividuals
	 }
      } else {
         if (preshuffle) {
	   print("Using f3 shuffleSource1FreelyAllPre")
	   shuffleSource1Function <- function() {shuffleSource1FreelyAllPre(curPermutationi) }
         } else {
	   print("Using f4 shuffleSource1FreelyAll")
	   shuffleSource1Function <- shuffleSource1FreelyAll
	 }
      }
    }
    else if (permutationMethod == "alleles") {
      if (max(num.alleles.in.src[[2]]) == 2) { # we have individuals
	print("a1 using shuffleSource1/iIndividuals");
        shuffleSource1Function <- shuffleSource1Individuals
        shuffleSourceiFunction <- shuffleSourceiAllelesIndividuals
	if (preshuffle) {
          print("a2 Using shuffleSource1IndividualsPre() and shuffleSourceiAllelesIndividualsPre()")
	  shuffleSource1Function <- function() {shuffleSource1IndividualsPre(curPermutationi) }
	  shuffleSourceiFunction <- function(i) {shuffleSourceiAllelesIndividualsPre(i, curPermutationi) }
	}
      } else  {
        print("a3 using shuffleSource1/i")
        shuffleSource1Function <- shuffleSource1
        shuffleSourceiFunction <- shuffleSourceiAlleles
	if (preshuffle) {
	  print("a4 using shuffleSource1Pre() and shuffleSourceiAllelesPre() ")
	  shuffleSource1Function <- function() {shuffleSource1Pre(curPermutationi) }
	  shuffleSourceiFunction <- function(i) {shuffleSourceiAllelesPre(i, curPermutationi) }
	}
      }
    } else { # "exact" and "bird" 
      if (max(num.alleles.in.src[[2]]) == 2) { # we have individuals
	print("e1 using shuffleSource1/iIndividuals");
        shuffleSource1Function <- shuffleSource1Individuals
        #shuffleSourceiFunction <- shuffleSourceiIndividuals
	if (preshuffle) {
          print("e2 Using shuffleSource1IndividualsPre() and shuffleSourceiIndividualsPre()")
	  shuffleSource1Function <- function() {shuffleSource1IndividualsPre(curPermutationi) }
	  #shuffleSourceiFunction <- function(i) {shuffleSourceiIndividualsPre(i, curPermutationi) }
	}
      } else  {
        print("e3 using shuffleSource1/i")
        shuffleSource1Function <- shuffleSource1
        #shuffleSource1Function <- shuffleSource1withc
        #shuffleSourceiFunction <- shuffleSourcei
	if (preshuffle) {
	  print("e4 using shuffleSource1Pre() and shuffleSourceiPre() ")
	  shuffleSource1Function <- function() {shuffleSource1Pre(curPermutationi) }
	  #shuffleSourceiFunction <- function(i) {shuffleSourceiPre(i, curPermutationi) }
	}
      }
    }
    #num.source[num.sources.of.variation+1] <- 1
num.alternate.alleles[[1]] <- 0;
num.reference.alleles[[1]] <- 1;
#num.alleles.in.src[[1]] <- 1;   # these are not correct, but will give us the correct values here.

    MS.total <<- SSW[num.sources.of.variation+1] / DF.total                                   # I'm invariant

  # if (max(num.alleles.in.src[[2]]) == 2 && min(num.alleles.in.src[[2]]) == 2) { # we have individuals

#  getCutoffs <- function(numReads) {
#    cutoffs <- rep(0,numReads+1)
#    heterorange <- round(.2*numReads):round(.8*numReads)
#    cutoffs[(round(.8*numReads)+1):numReads] <- 2
#    cutoffs[heterorange+1] <- 1 
#    cutoffs[1] <- 0    # At least one homozygote
#    cutoffs[numReads+1] <- 2    # At least one homozygote
#    return(cutoffs)
#  }
  getNumrefs <- function(r, N) {
    c1 <- round(.2*N)   # This could be stored?
    if (r <c1 || r == 0) return(0)  # Homo if no reads or < the 20% cutoff
    c2 <- round(.8*N)
    if (r >c2 || r == N) return(0)  # Homo if all reads or > the 80% cutoff
    return(1)    # otherwise it is a heterozygote
  }

  getNumrefs2 <- function(r, N) {
    results <- allhetero
    results[r < round(.2*N)] <- 0
    results[ r == 0] <- 0
    results[r > round(.8*N)] <- 2
    results[r == N] <- 2
    ##browser()
    return(results) 
  }
                            
  processSNP <- function(row, maxPermutations) {
    #if (row%%100 == 0) print(paste("processing", row, "/", numSNPs))
    print(paste("processing", row, "/", numSNPs))
    ##browser()
    
    if (save.distributions) {
      NullDistributionsTmp <- lapply(variationSourcesPlusAllelesList, function(x) vector(length = maxPermutations)) 
      NullDistributions <- list()
    }

    # Used for bootstrapping of reads
    # these only need to be calculated once so do them here.

    alt.reads <- sapply(alternate.reads.N[row,], function(x) if (x==0) x+1 else x)   # Get the vector of alternate reads for this SNP, no zeros
    ref.reads <- sapply(reference.reads.N[row,], function(x) if (x==0) x+1 else x)   # Get the vector of reference reads for this SNP, no zeros
    reads <- alt.reads+ref.reads                                                     # Get new total reads to calcuate percentages
  ##ReadsForSNP <- total.reads[row,]
    ReadsForSNP <- reads

    tube.ref.frequencies <- ref.reads/reads  

    
    # per SNP
    #ps <- reference.reads.N[row,]/total.reads[row,]
    if (data.is.individuals) {
      allhetero <<- rep(1, num.source[2])
      #reference.allleles.in.tube <- getNumrefs2(new.ref.reads,reads)
      #num.reference.alleles[[2]] <<- sapply(all.source1s, function(x) getNumrefs2(reference.reads.N[row,x], total.reads[row,x]))
      num.reference.alleles[[2]] <<- getNumrefs2(reference.reads.N[row,], total.reads[row,])
    } else {
      num.reference.alleles[[2]] <<- round((reference.reads.N[row,]/total.reads[row,])*as.vector(num.alleles.in.src[[2]]))
    }

    num.alternate.alleles[[2]] <<- num.alleles.in.src[[2]]-num.reference.alleles[[2]]
    alleles.state <<- rep(rep(c(0,1), length(num.alleles.in.src[[2]])), as.vector(rbind(num.reference.alleles[[2]], num.alternate.alleles[[2]])))   

#cat("Alleles state for SNP ", row)
#print(alleles.state)
    #distance.matrix <- abs(outer(alleles.state, alleles.state, FUN="-"))
    if (preshuffle && row == 1) {
      print("Setting up preshuffle")
      if (permutationMethod == "freely") {
	preshuffleSourcei.freely(1, maxPermutations)         # only need source 1 preshuffled
      } else {
	#for (i in variationSourcesPlusAllelesList ) 
	for (i in 1:(num.sources.of.variation-1) ) {
	  print(paste("preshuffle for src",i))
	  if (permutationMethod == "alleles")   
	    preshuffleSourcei.alleles(i, maxPermutations) # always shuffle alleles
	  else 
	    preshuffleSourcei(i, maxPermutations)         # shuffles srci-1 within src i
	}
      }
    }

    #num.alleles.in.src[[num.sources.of.variation+1]] <- length(alleles.state)  
    
    ## Look into this for rexact
    if (num.sources.of.variation > 2) {
      # These are used by shuffleSource1wreplacement()   used by the rexact method only
      alt.frequencies[[3]] <<- sapply(1:num.source[3], function(x) sum(num.alternate.alleles[[2]][srcX.in.srcY[[3]][[2]]==x])/
                                     (sum(num.alternate.alleles[[2]][srcX.in.srcY[[3]][[2]]==x])+sum(num.reference.alleles[[2]][srcX.in.srcY[[3]][[2]]==x])))
      alt.frequencies[[2]] <<- alt.frequencies[[3]][srcX.in.srcY[[3]][[2]]] # This is really just the source2 frequencies repeated
    }

# this one is slower, but the number of alleles often doesn't change, so we only count of the alternate 
# alleles, and then call this one instead of either summing up or subtracting to get ref counts.
calculateSSW <- function(levels) {
    vapply(levels, function(i) sum(num.alternate.alleles[[i]]*(num.alleles.in.src[[i]]-num.alternate.alleles[[i]])/num.alleles.in.src[[i]]),0)
}

# RMH: this is new since first submission, might be where multi.core is choking??
# shared library ccalcSSW is only available on process that compiled the function (so n-1 threads will fail)
# see: https://stackoverflow.com/questions/38518387/using-rcpp-inside-parlapply-within-the-parallel-r-package
# cppFunction('
#   NumericVector ccalcSSW(List Alts, List Ns, NumericVector levels) {
#     int n = levels.size();
#     NumericVector results(n);
#     
#     for(int i=0; i<n ; ++i) {
#       int l = levels[i]-1;   // C starts at 0, R at 1
#       NumericVector A(Alts[i]);
#       NumericVector N(Ns[i]);
#       int m=A.size();
#       results[l] = 0;
#       for (int j=0; j<m; j++) {
#         results[l] += A[j]*(N[j]-A[j])/N[j];
#       }
#     }
#     return results;
#   }
# ')

#print("num.alternate.alleles"); print(num.alternate.alleles)
#print("num.alleles.in.src"); print(num.alleles.in.src)

# Call this after changing alternate counts, so can't call this unless have recalculate ref counts as well.
calculateSSW2 <- function(levels) {
  vapply(levels, function(i) sum(num.alternate.alleles[[i]]*num.reference.alleles[[i]]/num.alleles.in.src[[i]]),0) 
}

#
# I'm broken.  So should either set the global num.alternate.alleles in here
# or return a list the correct size to set outside of the function.
# But then a little tricky when one list versus multiple lists
#
calculateAlternateCounts <- function(levels) {
  ac <- vector(mode="list", length=length(levels))
  i <- levels[1]-1
  ac[[i]]<-num.alternate.alleles[[i]]
  for (i in levels) {
    ac[[i]] <- vapply(srcX.in.srcY[[i-1]][[i]], function(x) sum(ac[[i-1]][x]),0)  #3% #10
  }
    #ac[[i]] <- vapply(srcX.in.srcY[[i-1]][[i]], function(x) sum(num.alternate.alleles[[i-1]][x]),0)  #3% #10
##browser()
  return(ac[levels])
}


##    SSW[1] <<- 0   # for alleles no SSW
#    for (i in 2:num.sources.of.variation) {
#      if (i > 2) {
#        num.alternate.alleles[[i]]<-vapply(srcX.in.srcY[[i-1]][[i]], function(x) sum(num.alternate.alleles[[i-1]][x]),0)
#        num.reference.alleles[[i]]<-vapply(srcX.in.srcY[[i-1]][[i]], function(x) sum(num.reference.alleles[[i-1]][x]),0)
#      }
##      SSW[i]<<-sum(num.alternate.alleles[[i]]*(num.alleles.in.src[[i]]-num.alternate.alleles[[i]])/num.alleles.in.src[[i]]) 
#    }

calculateAlternateAlleles <- function() {
  if (num.sources.of.variation > 2)
    for (i in 3:num.sources.of.variation) {
      num.alternate.alleles[[i]]<<-vapply(srcX.in.srcY[[i-1]][[i]], function(x) sum(num.alternate.alleles[[i-1]][x]),0)
    }
  num.alternate.alleles[[num.sources.of.variation+1]]<<- sum(num.alternate.alleles[[num.sources.of.variation]])
}

calculateReferenceAlleles <- function() {
  if (num.sources.of.variation > 2)
    for (i in 3:num.sources.of.variation) {
      num.reference.alleles[[i]]<<-vapply(srcX.in.srcY[[i-1]][[i]], function(x) sum(num.reference.alleles[[i-1]][x]),0)
    }
  num.reference.alleles[[num.sources.of.variation+1]]<<- sum(num.reference.alleles[[num.sources.of.variation]])
}

##browser()
    #num.alternate.alleles <- calculateAlternateAlleles()
    #num.reference.alleles <- calculateReferenceAlleles()
    calculateAlternateAlleles()
    calculateReferenceAlleles()
    SSW <<- calculateSSW(1:(num.sources.of.variation+1))
#-print("SSW is")
#-print(SSW)

## This function bootstraps the alleles
bootstrapAlleles <- function(alts) {
  freq = 1-alts/num.alleles.in.src[[2]]
  ref.reads<-rbinom(num.source[2], total.reads[1,], freq)
  num.reference.alleles <- round((ref.reads/total.reads[1,])*as.vector(num.alleles.in.src[[2]]))
  num.alternate.alleles <- num.alleles.in.src[[2]]-num.reference.alleles
  #alleles.state <<- rep(rep(c(0,1), length(num.alleles.in.src[[2]])), as.vector(rbind(num.reference.alleles[[2]], num.alternate.alleles[[2]])))   
  return(num.alternate.alleles)
}
## This function bootstraps the alleles
## alts is the new number of alts in each source 1
bootstrapReads <- function(alts) {
  freq = 1-alts/num.alleles.in.src[[2]]
  #ref.reads<-rbinom(num.source[2], total.reads[1,], freq)
  ref.reads<-rbinom(num.source[2], ReadsForSNP, freq)
  num.reference.alleles <- round((ref.reads/ReadsForSNP)*as.vector(num.alleles.in.src[[2]]))
  num.alternate.alleles <- num.alleles.in.src[[2]]-num.reference.alleles
  return(num.alternate.alleles)
}

# This function will bootstrap the number of reads.
# 
bootstrapReadsOld <- function(alts) { 
  # Be wary, these really shouldn't be calculated every time since they won't change

  # Now resample the tube using the original number of reads using the percentage of the original reads 

  #new.ref.reads <- rbinom(num.source[2], total.reads[row,], tube.ref.frequencies)
  #new.alt.reads <- total.reads[row,]-new.ref.reads
  ref.reads <- sapply(reference.reads.N[row,], function(x) if (x==0) x+1 else x)   # Get the vector of reference reads for this SNP, no zeros
  reads <- alt.reads+ref.reads                                                     # Get new total reads to calcuate percentages
  tube.ref.frequencies <- ref.reads/reads  

  new.ref.reads <- rbinom(num.source[2], reads, tube.ref.frequencies)
  new.alt.reads <- reads-new.ref.reads

  if (data.is.individuals) 
    #reference.allleles.in.tube <- sapply(all.source1s, function(x) getNumrefs(new.ref.reads[x],total.reads[row, x]))
    reference.allleles.in.tube <- getNumrefs2(new.ref.reads,reads)
  else
    reference.allleles.in.tube <- round((new.ref.reads/total.reads[row, ])*as.vector(num.alleles.in.src[[2]]))

##browser()
  if (identical(reference.allleles.in.tube, num.reference.alleles[[2]])) {
##browser()
    return()  # do nothing else
  } else { # number of alleles has changed so recalculate everything!
    num.reference.alleles[[2]] <<- round((ref.reads/reads)*as.vector(num.alleles.in.src[[2]]))
    num.alternate.alleles[[2]] <<- num.alleles.in.src[[2]]-num.reference.alleles[[2]]
    # The alleles.state have changed so recaculte!
    ##browser()
    alleles.state <<- rep(rep(c(0,1), length(num.alleles.in.src[[2]])), as.vector(rbind(num.reference.alleles[[2]], num.alternate.alleles[[2]])))
    # Do it this way for now, but maybe should do this elsewhere so don't need to do <<
    #num.alternate.alleles <<- calculateAlternateAlleles()
    #num.reference.alleles <<- calculateReferenceAlleles()
    calculateAlternateAlleles()
    calculateReferenceAlleles()
    SSW <<- calculateSSW(1:(num.sources.of.variation+1))

#print(cat("bootstrapping SSW is", SSW))
  
##browser()
    #  A typical call to amova amova(new.srcNs, tsg, tn, TRUE) 
    #  A typical call to amova amova(new.srcNs, tn, TRUE)    # no longer pass in sg, only need n's
    # So we have to update those variables.
  }
  return()
}

    
##browser()

    #
    #  First we get the Amova Table, this will have everything except the p-values and the counts for the p-values
    
    #amovaTable <- amova(num.alleles.in.src, sg.noshuffle, n.noshuffle)
    amovaTable <- amova(num.alleles.in.src, n.noshuffle)


    #
    # Now we need to get the p-values.  This involves shuffling samples around via a shuffle method.

    if (maxPermutations < 1)  return(amovaTable)    # just return table and be done!


# ---------------------------------------------------------------------------------------------------------------------
#
#   exact  Permutation
    if (permutationMethod == "exact") { 
      cnt <- rep(0, num.sources.of.variation-1)   # set all the cnts to 0
      pcnt <- rep(0, num.sources.of.variation-1)  # set all the cnts to 0
      SSW.old <- SSW
      old.alternate.alleles <- num.alternate.alleles  # Since shuffles are moved around, we have to restore levels 2 through N each time through

#cat("Going to permumte for ", 1:maxPermutations, " and shuffle for ", variationSourcesPlusAllelesList, "\n")
      for(permutationi in 1:maxPermutations) {
	curPermutationi <<- permutationi   # don't know how to get it to function otherwise....

	if (min(cnt) >= NresamplesToStop) break  # done with shuffling all, so lets be done!
	# note the number of alleles doesn't change when bootstrap just the distrubtion of alts and refs.
        new.srcNs <- num.alleles.in.src  # restore the old  srcNs since doing a new set of shuffles  

	# First bootstrap the reads in the tube
	if (do.bootstrapReadsOld) {  # Bootstrap reads
	  bootstrapReadsOld(alts.shuffled.N)
	} else if(!do.bootstrapReads) {
	  SSW <<- SSW.old           # restore SSW
	  num.alternate.alleles <- old.alternate.alleles
	}


	for (index.of.src.being.shuffled in variationSourcesPlusAllelesList) {
#cat("permutation", permutationi, " shuffling index ", index.of.src.being.shuffled, "\n")
	  pcnt[index.of.src.being.shuffled] <- pcnt[index.of.src.being.shuffled] +1
	  level <- index.of.src.being.shuffled+1;

	  if (do.bootstrapReads) {
	    if (index.of.src.being.shuffled == 1) { 
	      #alts.shuffled.N<-bootstrapReads(alts.shuffled.N)
	      #num.alternate.alleles[[2]] <- bootstrapReads(alts.shuffled.N)
	      #num.alternate.alleles[[2]] <- alts.shuffled.N
#cat("alleles before bootstrap:", num.alternate.alleles[[2]], "\n")
	      num.alternate.alleles[[2]] <- bootstrapReads(shuffleSource1Function())
#cat("alleles after bootstrap:", num.alternate.alleles[[2]], "\n")
	      calculateAlternateAlleles()
	      #calculateReferenceAlleles()
	      SSW <<- calculateSSW(1:(num.sources.of.variation+1))
	      #SSW <<- ccalcSSW(num.alternate.alleles, num.alleles.in.src,1:(num.sources.of.variation+1)) 
	    } else {
	      # for source i
	      #   1) shuffle i-1 withing i+1
	      #   2) bootstrap reads for level 1 (alleles)
	      #   3) convert to alleles

	      #   1)
	      shuffledsrcs <- getSourceiShuffle(index.of.src.being.shuffled)
#cat("shuffledsrcs:", shuffledsrcs, "\n")
	      new.srcNs[[index.of.src.being.shuffled]] <- num.alleles.in.src[[index.of.src.being.shuffled]][shuffledsrcs]
#cat("new.srcNs:"); print(shuffledsrcs)
	      #   2) 3)
#cat("alleles before bootstrap:", num.alternate.alleles[[2]], "\n")
	      num.alternate.alleles[[2]] <- bootstrapReads(num.alternate.alleles[[2]])
#cat("alleles after bootstrap:", num.alternate.alleles[[2]], "\n")
	      # since number of alleles change in level 1, they can change in all levels so recalculate
	      calculateAlternateAlleles()
#cat("num.alternate.alleles after bootstrap:"); print(num.alternate.alleles)
	      # SSW therefore changes.
	      SSW <<- calculateSSW(1:(num.sources.of.variation+1))
	      #SSW <<- ccalcSSW(num.alternate.alleles, new.srcNs,1:(num.sources.of.variation+1)) 
#cat("SSW: ", SSW, "\n");
	    }
	  } else {
	    # When we shuffle src1, the alleles   the number of alleles 
	    if (index.of.src.being.shuffled == 1) {   #shuffle src1  or alleles
		alts.shuffled.N <-shuffleSource1Function()   # Get shuffled # of alternate alleles for each s1
		SSW[level]<<-sum(alts.shuffled.N*(num.alleles.in.src[[level]]-alts.shuffled.N)/num.alleles.in.src[[level]]) 
	    } else {   ##  if (index.of.src.being.shuffled >= 2)  
	      # restore the previous SSW - only need to restore the one that has changed!
	      # Bird method doesn't do this.  Is that really the only difference????
	      SSW[index.of.src.being.shuffled]<<-sum(num.alternate.alleles[[index.of.src.being.shuffled]]*(num.alleles.in.src[[index.of.src.being.shuffled]]-num.alternate.alleles[[index.of.src.being.shuffled]])/num.alleles.in.src[[index.of.src.being.shuffled]])  
	      shuffledsrcs <- getSourceiShuffle(index.of.src.being.shuffled)
	      new.srcNs[[index.of.src.being.shuffled]] <- num.alleles.in.src[[index.of.src.being.shuffled]][shuffledsrcs]
	      num.alternate.alleles[[level]]<-vapply(srcX.in.srcY[[index.of.src.being.shuffled]][[level]],  
			   function(x) sum(num.alternate.alleles[[index.of.src.being.shuffled]][shuffledsrcs[x]]),0)  
	      SSW[level] <<- sum(num.alternate.alleles[[level]]*(new.srcNs[[level]]-num.alternate.alleles[[level]])/new.srcNs[[level]])   
	    }
	  }
	  
	  if (index.of.src.being.shuffled >= 2) {   # srcNs change
	      oldtsg <- tsg 
	      tsg <- calcSG(new.srcNs)     ## 
	      oldtn <- tn
	      if (!identical(tsg, oldtsg)) {
		  tn <- calcN(tsg)             ## 
	      } else {
		  tn <- oldtn
	      }
	      foo <- amova(new.srcNs, tn, TRUE) #  
	      ##browser()
	  } else {
	      # if we shuffled the alleles within src1, then the numbers cannot change.
	      foo <- amova(num.alleles.in.src, n.noshuffle, TRUE)  #  Ns don't change
	      ##browser()
	  }
	  #fval <- foo[num.sources.of.variation-index.of.src.being.shuffled, 5]
	  # 1 ==>  4-1 = foo[3,5]
	  # 2 ==>  4-2 = foo[2,5]
	  # 3 ==>  4-3 = foo[1,5]
	  fval <- foo[index.of.src.being.shuffled+1]
#cat("fval:", fval, "AT:", amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5],"\n")
	  if (is.na(fval) || is.nan(fval)) break
          if (save.distributions) {
	    NullDistributionsTmp[[index.of.src.being.shuffled]][permutationi] <- fval
	  }
  ##browser()
	  if (fval >= amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]) {  
	    cnt[index.of.src.being.shuffled] <- cnt[index.of.src.being.shuffled] +1
	  }
#if ((curPermutationi %% 100) == 0)
if (curPermutationi %in% permcntsToReportOn)
    cat("permutations:", curPermutationi, "pcnt:",  pcnt[index.of.src.being.shuffled], 
    "cnt :", cnt[index.of.src.being.shuffled], " curp:", cnt[index.of.src.being.shuffled]/ pcnt[index.of.src.being.shuffled], "\n") 


#cat("permutation", permutationi, " Shuffled index ", index.of.src.being.shuffled, " (level):", level, "pcnt[",index.of.src.being.shuffled,"]:",
	  #pcnt[index.of.src.being.shuffled], " perm fval:", fval, " fval:", 
	  #amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5], "cnt[", index.of.src.being.shuffled,
	  #"]:", cnt[index.of.src.being.shuffled], " curp:", cnt[index.of.src.being.shuffled]/ pcnt[index.of.src.being.shuffled], "\n") 
	}  # end  for indexed.being.shuffled
#print(foo)

      }  # end for 1..maxPermutations
	  ## Now lets fill in the p-values and # permutations columns
      for (which.src in variationSourcesPlusAllelesList) {
        if (save.distributions) {
	  NullDistributions[[which.src]] <- NullDistributionsTmp[[which.src]][1:pcnt[which.src]]  # remove non used slots  
	}
	amovaTable[num.sources.of.variation-which.src,6] <- cnt[which.src]/pcnt[which.src]
#cat("cnt", cnt, "pcnt", pcnt, "\n")
	#if (is.na(amovaTable[num.sources.of.variation-which.src,6])) #browser()
	amovaTable[num.sources.of.variation-which.src,7] <- pcnt[which.src]
      }
	#print(amovaTable)
    } #end if (permutationMethod == "exact")


# ---------------------------------------------------------------------------------------------------------------------
#
#  freely   Purmutation
#
    ##
    ##  The first shuffle methods is the "freely" method
    #
    #   For the freely method, the alleles are freely shuffled among all groups of the current source level being tested, that is
    #   there is no restriction on the shuffle.  Only one shuffle needs to be done for all sources, and it is reused.
    #
    else if (permutationMethod == "freely") {
      cnt <- rep(0, num.sources.of.variation-1)    # set all the cnts to 0
      for(permutationi in 1:maxPermutations) {
	curPermutationi <<- permutationi           # don't know how to get it to function otherwise....

	# First bootstrap the reads in the tube
	if (do.bootstrapReadsOld) {
	    num.alternate.alleles[[2]]<-bootstrapReadsOld(alts.shuffled.N)
        }

	# Next shuffle source 1

	alts.shuffled.N <-shuffleSource1Function()   # Get shuffled # of alternate alleles for each s1
##browser()
#cat(paste("shuffle #", permutationi, "alts.shuffled.N are now"), alts.shuffled.N,"\n")
	# update SSW for source 1
	SSW[2]<<-sum(alts.shuffled.N*(num.alleles.in.src[[2]]-alts.shuffled.N)/num.alleles.in.src[[2]]) 
#print(paste("setting SSW[", 2, "] to", SSW[2]))
	# now update all the others
	if (do.bootstrapAlleles) {
	    num.alternate.alleles[[2]]<-bootstrapAlleles(alts.shuffled.N)
	}
	else if (do.bootstrapReads) {
	    num.alternate.alleles[[2]]<-bootstrapReads(alts.shuffled.N)
	}
	 else {
	    num.alternate.alleles[[2]]<-alts.shuffled.N
	}

	# Calculate the number of alternate alleles at each level and calulate SSW
	for (index.of.src.being.shuffled in variationSourcesList ) {
          if (is.nan(amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5])) {print("NaN freely f AT");  break }
	  level <- index.of.src.being.shuffled+1;
	  num.alternate.alleles[[level]]<-vapply(srcX.in.srcY[[level-1]][[level]], function(x) sum(num.alternate.alleles[[level-1]][x]),0)
#-cat(paste("num.alternate.allles[[", level-1,"]] are now"), num.alternate.alleles[[level-1]], "\n")
#cat(paste("num.alternate.allles[[", level,"]] are now"), num.alternate.alleles[[level]], "\n")
	  SSW[level]<<-sum(num.alternate.alleles[[level]]*(num.alleles.in.src[[level]]-num.alternate.alleles[[level]])/num.alleles.in.src[[level]])   
#print(paste("setting SSW[", level, "] to", SSW[level]))
	}
#print("SSW:")
#print(SSW)
	#foo <- amova(num.alleles.in.src, sg.noshuffle, n.noshuffle, TRUE)  # group counts don't change at all, use original
	foo <- amova(num.alleles.in.src, n.noshuffle, TRUE)  # group counts don't change at all, use original
##browser()
#print(foo)
	
	for (index.of.src.being.shuffled in variationSourcesPlusAllelesList) {
          if (is.nan(amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5])) {print("NaN freely f foo");  break }
	  level <- index.of.src.being.shuffled+1;
	  #fval <- foo[num.sources.of.variation-index.of.src.being.shuffled, 5]
	  fval <- foo[index.of.src.being.shuffled+1]
	  if (is.na(fval) || is.nan(fval)) break
          if (save.distributions) {
	    NullDistributionsTmp[[index.of.src.being.shuffled]][permutationi] <- fval
	  }
	  if (fval >= amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]) {  
	    cnt[index.of.src.being.shuffled] <- cnt[index.of.src.being.shuffled] +1
	  }
	}
	if (min(cnt) >= NresamplesToStop) break  # done with shuffling since all ps have enough observations >= f
      }  # end of foreach permutation

      # Put p-vals into amovaTable  and store Null Distributions
      for (index.of.src.being.shuffled in variationSourcesPlusAllelesList) {
        amovaTable[num.sources.of.variation-index.of.src.being.shuffled,6] <- cnt[index.of.src.being.shuffled]/permutationi
        amovaTable[num.sources.of.variation-index.of.src.being.shuffled,7] <- permutationi
        if (save.distributions) {
          NullDistributions[[index.of.src.being.shuffled]] <- NullDistributionsTmp[[index.of.src.being.shuffled]][1:permutationi]  # remove non used slots
	}
      }
    }  # end of if permutationMethod == "freely" 

#
#  Deal with the purmutation method "alleles"
#
#  For this method, only alleles get shuffled.  
#  When testing source i, the alleles are moved around within the confines of src i+1

#  So for each permutation, the sizes (Ns) of the groupings do not change, but the alt/ref allele counts do.
#  So that means the SSW from src 1 to src i+1 will change.

# ---------------------------------------------------------------------------------------------------------------------
#
#  alleles   Purmutation
#
    else if (permutationMethod == "alleles") {
      for (index.of.src.being.shuffled in variationSourcesPlusAllelesList) {
        level <- index.of.src.being.shuffled+1;
#-print(paste("Shuffling src", index.of.src.being.shuffled, "level", level))
#if (is.nan(amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]) ) {  #Can we have Nans?
#}else 
	cnt <- 0;

	for(permutationi in 1:maxPermutations) {
	  curPermutationi <<- permutationi   # don't know how to get it to function otherwise....
	  # Level one can actually use the individual speed up, so lets do that.

	  #num.alternate.alleles[[2]] <- shuffleSourceiFunction(index.of.src.being.shuffled)  
	  alts.shuffled.N <- shuffleSourceiFunction(index.of.src.being.shuffled)  
##browser()
	  if (do.bootstrapAlleles) {
	    num.alternate.alleles[[2]]<-bootstrapAlleles(alts.shuffled.N)
	  } else if (do.bootstrapReads) {
	    num.alternate.alleles[[2]]<-bootstrapReads(alts.shuffled.N)
	  } else {
	    num.alternate.alleles[[2]]<-alts.shuffled.N
	  }

#-print(paste("shuffle for src", index.of.src.being.shuffled, "is:"))
#-print(num.alternate.alleles[[2]])

##browser()
	  if (level > 2) {
	     num.alternate.alleles[3:level] <- calculateAlternateCounts(3:level) 
#-print(paste("now alternates counts for others are 3:",level, "is:"))
#-print(num.alternate.alleles[3:level])
	     }
          SSW[2:(level+1)] <<- calculateSSW(2:(level+1))

#-print("SSW is:")
#-print(SSW)

            #foo <- amova(num.alleles.in.src, sg.noshuffle, n.noshuffle, TRUE)  #  Ns don't change
            foo <- amova(num.alleles.in.src, n.noshuffle, TRUE)  #  Ns don't change
#-print(foo)
##browser()
            
            
#print(paste("comparing if (",foo[num.sources.of.variation-index.of.src.being.shuffled,5],">=",amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]))
	    #fval <- foo[num.sources.of.variation-index.of.src.being.shuffled, 5]
	    fval <- foo[index.of.src.being.shuffled+1]
	    if (is.na(fval) || is.nan(fval)) break
            if (save.distributions) {
	      NullDistributionsTmp[[index.of.src.being.shuffled]][permutationi] <- fval
	    }
#print(paste("fval ", fval, " vs ", amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]))  
	    if (fval >= amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]) {  
              if ((cnt<-cnt+1) == NresamplesToStop) {break}     # if we have 100 above or equal, no point in continuing
	    }
          }
          if (save.distributions) {
            NullDistributions[[index.of.src.being.shuffled]] <- NullDistributionsTmp[[index.of.src.being.shuffled]][1:permutationi]  # remove non used slots
	  }
          
          amovaTable[num.sources.of.variation-index.of.src.being.shuffled,6] <- cnt/permutationi
          amovaTable[num.sources.of.variation-index.of.src.being.shuffled,7] <- permutationi
        }
        
    } # end if (permutationMethod == "alleles") 

# ---------------------------------------------------------------------------------------------------------------------
#
#  Rexact    Purmutation
#
    else {#  not freely, bird, or alleles or exact!!.      So rexact 
cat("here 3\n")
      for (index.of.src.being.shuffled in variationSourcesPlusAllelesList) {
        level <- index.of.src.being.shuffled+1;
#-print(paste("Shuffling src", index.of.src.being.shuffled, "level", level))
        #if (is.nan(amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]) ) {  #Can we have Nans?
        #}else 
	cnt <- 0;
	if (index.of.src.being.shuffled >= 2) {
	  # restore the previous SSW - only need to restore the one that has changed!
	  SSW[index.of.src.being.shuffled]<<-sum(num.alternate.alleles[[index.of.src.being.shuffled]]*(num.alleles.in.src[[index.of.src.being.shuffled]]-num.alternate.alleles[[index.of.src.being.shuffled]])/num.alleles.in.src[[index.of.src.being.shuffled]])  
	  #print(paste("setting SSW[", index.of.src.being.shuffled, "] to", SSW[index.of.src.being.shuffled]))
	}
	new.srcNs <- num.alleles.in.src  # restore the old  srcNs
	
	for(permutationi in 1:maxPermutations) {
	  curPermutationi <<- permutationi   # don't know how to get it to function otherwise....
	  
	  if (index.of.src.being.shuffled == 1) {   #shuffle src1  or alleles
	    alts.shuffled.N <-shuffleSource1Function()   # Get shuffled # of alternate alleles for each s1      ## 

	    SSW[level]<<-sum(alts.shuffled.N*(num.alleles.in.src[[level]]-alts.shuffled.N)/num.alleles.in.src[[level]]) 
	  } else if (index.of.src.being.shuffled >= 2) {   #shuffle src2
	    # Shuffle the src1 groups, across the source2s
### SAK need to deal with preshuffle for other levels here!!!
	    shuffledsrcs <- getSourceiShuffle(index.of.src.being.shuffled)

	    new.srcNs[[index.of.src.being.shuffled]] <- num.alleles.in.src[[index.of.src.being.shuffled]][shuffledsrcs]
	    ## alleles move around so numbers change, but just totals  SSWs for all above can change
	    num.alternate.alleles[[level]]<-vapply(srcX.in.srcY[[index.of.src.being.shuffled]][[level]],  
			 function(x) sum(num.alternate.alleles[[index.of.src.being.shuffled]][shuffledsrcs[x]]),0)  

	    SSW[level] <<- sum(num.alternate.alleles[[level]]*(new.srcNs[[level]]-num.alternate.alleles[[level]])/new.srcNs[[level]])   
	  }


	  if (do.bootstrapAlleles) {
	    alts.shuffled.N<-bootstrapAlleles(alts.shuffled.N)
	  } 
	  # ---------------------------------------------
	  # bootstrap:
	  #
	  # take shuffled alleles, and then bootstrap the reads

	  if (do.bootstrapReads) {
	    alts.shuffled.N<-bootstrapReads(alts.shuffled.N)
	  } 

	  if (index.of.src.being.shuffled >= 2) {   # srcNs change
	      oldtsg <- tsg 
	      tsg <- calcSG(new.srcNs)     ## 
	      oldtn <- tn
	      if (!identical(tsg, oldtsg)) {
		  tn <- calcN(tsg)             ## 
	      } else {
		  tn <- oldtn
	      }
	      foo <- amova(new.srcNs, tn, TRUE) #  
	  } else {
	      foo <- amova(num.alleles.in.src, n.noshuffle, TRUE)  #  Ns don't change
	  }
	  #fval <- foo[num.sources.of.variation-index.of.src.being.shuffled, 5]
	  # 1 ==>  4-1 = foo[3,5]
	  # 2 ==>  4-2 = foo[2,5]
	  # 3 ==>  4-3 = foo[1,5]
	  fval <- foo[index.of.src.being.shuffled+1]
	  if (is.na(fval) || is.nan(fval)) break
          if (save.distributions) {
	    NullDistributionsTmp[[index.of.src.being.shuffled]][permutationi] <- fval
	  }
##browser()
	  if (fval >= amovaTable[num.sources.of.variation-index.of.src.being.shuffled,5]) {  
	    if ((cnt<-cnt+1) == NresamplesToStop) {break}     # if we have 100 above or equal, no point in continuing
	  }
	}
        if (save.distributions) {
	  NullDistributions[[index.of.src.being.shuffled]] <- NullDistributionsTmp[[index.of.src.being.shuffled]][1:permutationi]  # remove non used slots
	}
	
	amovaTable[num.sources.of.variation-index.of.src.being.shuffled,6] <- cnt/permutationi
  cat( "cnt:", cnt, " permutationi:", permutationi,"\n")
	amovaTable[num.sources.of.variation-index.of.src.being.shuffled,7] <- permutationi
      }  # end of for loop for source to permute
    } # End of if permutation method, elseif...

##browser()    
    
    if (save.distributions) 
      return (list(amovaTable, NullDistributions))
    else 
      return(amovaTable)
  }  # end  ProcessSNP
  
  if(multi.core==F){
    AllAmovaTables = lapply(1:(numSNPs), function(s) processSNP(s, maxPermutations=maxPermutations))
  } else if(multi.core==T | is.numeric(multi.core)){
    library(parallel)
    num.cores <- ifelse(multi.core==T,detectCores()-1,multi.core)
    
    if(Sys.info()['sysname']=='Windows'){
      my.cluster <- makeCluster(num.cores, type="PSOCK")
      vars <- c("DF", "DF.total","MS.total")
      clusterExport(my.cluster, vars)
    } else {
      my.cluster <- makeCluster(num.cores, type="FORK")
    }
    
    AllAmovaTables = parLapply(my.cluster, 1:(numSNPs), function(s) processSNP(s, maxPermutations=maxPermutations))
    stopCluster(my.cluster)
    
  } else {
    stop('multi.core argument must be logical or numeric')
  }


#  if (multi.node) {
  #   library(Rmpi)
  #   library(snow)
  # 
  #   np <- 39       # This should be a paramter!!!
  #   cluster <- makeMPIcluster(np)
  # 
  # 
  #   vars <- c("DF", "DF.total","MS.total")
  # 
  #   clusterExport(cluster, vars)
  #   AllAmovaTables = parLapply(cluster, 1:(numSNPs), function(s) processSNP(s, maxPermutations=maxPermutations))
  # 
  #   stopCluster(cluster)
  #   mpi.exit()
  # } else if (multi.core) {
  #   library(parallel)
  #   num.cores <- detectCores()-1
  #   print(num.cores)
  #   my.cluster <- makeCluster(num.cores, type="FORK")
  #   print(my.cluster)
  #   AllAmovaTables = parLapply(my.cluster, 1:(numSNPs), function(s) processSNP(s, maxPermutations=maxPermutations))
  #   stopCluster(my.cluster)
  # } else {
  #   AllAmovaTables = lapply(1:(numSNPs), function(s) processSNP(s, maxPermutations=maxPermutations))
  # }
  
  addColnames<-function(x){
    colnames(x)<-c('DF','SSa','MS','sig.sq','f','p','permutations')
    rownames(x)<-c(names(expDesign)[num.sources.of.variation:2],'Residual','Total')
    x
  }
  
  AllAmovaTables<-lapply(AllAmovaTables,addColnames)
  
  return(AllAmovaTables)
  
} # End of runAMOVA


