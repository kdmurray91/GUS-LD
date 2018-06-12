#########################################################################
# Genotyping Uncertainty with Sequencing data and Linkage Disequilibrium (GUS-LD)
# Copyright 2017 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#########################################################################
##### File for computing LD estimates using high-throughtput sequencing data
### Author: Timothy P. Bilton
### Date: 07/03/18

######### Notes ###########
# We require:
#   > the data set is coded in terms of 0,1 and 2
#     representing the number of alleles present for the major allele
#   > The allele counts for the reference and alternate allele are known

require(boot)
require(foreach)
require(doFuture)


###### Define functions

## likelihood functions
ll_wrapper <- function(logit_x, depth_Ref, depth_Alt, pAB){
  pA1=inv.logit(logit_x[2])
  pA2=inv.logit(logit_x[3])
  ep =inv.logit(logit_x[4])
  r = 2*inv.logit(logit_x[1]) - 1
  D = r*sqrt(pA1*pA2*(1-pA1)*(1-pA2))
  return(ll(c(pA1,pA2,D,ep), depth_Ref, depth_Alt, pAB))
}
ll <- function(x, depth_Ref, depth_Alt, pAB){
  pA1=x[1]; pA2=x[2]; D=x[3]
  pB1 = 1-pA1; pB2 = 1-pA2
  C1 <- max(-pA1*pA2,-pB1*pB2); C2 <- min(pB1*pA2,pA1*pB2)
  if(!((C1 <= D) & (D <= C2))) return(NA)
  else{
    ep = x[4]
    ## compute the probability matrices
    pAA <- choose(depth_Ref+depth_Alt,depth_Ref) * (1-ep)^depth_Ref * ep^depth_Alt
    pBB <- choose(depth_Ref+depth_Alt,depth_Ref) * (ep)^depth_Ref * (1-ep)^depth_Alt

    h1 <- pA1*pA2+D; h2 <- pA1*pB2-D; h3 <- pB1*pA2-D; h4 <- pB1*pB2+D
    p11 <- h1^2
    p21 <- 2*h1*h3
    p31 <- h3^2
    p12 <- 2*h1*h2
    p32 <- 2*h3*h4
    p13 <- h2^2
    p23 <- 2*h2*h4
    p33 <- h4^2
    p22 <- 2*h1*h4+2*h2*h3

    return(-(sum( log(
      pAA[,1]*pAA[,2]*p11 + pAA[,1]*pAB[,2]*p12 + pAA[,1]*pBB[,2]*p13 +
        pAB[,1]*pAA[,2]*p21 + pAB[,1]*pAB[,2]*p22 + pAB[,1]*pBB[,2]*p23 +
        pBB[,1]*pAA[,2]*p31 + pBB[,1]*pAB[,2]*p32 + pBB[,1]*pBB[,2]*p33
    )))  )
  }
}

## function for estimating LD estimates for given MLE values
Iden_D <- function(pA1,pA2,D){
  return(D)
}
Dprime <- function(pA1,pA2,D){
  C1 = max(-prod(c(pA1,pA2)),-prod(1-c(pA1,pA2))); C2 = min((1-pA1)*pA2,pA1*(1-pA2))
  return(D/ifelse(D<0,C1,C2))
}
r2 <- function(pA1,pA2,D){
  return(D^2/(prod(c(pA1,pA2,1-c(pA1,pA2)))))
}

## needed for foreach loop
comb <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

############ Main Function ########################
#### Computing matrics of the GBS Dprime and r2 estimates
# Arguments are
#  - GBSdata: Matrix of GBS genotype calls. 0=BB, 1=AB, 2=AA, NA=missing
#  - depth_Ref: Matrix of counts for the reference allele
#  - depth_Alt: Matrix of counts for the alternate allele
#  - parallel: If TRUE then parallel computing is used to speed up the computations
#  - Nclust: Number of clusters the foreach loop using for the parallel computing
#            Note: make sure you don't use for than the number of cores available
#  - LDmeasure: A vector or function numbers which compute the required LD measure from the MLE of D
#               Note: Some predefined functions are above.


#' Title
#'
#' @param snp1 SNP index 1
#' @param genon Genotype matrix
#' @param depth_Ref Ref. allele depth
#' @param depth_Alt Alt. allele depth
#' @param freq Allele frequencies
#' @param LDmeasures LD measures (\in ("Dprime", "R2"))
GUSLD_core <- function(snp1, genon,depth_Ref,depth_Alt, freq, LDmeasures) {
  LDvec <- c(replicate(length(LDmeasure),numeric(nSnps),simplify=F),list(numeric(nSnps)))
  for(snp2 in seq_len(snp1-1)){
    ind <- c(snp1,snp2)
    ## extract data we want
    depth_Ref2 <- depth_Ref[,ind]
    depth_Alt2 <- depth_Alt[,ind]

    pAB <- choose(depth_Ref2+depth_Alt2,depth_Ref2) * (1/2)^(depth_Ref2 + depth_Alt2)
    N <- length(apply(depth_Ref2==0&depth_Alt2==0,1,any))

    MLE <- optim(logit(c(0.5,freq[ind],0.01)),ll_wrapper,
                 method="Nelder-Mead", control=list(reltol=1e-10, maxit=2000),
                 depth_Ref=depth_Ref2,depth_Alt=depth_Alt2,pAB=pAB)
    MLE2 <- optim(logit(c(0.5,freq[ind],0.0001)),ll_wrapper,
                  method="Nelder-Mead", control=list(reltol=1e-10, maxit=2000),
                  depth_Ref=depth_Ref2,depth_Alt=depth_Alt2,pAB=pAB)
    if(MLE$value < MLE2$value)
      MLE <- MLE$par
    else
      MLE <- MLE2$par
    pA1_hat = inv.logit(MLE[2])
    pA2_hat = inv.logit(MLE[3])
    D_hat <- (2*inv.logit(MLE[1])-1)*sqrt(pA1_hat*pA2_hat*(1-pA1_hat)*(1-pA2_hat))
    D_hat <- D_hat*(2*N/(2*N-1))
    ep_hat = inv.logit(MLE[4])
    ## Check that D is within its range
    C1hat <- max(-(pA1_hat*pA2_hat),-prod(1-c(pA1_hat,pA2_hat))); C2hat <- min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
    D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))

    for(meas in 1:length(LDmeasure)) {
      LDvec[[meas]][snp2] <- get(LDmeasure[[meas]])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
    }
    LDvec[[length(LDmeasure)+1]][snp2] <- ep_hat
  }
  return(lapply(LDvec,round,digits=8))

}

#' Estimate LD using the GUS-LD model
#'
#' @param genon 012-coded genotype matrix. SNPs are columns.
#' @param depth_Ref Reference/major allele depth, same dims as `genon`
#' @param depth_Alt Alternate/minor allele depth, same dims as `genon`
#' @param LDmeasure Measure of LD. One or both of 'Dprime' and 'r2' (case sensitive)
#'
#' @return list with slots for each LDmeasure
#' @export
GUS_LD <- function(genon, depth_Ref, depth_Alt, LDmeasure=c('Dprime','r2')){
  ## Do some checks
  if(!is.matrix(depth_Ref)|!is.matrix(depth_Alt))
    stop("Matrix objects required for depth inputs")
  if(any(na.omit(as.vector(depth_Ref)<0)) & any(na.omit(as.vector(depth_Alt)<0)))
    stop("Negative depths are invalid.")
  if(any(dim(depth_Ref)!=dim(depth_Alt)))
    stop("GBS data and depths matrices are different dimensions")
  if(!(length(LDmeasure) > 0))
    stop("No LD measure has been specified")
  # If to run it using parallel programming
  freq <- colMeans(genon,na.rm=T)/2
  if(any(freq > 0.99))
    freq[which(freq > 0.99)] <- 0.99
  if(any(freq < 0.01))
    freq[which(freq < 0.01)] <- 0.01
  nSnps <- ncol(genon)

  status <- try(
    LDmat <- foreach(snp1=iter(1:nSnps),.combine=comb,.multicombine=TRUE,
                     .export=c("ll","ll_wrapper",LDmeasure, "LDmeasure", "freq", "genon", "depth_Ref", "depth_Alt"), .packages=c("GUSLD", "boot")) %dopar% {

      nSnps <- ncol(genon)
      LDvec <- c(replicate(length(LDmeasure),numeric(nSnps),simplify=F),list(numeric(nSnps)))
      for(snp2 in seq_len(snp1-1)){
        ind <- c(snp1,snp2)
        ## extract data we want
        depth_Ref2 <- depth_Ref[,ind]
        depth_Alt2 <- depth_Alt[,ind]

        pAB <- choose(depth_Ref2+depth_Alt2,depth_Ref2) * (1/2)^(depth_Ref2 + depth_Alt2)
        N <- length(apply(depth_Ref2==0&depth_Alt2==0,1,any))

        MLE <- optim(logit(c(0.5,freq[ind],0.01)),ll_wrapper,
                     method="Nelder-Mead", control=list(reltol=1e-10, maxit=2000),
                     depth_Ref=depth_Ref2,depth_Alt=depth_Alt2,pAB=pAB)
        MLE2 <- optim(logit(c(0.5,freq[ind],0.0001)),ll_wrapper,
                      method="Nelder-Mead", control=list(reltol=1e-10, maxit=2000),
                      depth_Ref=depth_Ref2,depth_Alt=depth_Alt2,pAB=pAB)
        if(MLE$value < MLE2$value)
          MLE <- MLE$par
        else
          MLE <- MLE2$par
        pA1_hat = inv.logit(MLE[2])
        pA2_hat = inv.logit(MLE[3])
        D_hat <- (2*inv.logit(MLE[1])-1)*sqrt(pA1_hat*pA2_hat*(1-pA1_hat)*(1-pA2_hat))
        D_hat <- D_hat*(2*N/(2*N-1))
        ep_hat = inv.logit(MLE[4])
        ## Check that D is within its range
        C1hat <- max(-(pA1_hat*pA2_hat),-prod(1-c(pA1_hat,pA2_hat))); C2hat <- min((1-pA1_hat)*pA2_hat,pA1_hat*(1-pA2_hat))
        D_hat <- ifelse(D_hat>=0,min(D_hat,C2hat),max(D_hat,C1hat))

        for(meas in 1:length(LDmeasure)) {
          LDvec[[meas]][snp2] <- get(LDmeasure[[meas]])(pA1=pA1_hat,pA2=pA2_hat,D=D_hat)
        }
        LDvec[[length(LDmeasure)+1]][snp2] <- ep_hat
      }
      return(lapply(LDvec,round,digits=8))

                     }
                )

  if(class(status)=='try-error'){
    stop('Error: Problem with estimation proceedure using parallel computing.')
  }

  names(LDmat) <- LDmeasure

  for(i in 1:(length(LDmeasure))){
    diag(LDmat[[i]]) <- get(LDmeasure[i])(0.5,0.5,0.25)
    LDmat[[i]][upper.tri(LDmat[[i]])] <- t(LDmat[[i]])[upper.tri(t(LDmat[[i]]))]
  }
  diag(LDmat[[i+1]]) <- NA
  LDmat[[i+1]][upper.tri(LDmat[[i+1]])] <- t(LDmat[[i+1]])[upper.tri(t(LDmat[[i+1]]))]

  return(LDmat)
}
