##R4.1.0
rm(list = ls())
#source("R/CompCCAT.R")
#source("R/InferDMAPandRoot.R")
#source("R/InferPotencyStates.R")

library(igraph)
library(mclust)
library(parallel)
library(stats)
library(Matrix)

#' @title 
#' Computes potency estimates of single-cells 
#' using the signaling entropy rate
#' 
#' @aliases CompSRana
#'  
#' @description 
#' This is the main user function for computing signaling entropy of 
#' single cells. It takes as input the gene expression profile of 
#' single cells and the adjacency matrix of a connected network. These 
#' inputs will be typically the output of the \code{DoIntegPPI} function.
#' 
#' @param integ.l
#' The output from \code{DoIntegPPI} function.
#' 
#' @param local
#' A logical (default is FALSE). If TRUE, function returns the normalized 
#' local signaling entropies of each gene in the network alongside the unnormalized
#' local entropies. If FALSE, only unnormalized local entropies are returned.
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many child processes will 
#' be run simultaneously. The option is initialized from environment variable 
#' MC_CORES if set. Must be at least one, and parallelization requires at 
#' least two cores.
#' 
#' @return A list with four elements:
#' 
#' @return SR
#' The global signaling entropy rate. It is normalized by the 
#' maximum rate, hence a value between 0 and 1
#' 
#' @return inv
#' The stationary distribution of every sample
#' 
#' @return locS
#' The unnormalised local entropies of each gene in every cell
#' 
#' @return nlocS
#' The normalised local entropies of each gene, so that each value is 
#' between 0 and 1
#' 
#' @references 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell's transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Teschendorff AE, Banerji CR, Severini S, Kuehn R, Sollich P. 
#' \emph{Increased signaling entropy in cancer requires the scale-free 
#' property of protein interaction networks.}
#' Scientific reports 5 (2015): 9646.
#' doi:\href{https://doi.org/10.1038/srep09646}{
#' 10.1038/srep09646}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Intra-tumour signalling entropy determines clinical outcome 
#' in breast and lung cancer.}
#' PLoS computational biology 11.3 (2015): e1004115.
#' doi:\href{https://doi.org/10.1371/journal.pcbi.1004115}{
#' 10.1371/journal.pcbi.1004115}.
#' 
#' Teschendorff, Andrew E., Peter Sollich, and Reimer Kuehn.
#' \emph{Signalling entropy: A novel network-theoretical framework 
#' for systems analysis and interpretation of functional omic data.}
#' Methods 67.3 (2014): 282-293.
#' doi:\href{https://doi.org/10.1016/j.ymeth.2014.03.013}{
#' 10.1016/j.ymeth.2014.03.013}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Cellular network entropy as the energy potential in 
#' Waddington's differentiation landscape.}
#' Scientific reports 3 (2013): 3039.
#' doi:\href{https://doi.org/10.1038/srep03039}{
#' 10.1038/srep03039}.
#' 
#' @examples 
#'
#' 
#' @import parallel
#' 
#' @export
#'
CompSRana <- function(integ.l, local = FALSE, mc.cores=1){
  
  ### compute maxSR for SR normalization
  maxSR <- CompMaxSR(integ.l);
  
  idx.l <- as.list(seq_len(ncol(integ.l$expMC)));
  out.l <- mclapply(idx.l, CompSRanaPRL, 
                    exp.m=integ.l$expMC, 
                    adj.m=integ.l$adjMC,
                    local=local,
                    maxSR=maxSR,
                    mc.cores=mc.cores)
  SR.v <- sapply(out.l, function(v) return(v[[1]]))
  invP.v <- sapply(out.l, function(v) return(v[[2]]))
  S.v <- sapply(out.l, function(v) return(v[[3]]))
  NS.v <- sapply(out.l, function(v) return(v[[4]]))
  
  return(list(SR=SR.v,inv=invP.v,locS=S.v,nlocS=NS.v));
}

#### Auxilliary functions
#' @import igraph
#' 
CompMaxSR <- function(integ.l){
  
  adj.m <- integ.l$adjMC
  
  # find right eigenvector of adjacency matrix
  fa <- function(x,extra=NULL) {
    as.vector(adj.m %*% x)
  }
  ap.o <- igraph::arpack(fa,options=list(n=nrow(adj.m),nev=1,which="LM"), sym=TRUE)
  v <- ap.o$vectors
  lambda <- ap.o$values
  
  # maximum entropy
  maxSR <- log(lambda)
  
  return(maxSR);
}

CompSRanaPRL <- function(idx,
                         exp.m,
                         adj.m,
                         local=TRUE,
                         maxSR=NULL)
{
  
  # compute outgoing flux around each node
  exp.v <- exp.m[,idx];
  sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
  invP.v <- exp.v*sumexp.v;
  nf <- sum(invP.v);
  invP.v <- invP.v/nf;
  p.m <- t(t(adj.m)*exp.v)/sumexp.v;
  S.v <- apply(p.m,1,CompS);
  SR <- sum(invP.v*S.v);
  # if provided then normalise relative to maxSR
  if(is.null(maxSR)==FALSE){
    SR <- SR/maxSR;
  }
  if(local){
    NS.v <- apply(p.m,1,CompNS);
  }
  else {
    NS.v <- NULL;
  }
  return(list(sr=SR,inv=invP.v,locS=S.v,nlocS=NS.v));
}

CompNS <- function(p.v){
  
  tmp.idx <- which(p.v>0);
  if(length(tmp.idx)>1){
    NLS <- -sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
  }
  else {
    # one degree nodes have zero entropy, avoid singularity.
    NLS <- 0;
  }
  return(NLS);
}

CompS <- function(p.v){
  
  tmp.idx <- which(p.v>0);
  LS <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
  return(LS);
}
#' @title 
#' Integration of gene expression matrix and PPI network
#' 
#' @aliases DoIntegPPI
#'  
#' @description 
#' This function finds the common genes between the scRNA-Seq data matrix 
#' and the genes present in the PPI network, and constructs the maximally 
#' connected subnetwork and associated expression matrix for the computation 
#' of signaling entropy.
#' 
#' @param exp.m
#' The scRNA-Seq data matrix normalized for library size and log2-transformed
#' with a pseudocount of 1.1
#' 
#' @param ppiA.m
#' The adjacency matrix of a user-given PPI network with rownames and 
#' colnames labeling genes (same gene identifier as in \code{exp.m})
#' 
#' 
#' @return A list of two or four objects:
#' 
#' @return expMC
#' Reduced expression matrix with genes in the maximally connected subnetwork
#' 
#' @return adjMC
#' Adjacency matrix of the maximally connected subnetwork
#' 
#' 
#'  
#' @references 
#' Chen, Weiyan, et al.
#' \emph{Single-cell landscape in mammary epithelium reveals 
#' bipotent-like cells associated with breast cancer risk 
#' and outcome.}
#' Communications Biology 2 (2019): 306.
#' doi:\href{https://doi.org/10.1038/s42003-019-0554-8}{
#' 10.1038/s42003-019-0554-8}. 
#' 
#' Teschendorff AE, Tariq Enver. 
#' \emph{Single-cell entropy for accurate estimation of differentiation 
#' potency from a cell's transcriptome.}
#' Nature communications 8 (2017): 15599.
#' doi:\href{https://doi.org/10.1038/ncomms15599}{
#' 10.1038/ncomms15599}.
#' 
#' Teschendorff AE, Banerji CR, Severini S, Kuehn R, Sollich P. 
#' \emph{Increased signaling entropy in cancer requires the scale-free 
#' property of protein interaction networks.}
#' Scientific reports 5 (2015): 9646.
#' doi:\href{https://doi.org/10.1038/srep09646}{
#' 10.1038/srep09646}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Intra-tumour signalling entropy determines clinical outcome 
#' in breast and lung cancer.}
#' PLoS computational biology 11.3 (2015): e1004115.
#' doi:\href{https://doi.org/10.1371/journal.pcbi.1004115}{
#' 10.1371/journal.pcbi.1004115}.
#' 
#' Teschendorff, Andrew E., Peter Sollich, and Reimer Kuehn.
#' \emph{Signalling entropy: A novel network-theoretical framework 
#' for systems analysis and interpretation of functional omic data.}
#' Methods 67.3 (2014): 282-293.
#' doi:\href{https://doi.org/10.1016/j.ymeth.2014.03.013}{
#' 10.1016/j.ymeth.2014.03.013}.
#' 
#' Banerji, Christopher RS, et al. 
#' \emph{Cellular network entropy as the energy potential in 
#' Waddington's differentiation landscape.}
#' Scientific reports 3 (2013): 3039.
#' doi:\href{https://doi.org/10.1038/srep03039}{
#' 10.1038/srep03039}.
#' 
#' @examples 
#'
#' 
#' @import Matrix
#' @importFrom igraph graph.adjacency
#' @importFrom igraph clusters
#' @export
#'     
DoIntegPPI <- function(exp.m, 
                       ppiA.m)
{
  
  if(max(exp.m) > 100){ ### check if data has been log2-transformed or not and if not, then log2-transform with a pseudcount of 1.1, as for SR computation later we are not allowed 0's.
    exp.m <- log2(exp.m+1.1);
  }
  
  # set common gene IDs
  commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m));
  
  # check the consistency of gene IDs and size of overlap
  if ( length(commonEID.v) < 5000 ){
    stop(paste("The overlap of common genes between PPI and expression matrix is only ",length(commonEID.v)," so check that gene identifiers are correct. We don't recommend running CCAT on less than 5000 overlapping genes.",sep=""));
  }
  
  
  # get maximum connected network
  match(commonEID.v,rownames(ppiA.m)) -> map2.idx
  adj.m <- ppiA.m[map2.idx, map2.idx]
  
  gr.o <- igraph::graph.adjacency(adj.m,mode="undirected")
  comp.l <- igraph::clusters(gr.o)
  cd.v <- summary(factor(comp.l$member))
  mcID <- as.numeric(names(cd.v)[which.max(cd.v)])
  maxc.idx <- which(comp.l$member==mcID)
  adjMC.m <- adj.m[maxc.idx, maxc.idx]
  
  match(rownames(adjMC.m),rownames(exp.m)) -> map1.idx
  expMC.m <- exp.m[map1.idx ,]
  
  return(list(expMC=expMC.m,adjMC=adjMC.m));        
  
} ### EOF


load("data/net17Jan16.rda")
print(dim(net17Jan16.m))
load("../SLICE_run/data/hs_km.Rda")

###integrate
library(org.Hs.eg.db)
library(annotate)
mapping<-lookUp(rownames(net17Jan16.m), 'org.Hs.eg', 'SYMBOL')
mapping<-mapping[mapping %in% rownames(km)]
net17Jan16.m<-net17Jan16.m[names(mapping), names(mapping)]
km<-km[as.character(unlist(mapping)), as.character(unlist(mapping))]
km_geneID<-km; rownames(km_geneID)<-rownames(net17Jan16.m);colnames(km_geneID)<-colnames(net17Jan16.m);

####Let us check the mapping
head(km[, 1:5])
head(net17Jan16.m[, 1:5])
####net17Jan16.m.GO<-net17Jan16.m * (1+km)/2  #km    We can not use km because the network will be unconnected
####Filter outun-connected, at least have two connections
####net17Jan16.m.GO<-net17Jan16.m.GO[rowSums(net17Jan16.m.GO>0)>10,rowSums(net17Jan16.m.GO>0)>10]

load("data/dataChu.rda")
print(dim(scChu.m));

lscChu.m <- log2(scChu.m+1.1)
range(lscChu.m)
## [1]  0.1375035 16.5819775
lscChu0.m <- log2(scChu.m+1)
range(lscChu0.m)

integ.l <- DoIntegPPI(exp.m = lscChu.m, ppiA.m = net17Jan16.m)
str(integ.l)
maxSR<-CompMaxSR(integ.l)
result<-CompSRanaPRL(1, exp.m = integ.l$expMC, adj.m = integ.l$adjMC, local = FALSE, maxSR = maxSR)
###Weighted with GO similarity
integ.l$adjMC<-integ.l$adjMC * (1+km_geneID[rownames(integ.l$adjMC), colnames(integ.l$adjMC)])/2


head(integ.l$adjMC[,1:5])
#sr.o <- CompSRana(integ.l, local = FALSE, mc.cores = 1)  #4)
maxSR<-CompMaxSR(integ.l)
result_nonInteger<-CompSRanaPRL(1, exp.m = integ.l$expMC, adj.m = integ.l$adjMC, local = FALSE, maxSR = maxSR)

results<-list()
#for(ii in 1:dim(integ.l$expMC)[2]){
for(ii in c(1:479)){
  cat(ii);cat("\n")
  results[[ii]]<-CompSRanaPRL(ii, exp.m = integ.l$expMC, adj.m = integ.l$adjMC, local = FALSE, maxSR = maxSR)
}
save(results, file = "results_PPIWithGOSimilarity.RData")

resultsNormal<-list()
#for(ii in 1:dim(integ.l$expMC)[2]){
for(ii in c(1:479)){
  cat(ii);cat("\n")
  resultsNormal[[ii]]<-CompSRanaPRL(ii, exp.m = integ.l$expMC, adj.m = integ.l$adjMC, local = TRUE, maxSR = maxSR)
}
save(resultsNormal, file = "results_PPIWithGOSimilarity_normalized.RData")


SR.v <- unlist(sapply(results, function(v) return(v[[1]])))
boxplot(SR.v ~ phenoChu.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
library(pROC)
roc1<-roc(phenoChu.v, SR.v); plot(roc1)

SR.v.normal <- unlist(sapply(resultsNormal, function(v) return(v[[1]])))
boxplot(SR.v.normal ~ phenoChu.v, main = "SR potency estimates", xlab = "Cell Type", ylab = "SR")
roc1<-roc(phenoChu.v, SR.v.normal); plot(roc1)


