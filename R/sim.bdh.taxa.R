#' Simulate a Phylogenetic Network to a Specififed Number of Taxa
#'
#' @description Simulates a Phylogenetic Network under a birth-death-hybridization model. Simulates to a specified number of extant tips under the Simple Sampling Approach.
#'
#' @param n The number of taxa.
#' @param numbsim Number of netowrks to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction Rate
#' @param nu Hybridization Rate
#' @param hybprops Vector that represents the proportion of Hybridizations that are lineage generative, lineage degenerative, and lineage neutral respectively.
#' @param hyb.rate.function The probability of a successful hybridization as a function of genetic distance between taxa. The default value of NULL assumes that hybridizaion success is independent of genetic distance between taxa.
#' @param alpha The alpha parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param beta  The beta parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param complete If complete = TRUE, the tree with the extinct lineages is returned. If complete = FALSE, the extinct lineages are suppressed.
#' @param frac Sampling fraction: The proportion of extant tips included in the phylogeny (incomplete sampling).
#' @param stochsampling When stochsampling=TRUE: Each tip is included into the final tree with probability frac.
#'
#' @return out Returns a list of numbsim networks. Each network has an additional attirbute "inheritance" that represents the inheritance probabilities on the edges in the "reticulation" attribute.
#' @export
#'
#' @examples
sim.bdh.taxa.ssa <- function(n,numbsim,
                        lambda,mu,
                        nu, hybprops, hyb.rate.function=NULL,
                        alpha=1,beta=1,
                        frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE){
	out<-lapply(1:numbsim,FUN=sim.bdh.taxa.ssa.help,n=n,lambda=lambda,mu=mu,nu=nu,hybprops=hybprops,hyb.rate.function=hyb.rate.function,alpha=alpha,beta=beta,frac=frac,mrca=mrca,complete=complete,stochsampling=stochsampling)

	out
}



#' Simulate a Phylogenetic Network to a Specififed Number of Taxa
#'
#' @description Simulates a Phylogenetic Network under a birth-death-hybridization model. Simulates to a specified number of extant tips under the General Sampling Approach.
#'
#' @param m The number of taxa to simulate under the SSA
#' @param n The number of taxa.
#' @param numbsim Number of netowrks to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction Rate
#' @param nu Hybridization Rate
#' @param hybprops Vector that represents the proportion of Hybridizations that are lineage generative, lineage degenerative, and lineage neutral respectively.
#' @param hyb.rate.function The probability of a successful hybridization as a function of genetic distance between taxa. The default value of NULL assumes that hybridizaion success is independent of genetic distance between taxa.
#' @param alpha The alpha parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param beta  The beta parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param complete If complete = TRUE, the tree with the extinct lineages is returned. If complete = FALSE, the extinct lineages are suppressed.
#' @param frac Sampling fraction: The proportion of extant tips included in the phylogeny (incomplete sampling).
#' @param stochsampling When stochsampling=TRUE: Each tip is included into the final tree with probability frac.
#'
#' @return out Returns a list of numbsim networks. Each network has an additional attirbute "inheritance" that represents the inheritance probabilities on the edges in the "reticulation" attribute.
#' @export
#'
#' @examples

sim.bdh.taxa.gsa <-function(m,n,numbsim,
                        lambda,mu,
                        nu, hybprops, hyb.rate.function=NULL,
                        alpha=1,beta=1,
                        frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE){
  out<-lapply(1:numbsim,sim.bdh.taxa.gsa.help,m=m,n=n,lambda=lambda,mu=mu,nu=nu,hybprops=hybprops,hyb.rate.function=hyb.rate.function,alpha=alpha,beta=beta,frac=frac,mrca=mrca,complete=complete,stochsampling=stochsampling)

  out
}

