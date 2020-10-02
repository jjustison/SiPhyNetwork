#' Simulate a Phylogenetic Network to a Specififed Number of Taxa
#'
#' @description Simulates a Phylogenetic Network under a birth-death-hybridization model. Simulates to a specified number of extant tips.
#'
#' @param n The number of taxa.
#' @param numbsim Number of netowrks to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction Rate
#' @param nu Hybridization Rate
#' @param hybprops Vector that represents the proportion of Hybridizations that are lineage generative, lineage degenerative, and lineage neutral respectively.
#' @param alpha The alpha parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param beta  The beta parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param complete If complete = TRUE, the tree with the extinct and non-sampled lineages is returned. If complete = FALSE, the extinct and non-sampled lineages are suppressed.
#' @param frac Sampling fraction: The actual number of tips is n/frac. If Stochsampling=FALSE only n tips are included (incomplete sampling). If stochsampling=TRUE then only n tips are included on average but is not garaunteed.
#' @param stochsampling When stochsampling=TRUE: Each tip is included into the final tree with probability frac.
#'
#' @return out Returns a list of numbsim networks. Each network has an additional attirbute "inheritance" that represents the inheritance probabilities on the edges in the "reticulation" attribute.
#' @export
#'
#' @examples
sim.bdh.taxa <-function(n,numbsim, complete=TRUE,
                           lambda,mu,
                           nu, hybprops,
                           alpha=1,beta=1,
                           frac=1,stochsampling=FALSE){
	out<-lapply(1:numbsim,sim.bdh.taxa.help,n=n,lambda=lambda,mu=mu,nu=nu,hybprops=hybprops,alpha=alpha,beta=beta,frac=frac,complete=complete,stochsampling=stochsampling)
	out1<-lapply(out,function(x){ x[[1]][[1]]})
	out2<-sapply(out,function(x){ x[[2]]})
	for (i in 1:length(out1)){
		out1[[i]]$root.edge<-out2[i]-max(getx(out1[[i]],sersampling=1)[,1])
	}
	#out3<-list(out1,out2)
	#out3
	out1
	}
