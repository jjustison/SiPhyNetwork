#' Simulate a Phylogenetic Network to a Specififed Number of Taxa
#'
#' @description Simulates a Phylogenetic Network under a birth-death-hybridization model. Simulates to a specified ages.
#'
#' @param age The time for each simulation.
#' @param numbsim Number of networks to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction rate.
#' @param nu Hybridization rate.
#' @param hybprops ector that represents the proportion of Hybridizations that are lineage generative, lineage degenerative, and lineage neutral respectively.
#' @param alpha The alpha parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param beta The beta parameter of the beta distribution. Used to determine inheritance probabilities.
#' @param frac Sampling fraction: The actual number of tips is n/frac, but only n tips are included (incomplete sampling).
#' @param mrca If mrca=FALSE: age is the time since origin. If mrca=TRUE: age is the time since most recent common ancestor of the extant tips.
#' @param complete If complete = TRUE, the tree with the extinct and non-sampled lineages is returned. If complete = FALSE, the extinct and non-sampled lineages are suppressed.
#' @param stochsampling When stochsampling=TRUE: Each tip is included into the final tree with probability frac.
#'
#' @return out Returns a list of numbsim networks with the time since origin / most recent common ancestor being 'age.' If tree goes extinct or no tips are sampled, return value is '0'. If only one extant and no extinct tips are sampled, return value is '1'. Each network has an additional attirbute "inheritance" that represents the inheritance probabilities on the edges in the "reticulation" attribute.
#' @export
#'
#' @examples
sim.bdh.age <-function(age,numbsim,
                      lambda,mu,
                      nu, hybprops, hyb.rate.function=NULL,
                      alpha=1,beta=1,
                      frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE){
	out<-lapply(1:numbsim,sim.bdh.age.help,age=age,lambda=lambda,mu=mu,frac=frac,mrca=mrca,complete=complete,nu=nu, hybprops=hybprops,hyb.rate.function=hyb.rate.function,alpha=alpha,beta=beta,stochsampling=stochsampling)
	out
}

sim.bdh.age(age=1,numbsim = 100,
            lambda = 1,mu=0.5,
            nu=0.5,hybprops = c(0.33,0.33,0.33),  hyb.rate.function=NULL,
            alpha = 1,beta = 1,
            frac=0.75,mrca=FALSE,complete = T,stochsampling = F)
