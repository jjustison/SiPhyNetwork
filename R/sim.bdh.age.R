#' Simulate a Phylogenetic Network to a Specified Number of Taxa
#'
#' @description Simulates a Phylogenetic Network under a birth-death-hybridization model. Simulates to a specified ages.
#'
#' @param age The time for each simulation.
#' @param numbsim Number of networks to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction rate.
#' @param nu Hybridization rate.
#' @param hybprops Vector that represents the proportion of Hybridizations that are lineage generative, lineage degenerative, and lineage neutral respectively.
#' @param hyb.inher.fxn A function for drawing the hybrid inheritance probabilities.
#' @param frac Sampling fraction: The proportion of extant tips included in the phylogeny (incomplete sampling).
#' @param twolineages If `twolineages=TRUE`: The process originates with two lineages that share a common ancestor. If `twolineages=FALSE`: The process originates with two lineages.
#' @param complete If complete = TRUE, the tree with the extinct lineages is returned. If complete = FALSE, the extinct lineages are suppressed.
#' @param stochsampling When `stochsampling=TRUE`: Each extant tip is included into the final tree with probability frac.
#' @param hyb.rate.fxn The probability of a successful hybridization as a function of genetic distance between taxa. The default value of `NULL`` assumes that hybridization success is independent of genetic distance between taxa.
#' @param trait.model A list that dictates how a trait affects the hybridization process. The default value of `NULL` doesn't take a trait into account for simulation. See Details for more information.
#' @return out Returns a list of numbsim networks with the time since origin / most recent common ancestor being 'age.' If tree goes extinct or no tips are sampled, return value is '0'. If only one extant and no extinct tips are sampled, return value is '1'. Each network has an additional attribute "inheritance" that represents the inheritance probabilities on the edges in the "reticulation" attribute.
#' @details `hyb.inher.fxn` should return values between 0 and 1 and shouldn't require any arguments. E.g. \link{make.beta.draw} and \link{make.uniform.draw} create functions that fit these specifications
#'
#' `hyb.rate.fxn` should take one argument for the genetic distance. The function should be defined on the range \eqn{[0,Inf)} and return values between \eqn{[0,1]}
#'
#' `trait.model` is a list with the following named elements:
#' \itemize{
#'  \item{`initial`}{ The initial trait state on the phylogeny}
#'  \item{`hyb.event.fxn`}{ A function that denotes the trait of a hybrid child after a hybridization event. The function should have the arguments `parent_states` and `inheritance`. `parent_states` is vector with the ploidy states of the hybrid parents while `inheritance` is the inheritance probability of the first lineage denoted in `parent_states`.}
#'  \item{`hyb.compatibility.fxn`}{ A function that describes when hybridization events can occur between two taxa based on their traits. The function should have the argument `parent_states`, a vector with the trait states of the two parents to the hybrid child. The function should return `TRUE` for when a hybridization event is allowed to proceed and `FALSE` otherwise.}
#'  \item{`time.fxn`}{ A function that describes how traits change over time. The function should have the arguments `trait_states` and `timestep` in that order. `trait_states` is a vector containing the ploidy of all taxa while `timestep` is the amount of time given for trait evolution. The function should return a vector with the updated ploidy states of all taxa.}
#'  \item{`spec.fxn`}{ A function that describes how the trait changes at speciation events.The function should have the argument `tip_state` which has the state of the lineage just before speciation. The function should return a vector with two values, one denoting the trait of each of the two new species after the event.}
#' }
#' @export
#'
#' @examples
#' set.seed(17) ##smallest Quartan prime as seed
#' #Generate a tree with extinct leaves
#' net<-sim.bdh.age(1,1,5,2,1,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),complete=TRUE)[[1]]

sim.bdh.age <-function(age,numbsim,
                      lambda,mu,
                      nu, hybprops,
                      hyb.inher.fxn,
                      frac=1,twolineages=FALSE,complete=TRUE,stochsampling=FALSE,
                      hyb.rate.fxn=NULL,
                      trait.model=NULL){
	out<-lapply(1:numbsim,sim.bdh.age.help,
	            age=age,
	            lambda=lambda,mu=mu,
	            nu=nu, hybprops=hybprops,
	            hyb.rate.fxn=hyb.rate.fxn,
	            hyb.inher.fxn=hyb.inher.fxn,
	            frac=frac,mrca=twolineages,
	            complete=complete, stochsampling=stochsampling,
	            trait.model=trait.model)
	class(out)<-c('list','multiPhylo')
	out
}

