#' Simulate a Phylogenetic Network to a Specififed Number of Taxa
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
#' @param mrca If `mrca=FALSE`: age is the time since origin. If `mrca=TRUE`: age is the time since most recent common ancestor of the extant tips.
#' @param complete If complete = TRUE, the tree with the extinct lineages is returned. If complete = FALSE, the extinct lineages are suppressed.
#' @param stochsampling When `stochsampling=TRUE`: Each extant tip is included into the final tree with probability frac.
#' @param hyb.rate.fxn The probability of a successful hybridization as a function of genetic distance between taxa. The default value of `NULL`` assumes that hybridizaion success is independent of genetic distance between taxa.
#' @param trait.model A list that dictates how a trait affects the hybridization process. The default value of `NULL` doesn't take a trait into account for simulation. See Details for more information.
#' @return out Returns a list of numbsim networks with the time since origin / most recent common ancestor being 'age.' If tree goes extinct or no tips are sampled, return value is '0'. If only one extant and no extinct tips are sampled, return value is '1'. Each network has an additional attirbute "inheritance" that represents the inheritance probabilities on the edges in the "reticulation" attribute.
#' @details `hyb.inher.fxn` should return values between 0 and 1 and shouldn't require any arguments. E.g. \link{make.beta.draw} and \link{make.uniform.draw} create functions that fit these specifications
#'
#' `hyb.rate.fxn` should take one argument for the genetic distance. The function should be defined on the range \eqn{[0,Inf)} and return values between \eqn{[0,1]}
#'
#' `trait.model` is a list with the following named elements:
#' \itemize{
#'  \item{`initial`}{ The initial trait state on the phylogeny}
#'  \item{`hyb.event.fxn`}{ A function that denotes the trait of a hybrid child after a hybridization event.}
#'  \item{`hyb.compatability.fxn`}{ A function that describes when hybridization events can occur between two taxa based on their traits.}
#'  \item{`time.fxn`}{ A function that describes how traits change over time.}
#'  \item{`spec.fxn`}{ A function that describes how the trait changes at speciation events.}
#' }
#' See \link{make.polyploid.model} for more information about using and creating this list.
#' @export
#'
#' @examples
sim.bdh.age <-function(age,numbsim,
                      lambda,mu,
                      nu, hybprops,
                      hyb.inher.fxn,
                      frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE,
                      hyb.rate.fxn=NULL,
                      trait.model){
	out<-lapply(1:numbsim,sim.bdh.age.help,
	            age=age,
	            lambda=lambda,mu=mu,
	            nu=nu, hybprops=hybprops,
	            hyb.rate.fxn=hyb.rate.fxn,
	            hyb.inher.fxn=hyb.inher.fxn,
	            frac=frac,mrca=mrca,
	            complete=complete, stochsampling=stochsampling,
	            trait.model=trait.model)
	out
}

