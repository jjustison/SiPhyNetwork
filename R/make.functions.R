#' Make an exponential decay function
#'
#' @description Create an exponential decay function for genetic distance of two taxa and the probability of success of a hybridization event
#' @param t A numeric representing how quickly the hybridization success decays. Samaller values denote a quicker decay
#' @param s A numeric for the power that the genetic distance is raised.
#' @return An exponential decay function
#' @details The function computes: \deqn{e^{- \frac{d^{s}}{t}}}{ exp(-(d^s)/t) }
#' where d is the genetic distance between taxa
#' @export
#'
#' @examples
#' set.seed(17)
#' dist_func<- make.exp.decay(1,1)
#' net<-sim.bdh.age(1,1,5,2,2,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),
#' hyb.rate.fxn=dist_func,complete=TRUE)[[1]]
make.exp.decay<-function(t=1,s=1){
  if( !(is.numeric(t) && is.numeric(s)) ){
    stop("t and s both must be numbers")
  }
  if(t<=0){
    stop("t must be positive")
  }
  myfunction<-function(distance){
    return(exp( -((distance^s)/t) ))
  }
  return(myfunction)
}

#' Make a linear decay function
#'
#' @description Create a linear decay function for genetic distance of two taxa and the probability of success of a hybridization event
#' @param threshold A numeric representing how quickly the hybridization success decays. Smaller values denote a quicker decay
#' @return A linear decay function
#' @details The function computes: \deqn{1-\frac{d}{t} }{1-(d/t)}
#' where d is the genetic distance between taxa
#' @note a distance \eqn{d} greater than \eqn{t} will return 0
#' @export
#'
#' @examples
#' set.seed(17)
#' dist_func<- make.linear.decay(0.5)
#' net<-sim.bdh.age(1,1,5,2,2,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),
#' hyb.rate.fxn=dist_func,complete=TRUE)[[1]]
make.linear.decay<-function(threshold){
  if( !is.numeric(threshold) ){
    stop("threshold must be a number")
  }
  if(threshold<0){
    stop("threshold must be positive")
  }
  if(threshold==0){
    return(NULL)
  }
  myfunction<-function(distance){
    mult<-(distance<threshold)
    val<-mult*(1-(distance/threshold))
    return(val)
  }
}

#' Make a stepwise decay function
#'
#' @description Create a stepwise decay function for genetic distance of two taxa and the probability of success of a hybridization event
#' @param probs Vector of dimension k, where k is the number of different probabilities of success. An individual time between (distances\[i-1],distances\[i]] has probability of success prob\[i]
#' @param distances Vector of k, containing the end of each interval where success probabilities shift. The first interval where success is prob\[1] is [0,distances\[1]]. For all i>1, the probability of success is prob\[i] over the interval (distance\[i-1],distances\[i]].
#' @return An stepwise decay function
#' @export
#'
#' @examples
#' set.seed(17)
#' dist_func<- make.stepwise(probs=c(1,0.5,0),distances=c(0.25,0.5,Inf))
#' net<-sim.bdh.age(1,1,5,2,2,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),
#' hyb.rate.fxn=dist_func,complete=TRUE)[[1]]
make.stepwise<-function(probs,distances){
  if( !(is.numeric(probs) && is.numeric(distances))){
    stop("probs and distances both must be vectors of numbers")
  }
  if(length(probs) != length(distances) ){
    stop("the length of the probs vector should be equal to the length of the distance vector")
  }
  if(sum(probs<0)!=0){
    stop("probs should all be nonnegative")
  }
  if(sum(distances<=0)!=0){
    stop("Distances should all be positive")
  }

  myfunction<-function(distance){
    x<-which(distances>=distance)
    if(length(x)==0){
      stop(paste("Rates are defined on the range 0 to ",max(distances),". You gave a distance of ",distance,sep = "" ))
    }
    return(probs[ min(which(distances>=distance))])
  }
  return(myfunction)
}

#' Make a polynomial decay function
#'
#' @description Create a polynomial decay function for genetic distance of two taxa and the probability of success of a hybridization event
#' @param threshold A numeric denoting how quickly the polynomial function decays. Distances greater than the threshold will return a success probability of 0.
#' @param degree The degree of the polynomial
#' @return An polynomial decay function
#' @export
#' @details The function computes: \deqn{1- {\frac{d}{t}}^degree}{1-((d/t)^degree)}
#' Where d is the distance and t is the threshold
#' @examples
#' set.seed(17)
#' dist_func<- make.polynomial.decay(0.5,2)
#' net<-sim.bdh.age(1,1,5,2,2,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),
#' hyb.rate.fxn=dist_func,complete=TRUE)[[1]]
make.polynomial.decay<-function(threshold,degree=1){
  if( !is.numeric(threshold) ){
    stop("threshold must be a number")
  }
  if(threshold<=0){
    stop("threshold must be positive")
  }
  if(degree<=0){
    stop("degree must be positive")
  }
  myfunction<-function(distance){
    mult<-(distance<threshold)
    val<- mult*(1- ((distance/threshold)^degree))
    return(val)
  }
  return(myfunction)
}

#' Function that makes a draw from a uniform distribution
#'
#' @description Create a function that makes draws from a uniform distribution. Function calls `runif(1)`
#' @return A function that makes draws from a uniform distribution
#' @export
#' @examples
#' set.seed(17)
#' net<-sim.bdh.age(1,1,5,2,1,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),
#' complete=TRUE)[[1]]
make.uniform.draw <-function(){
  myfunc <- function(){
    return(runif(1))
  }
  return(myfunc)
}

#' Function that makes a draw from a beta distribution
#'
#' @description Create a function that makes draws from a beta distribution. Function calls `rbeta(1,alpha,beta)`
#' @return A function that makes draws from a beta distribution.
#' @param alpha The first shape parameter of the beta distribution
#' @param beta The second shape parameter of the beta distribution
#' @export
#' @examples
#' set.seed(17)
#' inher_func<-make.beta.draw(10,10)
#' net<-sim.bdh.age(1,1,5,2,1,c(1/3,1/3,1/3),hyb.inher.fxn = inher_func,
#' complete=TRUE)[[1]]
#'
make.beta.draw <-function(alpha,beta){
  myfunc <- function(){
    return(rbeta(1,shape1=alpha,shape2=beta))
  }
  return(myfunc)
}

#' Function that makes a draw from a categorical distribution
#'
#' @description Create a function that makes draws from a categorical distribution. Function calls `sample(x=inheritances,size=1,prob=weights)`
#' @return A function that makes draws from a categorical distribution.
#' @param inheritances A vector of inheritance probabilities
#' @param weights A vector of weights for each inheritance probability
#' @export
#' @examples
#' set.seed(17)
#' inher_func<-make.categorical.draw(inheritances=c(0.25,0.50,0.75),weights=c(0.25,0.5,0.25))
#' net<-sim.bdh.age(1,1,5,2,1,c(1/3,1/3,1/3),hyb.inher.fxn = inher_func,
#' complete=TRUE)[[1]]
make.categorical.draw <- function(inheritances, weights){
  myfunc <-function(){
    return( sample(x=inheritances,size=1,prob = weights ))
  }
  return(myfunc)
}

#' Model for trait evolution across the phylogeny
#'
#' @description Create a model that dictates how a discrete or continuous trait evolves and affects the diversification of the phylogeny. This function creates a list that dictates how the trait affects hybridizations, how the trait is changes over time, and how the trait is inherited across speciation and hybridization events.
#' @return A model for trait evolution to be used as the `trait.model` argument in a `sim.bdh function``
#' @param initial_states the initial state on the phylogeny. if simulating networks with `twolineages=TRUE` then a vector of length two will be required.
#' @param hyb.event.fxn A function that describes how the trait changes after hybridization events. See Details for more information
#' @param hyb.compatibility.fxn A function that describes whether hybridization events can occur between taxa based on their trait values. See Details for more information
#' @param time.fxn A function that describes how trait values changes over time. See Details for more information
#' @param spec.fxn A function that describes how trait values change at speciation events. See Details for more information
#' @details
#' `hyb.event.fxn` is a function that denotes the trait value of a hybrid child after a hybridization event. The function should have the argument `parent_states`, a vector with the trait states of the two parents to the hybrid child and `inheritance`. `parent_states` is vector with the states of the hybrid parents while `inheritance` is the inheritance probability of the first lineage denoted in `parent_states`. The function should return a single value for the trait state of the hybrid child.
#'
#' `hyb.compatibility.fxn` describes when hybridization events can occur between two taxa based on their trait values. The function should have the arguments `parent_states`.  The function should return `TRUE` for when a hybridization event is allowed to proceed and `FALSE` otherwise.
#'
#' `time.fxn` is a function that describes how trait values change over time. The function should have the arguments `trait_states` and `timestep` in that order. `trait_states` is a vector containing the ploidy of all taxa while `timestep` is the amount of time given for ploidy evolution. The function should return a vector with the updated ploidy states of all taxa.
#' The default value of `NULL` indicates that trait values will not evolve within a lineage over time. **NOTE:** Values can still change at speciation or hybridization events if allowed.
#'
#' `spec.fxn` is a function that describes how trait values change at speciation events. The function should have the argument `tip_state` which has the state of the lineage just before speciation. The function should return a vector with two values, one denoting the trait of each of the two new species after the event.
#' The default value of `NULL` causes the two children lineage to inherit the same trait value as the parental lineage
#'
#' @export
#' @examples
#' initial_val<-2 ## The root starts off at 2N
#'
#' ###function for what happens at hybridization event
#' hyb_e_fxn <- function(parent_states,inheritance){
#'  ##For allopolyploidy we add the ploidy of both parents
#'  return(sum(parent_states))
#'}
#'
#' ##Function for determining whether hybridization occurs
#'hyb_c_fxn <-function(parent_states,hybrid_state){
#'  ##Hybridization occurs only when the ploidy is the same
#'  return(parent_states[1]==parent_states[2])
#'}
#'
#'
#'##Function for how the trait changes over time
#'t_fxn <- function(trait_states,timestep){
#'  ##We assume that autopolyploidy occur exponentially with rate lambda
#'  lambda<- 2 ##Rate of autopolyploidy
#'
#'  ##The number of autopolyploidy events that occur on each lineage over the timestep
#'  nevents<-rpois(length(trait_states),timestep)
#'
#'  ##each event doubles the ploidy
#'  new_states<- trait_states * (2^nevents)
#'  return(new_states)
#'}
#'
#'##Function for how the trait changes at speciation events
#'s_fxn <-function(tip_state){
#'  ##Ploidy doesn't change at speciation events.
#'  ##Both daughter lineages have the same ploidy as the parent
#'  return(c(tip_state,tip_state))
#'}
#'
#'trait_model<-make.trait.model(initial_states = initial_val,
#'                              hyb.event.fxn = hyb_e_fxn,
#'                              hyb.compatibility.fxn = hyb_c_fxn,
#'                              time.fxn = t_fxn,
#'                              spec.fxn = s_fxn)
#'

make.trait.model <-function(initial_states,
                                hyb.event.fxn,
                                hyb.compatibility.fxn,
                                time.fxn=NULL,
                                spec.fxn=NULL){

  ##implement default functions if null
  if(is.null(spec.fxn)){
    spec.fxn<- function(tip_state){
      return(c(tip_state,tip_state))
    }
  }
  if(is.null(time.fxn)){
    time.fxn <-function(trait_states,timestep){
      return(trait_states)
    }
  }

  x<-list(initial=initial_states,
          hyb.event.fxn=hyb.event.fxn,
          hyb.compatibility.fxn=hyb.compatibility.fxn,
          time.fxn=time.fxn,
          spec.fxn=spec.fxn)
  x[['initial']]<-initial_states
  x[['hyb.event.fxn']]<-hyb.event.fxn
  x[['hyb.compatibility.fxn']]<-hyb.compatibility.fxn
  x[['time.fxn']]<-time.fxn
  x[['spec.fxn']]<-spec.fxn
  return(x)
}



