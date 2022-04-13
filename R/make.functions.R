#' Make an exponential decay function
#'
#' @description Create an exponential decay function for genetic distance of two taxa and the probability of success of a hybridization event
#' @param t A numeric representing how quickly the hybridization success decays. Samaller values denote a quicker decay
#' @param s A numeric for the power that the genetic distance is raised.
#' @return An exponential decay function
#' @details The function computes: \deqn{e^{- \frac{d^{s}}}{t} }
#' where d is the genetic distance between taxa
#' @export
#'
#' @examples
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
#' @details The function computes: \deqn{1-\frac{d}{t} }
#' where d is the genetic distance between taxa
#' @note a distance \eqn{d} greater than \eqn{t} will return 0
#' @export
#'
#' @examples
make.linear.decay<-function(threshold){
  if( !is.numeric(threshold) ){
    stop("threshold must be a number")
  }
  if(threshold<=0){
    stop("threshold must be positive")
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
#' @details The function computes: \deqn{1- {\frac{d}{t}}^degree}
#' Where d is the distance and t is the threshold
#' @examples
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
make.categorical.draw <- function(inheritances, weights){
  myfunc <-function(){
    return( sample(x=inheritances,size=1,prob = weights ))
  }
  return(myfunc)
}

#' Model for polyploid evolution across the phylogeny
#'
#' @description Create a model that dictates how ploidy evolves and affects the diversification of the phylogeny. This function creates a list that dictates how ploidy affects hybridizations, how ploidy is inherited over time, and how ploidy is inherited across speciation and hybridization events.
#' @return A model for polyploidy evolution to be used as the `polyploid.model` argument in a `sim.bdh function``
#' @param initial_states the initial ploidy state on the phylogeny. if simulating networks with `MRCA=TRUE` then a vector of length two will be required.
#' @param hyb.event.fxn A function that describes how ploidy changes after hybridization events. See Details for more information
#' @param hyb.compatability.fxn A function that describes whether hybridization events can occur between taxa of different ploidy. See Details for more information
#' @param time.fxn A function that describes how ploidy changes over time. See Details for more information
#' @param spec.fxn A function that describes how ploidy changes after speciation events. See Details for more information
#' @details
#' `hyb.event.fxn` is a function that denotes the ploidy of a hybrid child after a hybridization event. The function should have the argument `parent_states`, a vector with the trait states of the two parents to the hybrid child. The function should return a single value for the ploidy state of the hybrid child
#' The default value of `make.allopolyploid.event(1)` causes the ploidy of the hybrid child to be the sum of the ploidy of the two parental lineages.
#'
#' `hyb.compatability.fxn` describes when hybridization events can occur between two taxa based on their ploidy. The function should have the arguments `parent_states` and `inheritance`. `parent_states` is vector with the ploidy states of the hybrid parents while `inheritance` is the inheritance probability of the first lineage denoted in `parent_states`. The function should return `TRUE` for when a hybridization event is allowed to proceed and `FALSE` otherwise.
#' The default value of `make.states.incompatible()` assumes only taxa of the same ploidy can hybridize.
#'
#' `time.fxn` is a function that describes how ploidy changes over time. The function should have the arguments `poly_states` and `timestep` in that order. `poly_states` is a vector containing the ploidy of all taxa while `timestep` is the amount of time given for ploidy evolution. The function should return a vector with the updated ploidy states of all taxa.
#' The default value of `NULL` indicates that ploidy will not evolve within a lineage over time. **NOTE:** ploidy can still change at speciation or hybridization events if allowed.
#'
#' `spec.fxn` is a function that describes how ploidy changes at speciation events. The function should have the argument `tip_state` which has the state of the lineage just before speciation. The function should return a vector with two values, one denoting the ploidy of each of the two new species after the event.
#' The default value of `NULL` causes the two children lineage to inherit the same ploidy as the parental lineage
#'
#' @export
#' @examples
make.polyploid.model <-function(initial_states=1,
                                hyb.event.fxn=make.allopolyploid.event(1),
                                hyb.compatability.fxn=make.states.incompatible(),
                                time.fxn=NULL,
                                spec.fxn=NULL){

  if(is.null(spec.fxn) && is.null(hyb.event.fxn) && is.null(time.fxn) ){
    warning('All components of the polyploid model are NULL. Returning NULL')
    return(NULL)
  }

  ##implement default functions if null
  if(is.null(spec.fxn)){
    spec.fxn<- function(tip_state){
      return(c(tip_state,tip_state))
    }
  }
  if(is.null(time.fxn)){
    time.fxn <-function(poly_states,timestep){
      return(poly_states)
    }
  }

  x<-list(initial=initial_states,
          hyb.event.fxn=hyb.event.fxn,
          hyb.compatability.fxn=hyb.compatability.fxn,
          time.fxn=time.fxn,
          spec.fxn=spec.fxn)
  x[['initial']]<-initial_states
  x[['hyb.event.fxn']]<-hyb.event.fxn
  x[['hyb.compatability.fxn']]<-hyb.compatability.fxn
  x[['time.fxn']]<-time.fxn
  x[['spec.fxn']]<-spec.fxn
  return(x)
}


#' Allopolyploidy at Hybridization Events
#'
#' @description Create a function that creates rules for allopolyploidy at hybridization events. The function creates an argument for the `hyb.event.fxn` element of the `trait.model` argument of `sim.bdh` style functions.
#' @return A function that simulates allopolyploidy at hybridization events
#' @param prob The probability that a hybridization event is an allopolploidy event
#' @export
#' @details
#' If a hybridization event is an allopolyploidy event then the ploidy of the hybrid child is the sum of the ploidy of the two parental lineages. Otherwise, the ploidy of the hybrid child is randomly selected from one of the two parental lineages.
#' @examples
make.allopolyploid.event <-function(prob){
  myfunc<-function(parent_states,inheritance){
    if(runif(1)<=prob){
      return(sum(parent_states))
    }else{
      return(sample(parent_states,size = 1))
    }
  }
  return(myfunc)
}


#' Incompatability at Hybridization Events
#'
#' @description Create a function that makes hybridizations unsuccessul if the two parental taxa don't have the same state. The function creates an argument for the `hyb.compatability.fxn` element of the `trait.model` argument of `sim.bdh` style functions.
#' @return A function that creates a rule for whether hybridization events occur
#' @export
#' @examples
make.states.incompatible <-function(){
  myfunc <- function(parent_states){
    if(parent_states[1]!=parent_states[2]){
      return(F)
    }else{
      return(T)
    }
  }
}

