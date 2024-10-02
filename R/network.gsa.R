#' Sample under the Generalized Sampling Approach
#'
#' @description Takes a phylogeny and samples a period where n lineages exist. This method properly samples n taxa under the GSA
#'
#' @param net A network to sample phylogenies.
#' @param ntaxa The number of desired taxa.
#' @param frac Sampling fraction: The proportion of extant tips included in the phylogeny (incomplete sampling).
#' @param complete If complete = TRUE, the tree with the extinct lineages is returned. If complete = FALSE, the extinct lineages are suppressed.
#' @param stochsampling When `stochsampling=TRUE`: Each extant tip is included into the final tree with probability frac.
#' @return A network with n extant taxa
#'
#' @export
#' @examples
#' set.seed(10)
#' ssa_net <-sim.bdh.taxa.ssa(n=20,numbsim=1,
#' lambda=1,mu=0.2,
#' nu=0.25, hybprops = c(1/3,1/3,1/3),
#' hyb.inher.fxn = make.beta.draw(1,1),
#' )[[1]]
#' gsa_net<-network.gsa(ssa_net,5)
#'
#'
network.gsa <- function(net,ntaxa,complete=TRUE,frac=1,stochsampling=FALSE){
  effective_n<-getEffectiveN(n=ntaxa,frac = frac,stochsampling = stochsampling)

  node_times<-node.depth.edgelength(net)
  ltt<-ltt.network(net,node_times)
  times<-ltt[(ltt[,3])==effective_n,]
  if(nrow(times)==0){
    warning(paste('After accounting for incomplete sampling, we are looking for intervals with ',effective_n, ' lineages. This never occurs on the phylogeny. returning 0',sep=''))
    return(0)
  }
  times <- times[,-3]

  extinct_labels <- get.extinct(net)
  phy<-internal.network.gsa(net = net,times=times,timecreation = node_times,extinct_labels = extinct_labels)

  ##handle complete and sampling frac
  phy<-handleTipsTaxa(phy=phy,complete=complete,target_ntaxa=ntaxa,current_n=effective_n)
  return(phy)
}



internal.network.gsa <- function(net,times,timecreation,extinct_labels){

  ##get weights of each time period with the correct number of tips
  times<-cbind(times,times[,2]-times[,1])
  times[,3]<-times[,3]/sum(times[,3])

  ##pick one of the periods and a time within that period
  time_rw<-sample(x=1:nrow(times),size = 1,prob = times[,3])
  timeslice_time <- runif(1,times[time_rw,1],times[time_rw,2])

  ##slice the network at that time
  net<-internalTimeSliceNetwork(net,timeslice_time,timecreation,extinct_labels)
  return(net)
}
