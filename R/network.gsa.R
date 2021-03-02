network.gsa <- function(net,ntaxa,complete=T,frac=1,stochsampling=F){
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
