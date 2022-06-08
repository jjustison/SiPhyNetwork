#' Remove Extinct Lineages from a Phylogenetic Network
#'
#' @description This function removes all extinct tips from a phylogenetic network, returning the reconstructed network.
#'
#' @param net An object of class 'evonet.'
#'
#' @return net The reconstructed network with all extinct tips removed.
#' @export
#'
#' @examples
#' set.seed(17) ##smallest Quartan prime as seed
#' ##Generate a tree with extinct leaves
#' net<-sim.bdh.age(1,1,5,2,1,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),complete=TRUE)[[1]]
#' recon_net<-reconstructedNetwork(net)
#' plot(net)
#' plot(recon_net)
reconstructedNetwork<-function(net){
  extinct_tips<- getReconstructedTips(net$tip.label,get.extinct(net))

  if(length(extinct_tips)==0){
    warning('No extinct tips found. Returning Network unchanged')
    return(net)
  }

  net<-deleteTips(net,extinct_tips)
  return(net)
}

getReconstructedTips<- function(tip.labels,extinct_labels){
  return( which(tip.labels %in% extinct_labels) )
}

#' Sample Tips on a Phylogenetic Network
#'
#' @description This function samples tips from a network. Only extant tips are downsampled from the network. Extinct tips, if they are present, will be unchanged.
#'
#' @param net An object of class 'evonet.'
#' @param rho The sampling probability.
#' @param stochastic If stochastic=FALSE then for a network with n tips we sample n*rho tips. If stochastic=TRUE then each tip porbability rho of being sampled.
#'
#' @return net A network with sampled tips
#' @export
#'
#' @examples
#' set.seed(23) ##set seed with the smallest Pillai prime
#' net<-sim.bdh.age(1,1,3,2,0.125,c(1/3,1/3,1/3),
#' hyb.inher.fxn = make.uniform.draw(),complete = FALSE)[[1]]
#' net<-incompleteSampling(net,0.5,stochastic=FALSE) ##randomly sample half of the extant taxa
incompleteSampling<-function(net,rho,stochastic=FALSE){
  extinct_tips <- getReconstructedTips(1:length(net$tip.label),get.extinct(net))
  if(length(extinct_tips)!=0){
    extant_tips <- (1:length(net$tip.label))[-extinct_tips]
  }else{
    extant_tips <- 1:length(net$tip.label)
  }


  deltips<-getSamplingTips(extant_tips,rho,stochastic)
  if(length(deltips)>0){ ##only delete tips if there are tips to delete
    net<-deleteTips(net,deltips)
  }

  return(net)
}

getEffectiveN <-function(n,frac,stochsampling){
  if(frac!=1){
    if(stochsampling){
      effective_n<-rnbinom(1,size = n,prob = frac)
    } else{
      effective_n<-round(n/frac)
    }
  }else{
    effective_n <- n
  }
  return(effective_n)
}

getSamplingTips <- function(tips,rho,stochastic){
  prop<- 1-rho
  ###Determine which tips need to be deleted###
  if(stochastic){
    n_deltips<-rbinom(1,length(tips),prop)
  }else{
    n_deltips<-round(length(tips)*prop)
  }
  deltips<-sample(tips,n_deltips)
  return(deltips)
}



#' Remove tips from a phylogenetic Network
#'
#' @description This function removes certain tips from a phylogenetic network, returning the pruned network.
#'
#' @param net An object of class 'evonet.'
#' @param tips A numeric vector specifying the tip numbers to delete
#'
#' @return net The network tips removed.
#' @export
#'
#' @examples
#' set.seed(17) ##Set seed with smallest Quartran Prime
#' net<-sim.bdh.age(1,1,5,2,1,c(1/3,1/3,1/3),hyb.inher.fxn = make.uniform.draw(),complete=FALSE)[[1]]
#' net<- deleteTips(net,c(1,6)) ##drop tips 1 and 6
#'
#'
deleteTips<-function(net,tips){

  if( sum(tips>length(net$tip.label))>0 ){
    stop("Some of the tip numbers given are not actually tips.")
  }
  if(length(tips)==length(net$tip.label)){
    warning("removing all tips. Returning NULL")
    return(NULL)
  }
  nodes<-as.list(tips) ##vector of all the nodes we want to delete
  Nnodes<-length(nodes)

  Ntips_original<-length(net$tip.label)
  Ntips_remaining<-Ntips_original-Nnodes
  deleted_nodes<-tips ##keep track of the nodes we will delete
  while(Nnodes!=0){
    # lapply(nodes, function(x) {if(length(x)==0) stop("we have a numeric(0)")})
    nd<-nodes[[Nnodes]] ##the node we are deleting
    nodes[[Nnodes]]<-NULL
    Nnodes<-Nnodes-1

    kid_edges<-which(net$edge[,1]==nd)
    hybkid_edges<-which(net$reticulation[,1]==nd)
    no_kids<-length(kid_edges)
    no_hybkids<-length(hybkid_edges)

    #################
    ###No Kid Case###
    #################
    if( (no_kids+no_hybkids)==0){


      par_e<-which(net$edge[,2]==nd)
      par_nd<-net$edge[par_e,1] ##store the parent node


      ##check if nd is a reticulate node
      hyb_e<-which(net$reticulation[,2]==nd)
      if(length(hyb_e)!=0){##if there is a parental hybrid edge it is a reticulate nd
        hyb_par_nd <- net$reticulation[hyb_e,1] ##The parent of the hybrid edge
        net$reticulation<-subset(net$reticulation,subset= !((1:nrow(net$reticulation))==hyb_e) ) ##delete the hybrid edge
        net$inheritance<-net$inheritance[-hyb_e]
      }

      ##delete nd and par_e. Then update attributes of net
      net$edge<-subset(net$edge,subset= !((1:nrow(net$edge))==par_e) )
      net$edge.length<-net$edge.length[-par_e]

      ###Recursive Calls###
      ##If present, check parent node of the hybrid edge to see if it needs deleting/fusing
      if(length(hyb_e)!=0){
        if(!(hyb_par_nd %in% deleted_nodes)){ ##only add the node if it isn't on the list
          ##Check par_nd to see if that needs deleting/fusing
          nodes[[Nnodes+1]]<-hyb_par_nd
          Nnodes<-Nnodes+1
          deleted_nodes<-c(deleted_nodes,hyb_par_nd)
        }
      }
      if(!(par_nd %in% deleted_nodes)){
        ##Check par_nd to see if that needs deleting/fusing
        nodes[[Nnodes+1]]<-par_nd
        Nnodes<-Nnodes+1
        deleted_nodes<-c(deleted_nodes,par_nd)
      }
    ##################
    ###One Kid Case###
    ##################
    }else if( (no_kids+no_hybkids)==1){
      ##we want to fuse the parent and child edges of nd
      ##determine if the edge to the kid is a hyb edge or normal one

      #####################
      ###Normal Edge Kid###
      #####################
      if(no_kids==1){

        rt<-length(net$tip.label)+1 ##The number of the root
        if(rt!=nd){ ##We can only fuse edges if the nd isn't the root
          par_e<-which(net$edge[,2]==nd)
          par_nd<-net$edge[par_e,1] ##store the parent node
          par_e_length<-net$edge.length[par_e]
          kid_nd<-net$edge[kid_edges,2]

          ##Update Network
          net$edge.length[kid_edges]<-net$edge.length[kid_edges]+par_e_length  ##add the parent edge length to the child edge
          net$edge[kid_edges,1]<-par_nd##have the child edge point from par_nd to kid_nd
          net$edge<-subset(net$edge,subset= !((1:nrow(net$edge))==par_e) ) ##remove  parent edge
          net$edge.length<-net$edge.length[-par_e]
        }else{ ##nd is the root, do nothing
        }
      #####################
      ###Hybrid Edge Kid###
      #####################
      }else if(no_hybkids==1){
        par_e<-which(net$edge[,2]==nd)
        par_nd<-net$edge[par_e,1] ##store the parent node
        kid_nd<-net$reticulation[hybkid_edges,2]

        ##Update Network
        net$reticulation[hybkid_edges,1]<-par_nd ##have the child edge point from par_nd to kid_nd
        net$edge<-subset(net$edge,subset= !((1:nrow(net$edge))==par_e) ) ##remove  parent edge
        net$edge.length<-net$edge.length[-par_e]
      }else{
        stop("Something went awry")
      }
    }else{
      stop(paste("nd had ",(no_kids+no_hybkids)," kids. We should never be at a node with more than 1 kid",sep="") )
    }
  }

  ##adjust node numberings now that we've deleted things
  tip.labels<-rep(NA,Ntips_remaining)
  net$edge<-net$edge*(-1)
  net$reticulation<-net$reticulation*(-1)

  Nnode<-0 ##used to keep track of the internal nodes
  leaf=1 ##Start numbering for leaves
  interior=Ntips_remaining+1 ##Start numberings. start at n+1 because 1:n is reserved for tips
  for (j in sort(unique(as.vector(net$edge)),decreasing = T) ){
    if( (-j) > Ntips_original) {
      replaced_value <- interior
      interior <- interior +1
      Nnode<-Nnode+1
    } else {
      replaced_value <- leaf
      tip.labels[leaf]<-net$tip.label[-j]
      leaf <- leaf +1
    }
    net$reticulation[net$reticulation ==  j]<-replaced_value
    net$edge[net$edge == j]<-replaced_value
  }
  net$Nnode<-Nnode
  net$tip.label<-tip.labels

  ##################################
  ####Hotfix of Github Issue #1 ####
  ##################################
  large_left<-max(net$edge[,1])
  large_right<-max(net$edge[,2])

  if(large_right>large_left){ ##we only care if the largest node number only appears in the second column

    ##do a switcheroo of the node numbering
    net$edge[net$edge==large_left]<- -Inf
    net$edge[net$edge==large_right] <-large_left
    net$edge[net$edge== -Inf] <-large_right

    net$reticulation[net$reticulation==large_left]<- -Inf
    net$reticulation[net$reticulation==large_right] <-large_left
    net$reticulation[net$reticulation== -Inf] <-large_right


  }
  ##################################
  ######## End hotfix Code #########
  ##################################



  return(net)
}


timeSliceNetwork<- function(net,time){
  node_times<-node.depth.edgelength(net)
  net<- internalTimeSliceNetwork(net,time,node_times,extinct_labels = NULL)
  return(net)
}


internalTimeSliceNetwork<-function(net,time,node_times,extinct_labels=NULL){

  Nnode<-max(c(net$edge,net$reticulation))
  nd_rm<- which(time <= node_times) ##nodes to be removed, they exist past the timeslice
  ntips_original<-length(net$tip.label)
  tips<-1:ntips_original
  tips_remaining<- tips[!(tips %in% nd_rm)]

  if(!is.null(extinct_labels)){ ##update the remaining extinct tips if we have that label
    remaining_tip_labels<-net$tip.label[tips_remaining]
    net$extinct<- remaining_tip_labels[remaining_tip_labels %in% extinct_labels]
  }

  ##delete rows that have a start time after the slice time
  keep_edges    <- !(net$edge[,1] %in% nd_rm)
  keep_hyb_edges<- !(net$reticulation[,1] %in% nd_rm)
  net$edge<-subset(net$edge,subset = keep_edges)
  net$edge.length<-net$edge.length[keep_edges]
  net$reticulation<-subset(net$reticulation,subset= keep_hyb_edges)
  net$inheritance<-net$inheritance[keep_hyb_edges]

  ##Deal with edges that cross the timeslice
  cross_edges_ind<-  ( !(net$edge[,1] %in% nd_rm) & (net$edge[,2] %in% nd_rm))
  new_tips<-net$edge[cross_edges_ind,2]
  tips_remaining<-c(tips_remaining,new_tips[new_tips %in% tips])
  new_tips<- new_tips[!(new_tips %in% tips)]
  net$edge.length[cross_edges_ind]<- time- node_times[net$edge[cross_edges_ind,1]]

  ##Deal with the hybrid edges that cross the timeslice
  cross_hyb_edges_ind<-  ( !(net$reticulation[,1] %in% nd_rm) & (net$reticulation[,2] %in% nd_rm))
  if(sum(cross_hyb_edges_ind) > 0){
    cross_hyb_e_lengths <- time - node_times[net$reticulation[cross_hyb_edges_ind,1]]
    new_hyb_tips<- (Nnode+1):(Nnode+sum(cross_hyb_edges_ind))
    new_tips <- c(new_tips,new_hyb_tips)
    net$reticulation[cross_hyb_edges_ind,2] <- new_hyb_tips ##renumber these nodes as they may match some of the tips in new_tips
    net$edge<-rbind(net$edge,net$reticulation[cross_hyb_edges_ind,]) ##add these hyb edges as normal edges
    net$edge.length <- c(net$edge.length,cross_hyb_e_lengths)
    net$reticulation<- subset(net$reticulation,subset = !cross_hyb_edges_ind)
    net$inheritance<-net$inheritance[!cross_hyb_edges_ind]
  }

  ##reassign node numberings
  leaf<-1
  interior<-length(tips_remaining)+length(new_tips)+1
  nds <- sort(unique(c(as.vector(net$edge),as.vector(net$reticulation))))
  new_tip.label<-rep(NA,interior-1)
  new_edge<-net$edge
  new_reticulation<-net$reticulation
  for(i in nds){
    is_remaining_tip<- i %in% tips_remaining
    is_new_tip<- i %in% new_tips
    if(is_remaining_tip || is_new_tip){##This is a tip
      new_edge[net$edge==i]<-leaf
      new_reticulation[net$reticulation==i]<-leaf
      if(is_remaining_tip){
        new_tip.label[leaf]<- net$tip.label[i]
      }else{
        new_tip.label[leaf]<-paste('l',leaf,sep = '')
      }
      leaf<-leaf+1
    } else{
      new_edge[net$edge==i]<-interior
      new_reticulation[net$reticulation==i]<-interior
      interior<-interior+1
    }
  }
  net$Nnode<-interior-leaf
  net$edge<-new_edge
  net$reticulation<-new_reticulation
  net$tip.label<-new_tip.label
  return(net)
}

handleTipsTaxa<-function(phy,complete,target_ntaxa,current_n){

  ##delete nessecary tips
  deltips<-phy$hyb_tips ##These are tips that result form the bdh process and should be fused
  extinct_tips<-getReconstructedTips(phy$tip.label,phy$extinct)
  if(!complete){
    deltips<-c(deltips,extinct_tips) ##add extinct tips if we want the reconstructed tree

  }
  extant_tips<-setdiff(1:length(phy$tip.label),c(phy$hyb_tips,extinct_tips))
  if(current_n!=length(extant_tips)){
    stop("the extant tips don't match the current_n")
  }
  sampling_tips <- sample(x = extant_tips, size = current_n-target_ntaxa)
  deltips<-c(deltips,sampling_tips)
  deltips<-rev(deltips)##we want to get rid of hyb_tips first

  if(length(deltips)==0){
    ##do nothing we don't want to get rid of anything
  }else if(length(deltips)==length(phy$tip.label)){
    #warning('deleting all tips, returning 0')
    return(0)
  }else if(length(deltips) == (length(phy$tip.label)-1) ){
    #warning('Only 1 tip remaining, returning 1')
    return(1)
  }else{
    phy<-deleteTips(phy,deltips)
  }
  phy$hyb_tips <- NULL
  phy$extinct <- NULL



  return(phy)
}



#' Lineages thru time on a network
#'
#' @description This function Computes the number of lineages thru time on a network
#'
#' @param phy An object of class 'evonet.'
#' @param node_times A numeric vector specifying times of each node. If left NULL then the function will use the output from node.depth.edgelength(phy)
#'
#' @return A dataframe that consists of intervals. The first column denotes the start time of the interval while the second column denotes the end time. The third column depicts the number of lineages present in that interval.
#'  NOTE: due to computational precision, two nodes that appear to occur on the same time (as in the case of lineage neutral and generative hybridization) may be part of different intervals in the output data frame.
#' @export
#'
#' @examples
ltt.network<-function(phy,node_times=NULL){
  if(is.null(node_times)){
    node_times<-node.depth.edgelength(phy)
  }

  edges<-rbind(phy$edge,phy$reticulation) ##we want to look at regular edges and hyb edges

  time_intervals<-sort(unique(node_times))
  n_intervals<-length(time_intervals)-1

  start <- time_intervals[1:(n_intervals)]
  end   <- time_intervals[2:(n_intervals+1)]

  intervals<-data.frame(start,end)
  n_lineages<-rep(NA,n_intervals)
  for(i in 1:n_intervals){
    rw<-unlist(intervals[i,])
    ave_time <- mean(rw) ##we get a time between the start and end of the interval. This guarantees we are between events
    x<- sum((node_times[edges[,1]] < ave_time) & (ave_time < node_times[edges[,2]]))

    n_lineages[i] <- x ##the number of lineages is the number of edges that start before the time and end after it
  }
  intervals<-cbind(intervals,n_lineages)
  return(intervals)
}

#' Create a more Plotting-Friendly phylogenetic Network
#'
#' @description This function creates a more plotting-friendly `evonet` object to be used with the `plot()` function
#'
#' @param net An object of class `evonet`.
#' @param tol a tolerance to account for floating point imprecision. any values smaller than `tol` are considered to be zero.

#' @return a network to be used with the `plot()` function
#' @export
#'
#' @examples
plottable.net<-function(net,tol=1e-8){

  node_times<-node.depth.edgelength(net)

  nhybs<-nrow(net$reticulation)

  nd_num<- max(net$edge)
  bad_hybs<-c()
  new_elength<-c()
  for(i in 1:nhybs){
    rw<-net$reticulation[i,]
    from_time<-node_times[rw[1]]
    to_time<-node_times[rw[2]]
    if(abs(from_time-to_time)>tol){
      bad_hybs<-c(bad_hybs,i)
      new_elength<-c(new_elength,(to_time-from_time))
    }
  }
  if(length(bad_hybs)>0){

    new_edges<-matrix(nrow = length(bad_hybs),ncol=2)

    for(i in 1:length(bad_hybs)){
      bad_e<-bad_hybs[i]
      new_edges[i,]<-net$reticulation[bad_e,]
      new_edges[i,2]<-nd_num+i
      net$reticulation[bad_e,1]<-nd_num+i
    }

    net$edge<-rbind(net$edge,new_edges)
    net$edge.length<-c(net$edge.length,new_elength)

    ##adjust node numberings now that we've deleted things
    tip.labels<-rep(NA,length(net$tip.label)+length(bad_hybs))
    orig_tips<-(1:length(net$tip.label))
    new_tips <-nd_num+(1:length(bad_hybs))
    tips<-c(orig_tips,new_tips)
    net$edge<-net$edge*(-1)
    net$reticulation<-net$reticulation*(-1)

    leaf=1 ##Start numbering for leaves
    interior=length(tip.labels)+1 ##Start numberings. start at n+1 because 1:n is reserved for tips
    for( j in (unique(as.vector(net$edge)))  ){
      if( -j %in% tips ) {
        replaced_value <- leaf
        if(-j %in% orig_tips){
          tip.labels[leaf]<-net$tip.label[-j]
        }else{
          tip.labels[leaf]<- ''
        }
        leaf <- leaf +1
      } else {
        replaced_value <- interior
        interior <- interior +1
      }
      net$reticulation[net$reticulation ==  j]<-replaced_value
      net$edge[net$edge == j]<-replaced_value
    }
    net$Nnode<-interior-leaf
    net$tip.label<-tip.labels

  }
  return(net)
}

get.extinct<-function(net,tol=1e-8){
  nd_times<-node.depth.edgelength(net)
  extant_time<-max(nd_times)
  tips<-1:length(net$tip.label)
  tip_times<-nd_times[tips]
  extinct_tips<-which(abs(extant_time-tip_times)>tol)
  extinct_labels<-net$tip.label[extinct_tips]
  return(extinct_labels)
}



