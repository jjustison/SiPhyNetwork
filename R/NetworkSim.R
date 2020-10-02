##This is where we simulate the lineage neutral hybridization events on the phylogeny
sim.hyb.neutral<-function(phy,total_time,events,nu,alpha,beta){
  ##the dynamics of gene flow are dependant on the number of lineages.
  ##Thus we will go thru the intervals denoted in events and simulate gene flow between each speciation/extinction event

  ##TODO be smarter about indexing and where things are in the loop so I don't have to add an arbitrary row


  tot_time<- total_time ##the total length of the tree
  curr_time<-0
  rw<-1 ##row index of where we currently are along the ltt coordinates
  start_t<-events$times[rw]
  if(nrow(events)!=1){
    end_t<-events$times[rw+1]
  } else{ ##If there is only one event then we take the end time to be total_time
    end_t<-total_time
  }

  dur<-end_t-start_t ##events occur when rexp() is less than dur
  n_lin<-events$N[rw] ##number of lineges at current time

  ###this is where we store the details of hybridization events
  from<-numeric()
  to<-numeric()
  time_event<-numeric()
  prop<-numeric()
  ##TODO if building these vectors becomes prohibitively slow  I could go thru each interval and draw from rpois() to get the number of events in each interval and initialize a dataframe with that many rows

  while(rw<=nrow(events)){
    if(n_lin==1){##Check the number of lineages
      rand<-Inf ##We can't have a hybridization with only one lineage
    }else{
      binom_coef<-choose(n_lin,2)
      rand<-rexp(1,(binom_coef*nu)) ## the rate of a gene flow event is (choose(n,2)*nu) for the interval
    }

    if(rand<dur){##an event occurs...

      ##First we pick two lineages for gene flow - we'll call these donor lineages
      edges<-edges.in.time(phy,curr_time) ##These are the edges that are 'extant' at our current time
      donors<-sample(edges,size=2) ##pick two of these lineages

      curr_time<-curr_time+rand ##update time
      dur<-end_t-curr_time ##update duration

      ##record the event
      from<-c(from,donors[1])
      to<-c(to,donors[2])
      time_event<-c(time_event,curr_time)
      prop<-c(prop,rbeta(1,alpha,beta)) ##this could be done afterwards as it doesn't depend on anything. Could speed things up

    } else{ ##an event doesn't occur in this interval
      ##move on to the next interval
      rw<-rw+1

      ##Set up the curr_time and interval time
      if(rw<=nrow(events)){
        start_t<-events$times[rw]
        if(rw!=nrow(events)){
          end_t<-events$times[rw+1]
        } else{
          end_t<- tot_time
        }
        dur<-end_t-start_t ##events occur when rexp() is less than dur
        n_lin<-events$N ##number of lineges at current time
        curr_time<-start_t
        n_lin<-max(events$N[which(events$times==curr_time)]) ##This gets around the odd case where two events happen at the same time
      }
    }
  } ##end simulating events

  ##Make a from with all the event information
  hyb_frame<-data.frame(from,to,time_event,prop)

  ##From and To are currently edges. we will need to make nodes at these places and we will store the nodes here
  from_nd<-rep(-Inf,length(from))
  to_nd<-rep(-Inf,length(to))
  hyb<-data.frame(from_nd,to_nd,time_event,prop) ##TODO could probably be smarter than making a whole new dataframe


  if(nrow(hyb_frame)!=0){##only do stuff if there are events

    node_times<-node.depth.edgelength(phy) ##TODO this is redundant given 'events'?
    for( i in 1:nrow(hyb_frame)){ #go thru each event and make appropriate nodes at each one
      edges<-hyb_frame[i,c(1,2)] ##these are the edges we will need to put nodes in.
      ord<-order(edges,decreasing=T) #to avoid indexing issues we do the later edge first

      for(j in ord){
        e<-unlist(edges[j])

        ##get the time that we make the node. it should be (time_event-time_parent_node)
        par_node<-phy$edge[e,1]
        child_node<-phy$edge[e,2]
        par_node_time<-node_times[par_node]
        t_event<-hyb_frame[i,3]-par_node_time

        items<-make.node(phy,e,t_event)
        phy<-items[[1]] ##record changes in the tree

        ##TODO use doBookKeeping() to do this bookkeeping
        ##update the nodes in hyb since we added a node
        hyb[hyb[,1]>=items[[2]],1]<-hyb[hyb[,1]>=items[[2]],1]+1
        hyb[hyb[,2]>=items[[2]],2]<-hyb[hyb[,2]>=items[[2]],2]+1


        ##we need to update the ev_node in events since we added a node
        events$ev_node[items[[2]]<=events$ev_node]<-events$ev_node[items[[2]]<=events$ev_node]+1

        hyb[i,j]<-items[[2]] ##record the node we just added

        ##update node_times since we just added a node
        node_times<-append(node_times,hyb_frame[i,3],par_node)

        ##TODO use doBookKeeping() to do this bookkeeping
        ##we need to update the edges in hyb_frame since we added an edge
        hyb_frame[e<=hyb_frame[,1],1]<-hyb_frame[e<=hyb_frame[,1],1]+1
        hyb_frame[e<=hyb_frame[,2],2]<-hyb_frame[e<=hyb_frame[,2],2]+1
      }
    }
  }
  return( list(phy,hyb,events) )
}



##returns the index of edges that cross thru time t
edges.in.time<-function(phy,t){
  ##NOTE: t is in terms of time from the root

  bts<-node.depth.edgelength(phy) ##get the timings of each node
  e<-which( (bts[phy$edge[,1]]<=t)  &  (t<bts[phy$edge[,2]]) ) ## get edges with parent nodes before t and child nodes after t
  return(e)
}


##makes a new node by splitting an edge into two separate ones with a node between
##TODO do bookkeeping within this function? make a VERBOSE option?
make.node<-function(phy,e,time){

  if( time<0 | phy$edge.length[e]<time  ){
    stop("YOU GOOFBALL! The time needs to be within the length of the edge e")
  }

  nds<-phy$edge[e,] ##the parent and child nodes respectively of edge e
  par_nd<-nds[1]

  ##TODO use doBookKeeping() to do this bookkeeping
  ##Increase all nodes greater than par_nd by 1
  phy$edge[phy$edge[,1]>par_nd,1]<-phy$edge[phy$edge[,1]>par_nd,1]+1
  phy$edge[phy$edge[,2]>par_nd,2]<-phy$edge[phy$edge[,2]>par_nd,2]+1

  new_nd<-par_nd+1

  phy$edge<-insertRow(phy$edge,c(new_nd,phy$edge[e,2]),e+1)   ##add new edge in
  phy$edge[e,2]<-new_nd  ##adjust current edge to point to new_nd instead of child_nd

  ##adjust edgelength and add new edgelength in
  phy$edge.length<-append(phy$edge.length,phy$edge.length[e]-time,after=e)
  phy$edge.length[e]<-time

  phy$Nnode<-phy$Nnode+1 ##we added a node
  return(list(phy,new_nd)) ##return the altered phylogeny and the location of the new node
}

##I obtained this function from Roland at stackoverflow: https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended
insertRow <- function(existingDF, newrow, r) {
  existingDF <- rbind(existingDF,newrow)
  existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
  rownames(existingDF)<-c()
  return(existingDF)
}


reduce.reuse.recycle<-function(phy,hyb){

  ##first get the reticulate nodes
  ret_nds<-hyb$to_nd

  ##get the parent nodes of the reticulate nodes
  nds<-c(phy$edge[(phy$edge[,2] %in% ret_nds),1],hyb$from_nd)
  nds<-nds[!(nds %in% ret_nds)] ##we don't want to consider reticulate nodes. i.e. we don't want nodes with an in degree >1

  ##We want to check the out degree of each of these nodes
  nd_degree<-rep(NA,length(nds))
  for(i in 1:length(nds)){
    nd<-nds[i]
    ##check the out-ness degree in phy$edge and add the out-ness from hyb$from
    nd_degree[i]<-sum(phy$edge[,1]==nd) + sum(hyb$from_nd==nd)
  }

  nds<-nds[which(nd_degree==1)] ##We will reduce, reuse, and recycle the nodes with an out degree of 1
  for(i in seq_along(nds)){
    nd<-nds[i]

    par_edge<-which(phy$edge[,2]==nd) ##TODO this is somewhat inefficient as I'm doing something similar above
    child_edge<-which(hyb$from_nd==nd)
    par_nd<-phy$edge[par_edge,1] ## the ancestral node of par_edge

    ###########
    ##REDUCE###
    ###########
    phy$edge<-phy$edge[-par_edge,] ##delete par_edge
    phy$edge.length<-phy$edge.length[-par_edge]
    phy$tip.label<-phy$tip.label[-nd]


    ###########
    ###REUSE###
    ###########
    hyb$from_nd[child_edge]<-par_nd ##remap the child edge so we have par_nd going to child_node

    ###########
    ##RECYCLE##
    ###########
    ##TODO use doBookKeeping() to do this bookkeeping
    ##update phy since we effectively deleted par_edge and nd
    phy$edge[(phy$edge[,1]>nd),1]<-phy$edge[(phy$edge[,1]>nd),1] -1
    phy$edge[(phy$edge[,2]>nd),2]<-phy$edge[(phy$edge[,2]>nd),2] -1


    ##update nds
    nds[nds>=nd]<-nds[nds>=nd]-1

    ##update hyb
    hyb$from_nd[hyb$from_nd>nd]<-hyb$from_nd[hyb$from_nd>nd]-1
    hyb$to_nd[hyb$to_nd>nd]<-hyb$to_nd[hyb$to_nd>nd]-1
  }
  return(list(phy,hyb))
}

determineNetwork<-function(phy,n,lambda,mu,nu,hybprops,alpha,beta,frac=1,complete=TRUE,stochsampling=FALSE){
  ##Construct a dataframe that stores the speciation, extinction, and (de)generative hybridization events
  ##The data frame stores:
  ##The time of the event
  ##number of lineages after the event
  ##The type of event - Speciation (S), Extinction (E), Hybrid Speciation (HS), Hybrid Extinction (HE)
  ##Node responsible for the respective event

  times<-node.depth.edgelength(phy) ##time of events
  tot_time<-max(times)
  Nevents<-length(times)-n
  ev_node<-order(times)[1:Nevents] ##The node responsible for the event
  times<-times[ev_node] ##get non-tip times
  if(length(times)!=length(unique(times)) ){ ##this should happen with Pr()=0 but you never know with computer precision
    warning("Two events happen at the same time. Things may not work properly")
  }
  ev_t<-(ev_node %in% phy$edge[,1]) ##the type of event. will be made more explicit later

  N<-ev_t ##This will store the number of lineages after an event
  N[!N]<- -1 ##speciations add 1 to the number of lineages while extinctions take one away
  N[1]<- 2 ##start with two lineages
  N<-cumsum(N)

  ##Now we need to determine whether a given speciation or extinction event is associated with a hybridization
  hyb_lambda<-  hybprops[1]*nu
  hyb_mu<-      hybprops[2]*nu
  ev_type<-rep(NA,Nevents) ##This will store the more explicit type of event
  ev_type[1]<-"S" ##First event is always a regular speciation
  if(length(ev_t)>1){
    for(rw in 2:Nevents){
      rand<-runif(1,0,1)
      N0<-N[rw-1]
      if(ev_t[rw]){ ##Speciation
        if( rand <= ( (choose(N0,2)*hyb_lambda)/( (choose(N0,2)*hyb_lambda)+(N0*lambda)) ) ){ #hybrid speciation
          ev_type[rw]<-"HS"
        }else{ #regular speciation
          ev_type[rw]<-"S"
        }
      }else{ ##Extinction
        if( rand <= ( (choose(N0,2)*hyb_mu)/( (choose(N0,2)*hyb_mu)+(N0*mu)) ) ){ #hybrid extinction
          ev_type[rw]<-"HE"
        }else{ #regular extinction
          ev_type[rw]<-"E"
        }
      }
    }
  }
  events<-data.frame(times,N,ev_type,ev_node)

  ##still debating whether I should subset events prior to the for loop
  #events<-events[(events$ev_type=="HS" | events$ev_type=="HE"),]## subset the dataframe to "HS" or "HE"
  #hyb<-data.frame(from_nd=rep(-Inf,nrow(events)),to_nd=rep(-Inf,nrow(events)),time_event=events$time,prop=rbeta(nrow(events),alpha,beta)) ##TODO could probably be smarter than making a whole new dataframe
  Nevents<-sum(events$ev_type=="HS" | events$ev_type=="HE") ##number of hybridization events
  new_hybs<-data.frame(from_nd=rep(-Inf,Nevents),to_nd=rep(-Inf,Nevents),time_event=rep(-Inf,Nevents),prop=rbeta(Nevents,alpha,beta))


  ##simulate lineage neutral hybrid events
  if(nu*hybprops[3]>0){ ##TODO this type of checking and case handling within sim.hyb.neutral
    items<-sim.hyb.neutral(phy,tot_time,events,nu*hybprops[3],alpha,beta)
    phy<-items[[1]]
    hyb<-items[[2]]
    events<-items[[3]]
    hyb_rw<-nrow(hyb)+1 ##Row of hyb where the new hybridization events start

    if(nrow(hyb)!=0){
      if(nrow(new_hybs)!=0){
        hyb<-insertRow(hyb,new_hybs,(nrow(hyb)+nrow(new_hybs)))
      }
    } else{
      hyb<-new_hybs
    }
  } else{ ##skip the function if the rate is 0
    hyb<-new_hybs
    hyb_rw<-1 ##Row of hyb where the new hybridization events start
  }

  ##deal with (de)generative hybridization events
  ##We Add the appriopiate nodes to for each event and record where hybrid edges get added. Things are stored in 'hyb'
  node_times<-node.depth.edgelength(phy) ##TODO. This vector is a hacky fix. This information is mostly in 'events' and leads to somewhat redundant computation and bookkeeping
  if(length(events)>1){
    for(rw in 2:nrow(events)){
      #for(rw in 2:4){
      if(events$ev_type[rw] %in% c("HS","HE")){
        t<-events$time[rw] ##time of event
        nd<-events$ev_node[rw] ##node of the event

        if(events$ev_type[rw]=="HS"){ ##speciating Hybridizaiton
          event_edge<-which(phy$edge[,2]==nd) ##the edge that leads to the event node

          ##This needs rethinking
          # t_before_edge<-events$times[rw-1]
          # t_before_hyb<-hyb$time_event[hyb$time_event<t]
          # t_before<-(max(t_before_hyb,t_before_edge)+t)/2 ##merely max() should work but we take the average for good measure
          t_before<-(max(node_times[node_times<t])+t)/2 ## 'node_times[rw-2]' should be sufficient but we take a time between two events for good measure

          edges_t<-edges.in.time(phy,t_before )
          edges_t<-edges_t[!(edges_t %in% event_edge)] ##we exlude the edge that leads to the event

          if(length(edges_t)==1){
            to_edge<-edges_t
          } else{
            to_edge<-sample(x=edges_t,1) ## pick an edge to make a node
          }
          par_nd<-phy$edge[to_edge,1] ##the parent node of to_edge
          items<-make.node(phy=phy,e=to_edge,t-node_times[par_nd])
          phy<-items[[1]]
          ##TODO use doBookKeeping() to do this bookkeeping
          ##update the nodes in hyb since we added a node
          hyb[hyb[,1]>=items[[2]],1]<-hyb[hyb[,1]>=items[[2]],1]+1
          hyb[hyb[,2]>=items[[2]],2]<-hyb[hyb[,2]>=items[[2]],2]+1

          ##we need to update the ev_node in events since we added a node
          events$ev_node[items[[2]]<=events$ev_node]<-events$ev_node[items[[2]]<=events$ev_node]+1
          nd<-events$ev_node[rw] ##node of the event
          events$ev_node[rw]<-items[[2]] ##the node we just created is the 'from' node

          ##update node_times since we just added a node
          node_times<-append(node_times,t,par_nd)

          edges_children<-which(phy$edge[,1]==nd) ##the edges that are children to nd
          if(length(edges_children)==1){ ## pick an edge to make a node
            e<-edges_children
            error("This shouldn't happen(?)")
          } else{
            e<-sample(x=edges_children,1)
          }

          par_nd<-phy$edge[e,1]
          items<-make.node(phy,e=e,0)

        } else if(events$ev_type[rw]=="HE"){ ##Lineage degenerative hybridization
          ##find edges that are present after the extinction
          edges_t<-edges.in.time(phy,t)
          if(length(edges_t)==1){ ## pick an edge to make a node
            from_edge<-edges_t
          } else{
            from_edge<-sample(x=edges_t,1)
          }
          par_nd<-phy$edge[from_edge,1]
          items<-make.node(phy=phy,e=from_edge,t-node_times[par_nd])
        }
        phy<- items[[1]]
        ##TODO use doBookKeeping() to do this bookkeeping
        ##update the nodes in hyb since we added a node
        hyb[hyb[,1]>=items[[2]],1]<-hyb[hyb[,1]>=items[[2]],1]+1
        hyb[hyb[,2]>=items[[2]],2]<-hyb[hyb[,2]>=items[[2]],2]+1

        ##we need to update the ev_node in events since we added a node
        events$ev_node[items[[2]]<=events$ev_node]<-events$ev_node[items[[2]]<=events$ev_node]+1

        ##update node_times since we just added a node
        node_times<-append(node_times,t,par_nd)

        ##record the gene flow event in hyb
        hyb[hyb_rw,1:3]<-c(events$ev_node[rw],items[[2]],events$time[rw])
        hyb_rw<-hyb_rw+1
      }
    }
  }

  x<-reduce.reuse.recycle(phy,hyb)
  netphy<-x[[1]]
  nhyb<-x[[2]]

  net<-evonet(netphy,nhyb$from_nd,nhyb$to_nd)
  net$inheritance<- hyb$prop

  ##Deal with net reconstruction and incomplete sampling
  if(!complete){
    net<-reconstructNetwork(net)
  }
  net<-incompleteSampling(net,rho=frac,stochastic = stochsampling)

  return(net)
}


#' Remove Extinxt Lineages from a Phylogenetic Tree
#'
#' @description This function removes all extinct tips from a phylogenetic network, returning the reconstructed network.
#'
#' @param net An object of class 'evonet.'
#'
#' @return net The reconstructed network with all extinct tips removed.
#' @export
#'
#' @examples
reconstructNetwork<-function(net){
  extinct_tips<-which( net$tip.label %in% getExtinct(net))
  tracked<-list(extinct=extinct_tips)

  items<-list(net=net,tracked=tracked)
  for(i in 1:length(tracked$extinct)){
    items<-delNode(items$net,items$tracked$extinct[i],tracked=items$tracked)
  }

  if( class(items$net$edge)!="matrix" ){
    items$net$edge<-matrix(items$net$edge,ncol=2,byrow=TRUE)
    colnames(items$net$edge)<-NULL
  }
  return(items$net)
}

#' Sample Tips on a Phylogenetic Network
#'
#' @description This function samples tips from a network. Both extant and extinct tips are sampled from the network.
#'
#' @param net An object of class 'evonet.'
#' @param rho The sampling probability.
#' @param stochastic If stochastic=FALSE then for a network with n tips we sample n*rho tips. If stochastic=TRUE then each tip porbability rho of being sampled.
#'
#' @return net A network with sampled tips
#' @export
#'
#' @examples
incompleteSampling<-function(net,rho,stochastic=F){
  prop<- 1-rho
  ###Determine which tips need to be deleted###
  if(stochastic){
    n_deltips<-rbinom(1,length(net$tip.label),prop)
  }else{
    n_deltips<-round(length(net$tip.label)*prop)
  }
  deltips<-sample.int(length(net$tip.label),n_deltips)

  ###Delete specified tips###
  tracked<-list(deltips=deltips)
  items<-list(net=net,tracked=tracked)

  if(length(tracked$deltips)>0){ ##only delete tips if there are tips to delete
    for(i in 1:length(tracked$deltips)){
      items<-delNode(items$net,items$tracked$deltips[i],tracked=items$tracked)
    }
  }

  if( class(items$net$edge)!="matrix" ){ ##if only one tip left then net$edge becomes a vector
    items$net$edge<-matrix(items$net$edge,ncol=2,byrow=TRUE) ##turn it back into a data.frame
    colnames(items$net$edge)<-NULL
  }
  return(items$net)
}

delNode<-function(net,nd,tracked=list() ){
  kid_edges<-which(net$edge[,1]==nd)
  hybkid_edges<-which(net$reticulation[,1]==nd)
  no_kids<-length(kid_edges)
  no_hybkids<-length(hybkid_edges)

  #################
  ###No Kid Case###
  #################
  if( (no_kids+no_hybkids)==0){

    par_e<-which(net$edge[,2]==nd)
    tracked[[length(tracked)+1]]<-net$edge[par_e,1] ##store the parent node in tracked
    par_nd_no<-length(tracked) ##record where the parent node is in tracked

    ##check if nd is a reticulate node
    hyb_e<-which(net$reticulation[,2]==nd)
    if(length(hyb_e)!=0){##if there is a parental hybrid edge it is a reticulate nd
      tracked[[length(tracked)+1]] <- net$reticulation[hyb_e,1] ##The parent of the hybrid edge
      par_hyb_nd_no<-length(tracked) ##record where the parent node is in tracked
      net$reticulation<-net$reticulation[-hyb_e,]##delete the hybrid edge
      net$inheritance<-net$inheritance[-hyb_e]
    }

    ##delete nd and par_e. Then update attributes of net
    net$edge<-net$edge[-par_e,]
    if(nd %in% (1:length(net$tip.label))){ ##Delete the node from tip.label if it is a tip
      net$tip.label<-net$tip.label[-nd]
    }
    items<-doBookKeeping(net,nd,added=F,tracked) ##items is a list containing updated net and tracked

    ###Recursive Calls###
    ##If present, check parent node of the hybrid edge to see if it needs deleting/fusing
    if(length(hyb_e)!=0){
      if( tracked[[par_nd_no]]!=tracked[[par_hyb_nd_no]] ){ ## only do stuff if the parents aren't the same node
        ##Check par_nd to see if that needs deleting/fusing
        items<-delNode(net=items$net,nd=items$tracked[[par_hyb_nd_no]], tracked=items$tracked)
      }
    }
    ##Check par_nd to see if that needs deleting/fusing
    items<-delNode(net=items$net, nd=items$tracked[[par_nd_no]], tracked=items$tracked)

  ##################
  ###One Kid Case###
  ##################
  }else if( (no_kids+no_hybkids)==1){
    ##we want to fuse the parent and child edges of nd
    ##determine if the edge to the kid is a hyb edge or normal one

    ###Kid is a Normal Edge###
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
        net$edge<-net$edge[-par_e,] ##remove  parent edge
        net$edge.length[-par_e]
        net$Nnode<-net$Nnode-1
        items<-doBookKeeping(net=net,nd=nd,added=F,tracked=tracked)
      }else{ ##nd is the root, do nothing
        items<-list(net=net,tracked=tracked)
      }

    ###Kid is a Hybrid Edge###
    }else if(no_hybkids==1){
      par_e<-which(net$reticulation[,2]==nd)
      par_nd<-net$reticulation[par_e,1] ##store the parent node
      kid_nd<-net$edge[kid_edges,2]

      ##Update Network
      net$reticulation[hybkid_edges,1]<-par_nd ##have the child edge point from par_nd to kid_nd
      net$edge<-net$edge[-par_e,] ##remove  parent edge
      net$edge.length[-par_e]
      net$Nnode<-net$Nnode-1
      items<-doBookKeeping(net=net,nd=nd,added=F,tracked=tracked)
    }else{
      stop("Something went awry")
    }
  }else{
    stop(paste("nd had ",(no_kids+no_hybkids)," kids. We should never be at a node with more than 1 kid",sep="") )
  }
  return(items)
}

doBookKeeping<-function(net,nd,added=T,tracked=list()){
  if(added){
    ##starting at nd, increment all nodes up by 1
    increment<- 1
  }else{
    ##starting at nd, increment all nodes down by 1
    increment<- -1
  }
  net$edge[net$edge>=nd]<- net$edge[net$edge>=nd]+increment
  if("evonet" %in% class(net)){
    net$reticulation[net$reticulation>=nd]<- net$reticulation[net$reticulation>=nd]+increment
  }
  if(!is.null(tracked)){
    mylilfun<-function(x,increment,nd){
      x[x>=nd]<-x[x>=nd]+increment
      return(x)
    }

    if(typeof(tracked)=="list"){
      tracked<-lapply(tracked,mylilfun,increment=increment,nd=nd)
    } else{
      tracked<-mylilfun(tracked,increment,nd)
    }
  }
  return(list(net=net,tracked=tracked ))
}


