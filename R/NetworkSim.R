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
reconstructNetwork<-function(net){

  if(is.null(net$extinct)){
    extinct_tips<-which( net$tip.label %in% getExtinct(net))
  } else{
    extinct_tips<-getReconstructedTips(net$tip.label,net$extinct)
  }
  if(length(extinct_tips)==0){
    warning('No extinct tips found. Returning NA')
    return(NA)
  }

  net<-deleteTips(net,extinct_tips)
  net$extinct<-NULL
  return(net)
}

getReconstructedTips<- function(tip.labels,extinct){
  return( which(tip.labels %in% extinct) )
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
  deltips<-getSamplingTips(1:length(net$tip.label),rho,stochastic)
  if(length(tracked$deltips)>0){ ##only delete tips if there are tips to delete
    net<-deleteTips(net,deltips)
  }

  return(net)
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
  return(net)
}





