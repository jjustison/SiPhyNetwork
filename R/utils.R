#' Determine whether a network is tree-child
#'
#' @description This function determines whether a network is tree-child
#'
#' @param net A phylogenetic network of class `evonet`.

#' @return A logical that is `TRUE` if the network is tree-child
#' @details A phylogenetic network is said to be tree-child if all internal nodes have at least one tree-like or leaf node as children.
#' @export
#'
#' @examples

isTreeChild <-function(net){
  parent_nds <- net$reticulation[,1]
  hyb_nds <- net$reticulation[,2]

  for(nd in parent_nds){##Make sure all parent nodes have at least one tree node child
    children <- net$edge[nd==net$edge,2]
    if( sum(!(children %in% hyb_nds))==0){ ##all children are hybrid nodes
      return(F)
    }
  }
  return(T)
}



#' Determine whether a network is tree-based
#'
#' @description This function determines whether a network is tree-based
#'
#' @param net A phylogenetic network of class `evonet`.

#' @return A logical that is `TRUE` if the network is tree-based
#' @details A phylogenetic network is said to be tree-based if it can be constructed with a base tree that has additional linking arcs added.See jetten 2016 Corollary 2.11 for the algorithm used to determine whether the network is tree-based
#'
#' @export
#'
#' @examples

isTreeBased <- function(net){
  ## Corollary 2.11 of Jetten & Iersel 2016
  hyb_nds <- net$reticulation[,2]

  edges<- rbind(net$edge,net$reticulation)

  ret_omnians <-c() ##reticulate nodes where the child is also reticulate
  for(nd in hyb_nds){
    nd_child <- edges[edges[,1]==nd,2]
    if( nd_child %in% hyb_nds ){
      ret_omnians<-c(ret_omnians,nd)
    }
  }

  for(nd in ret_omnians){

    zig<-T
    curr_nd<- edges[edges[,1]==nd,2] ##the child of nd
    prev_nd<- nd
    path_finished <-F
    while(!path_finished){ ##Look at the zigzag path of each omnian
      if(zig){
        ##The node should be a reticulate node
        if( !(curr_nd %in% hyb_nds) ){
          path_finished <- T
        }else{ ## move on to the next node

          nd_parents<- edges[edges[,2]==curr_nd,1]
          temp_nd<-curr_nd
          curr_nd<-nd_parents[nd_parents!=prev_nd]

          prev_nd<-temp_nd
          zig<-F
        }
      }else{ ##zag
        ##The node should be an omnian. i.e. all children should be reticulate
        nd_children <- edges[edges[,1]==curr_nd,2]
        if(!all(nd_children %in% hyb_nds)){
          path_finished <- T
        }else if(curr_nd %in% hyb_nds){ ##check if the node itself is reticulate
          return(FALSE) ##two reticulate omnians are connected by a zigzag path
        }else{ ##move on to the next node
          temp_nd<-curr_nd
          curr_nd<-nd_children[nd_children!=prev_nd]
          prev_nd<-temp_nd
          zig<-T
        }
      }
    }
  }
  return(TRUE)
}


#' Get the level of a Network
#'
#' @description This function gets the level of the network
#'
#' @param net A phylogenetic network of class `evonet`.

#' @return A numeric with the level of the network
#' @export
#'
#' @examples
getNetworkLevel <- function(net){
  rt<-as.integer(length(net$tip.labels)+1)
  edges<- rbind(net$edge,net$reticulation)
  mode(edges)<-'integer'
  nNode <- length(net$tip.label)+net$Nnode ##The total number of nodes

  blobs <- biconnectedComponents(edges,rt,nNode)
  return(blobs)
}

isStable <- function(net){

  hyb_nds<-net$reticulation[,2]
  edges<-rbind(net$edge,net$reticulation)
  ## First check to see if compressed
  children<- edges[edges[,1] %in% hyb_nds,2] ##children of the hyb_nodes
  if(any(children %in% hyb_nds)){ ##none of the children should be reticulate
    return(FALSE)
  }

  ##Now check to see if any of the tree vertices have the same children
  internal_nds <- unique(as.vector(net$edge))
  internal_nds <- internal_nds[internal_nds > length(net$tip.label)]
  tree_nds <- internal_nds[!(internal_nds %in% hyb_nds)]
  tree_children <- list() ##store the children of all the tree nodes
  index_added<-0 ##Keep track of the index of the tree nodes we've added
  for(nd in tree_nds){
    tree_children[[nd]] <- edges[edges[,1]==nd,2]

    for(nd2 in tree_children[seq_along(rep(NA,index_added))]){
      ##TODO compare the children of nd and nd2
    }
  }


}

