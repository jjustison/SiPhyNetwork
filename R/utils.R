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



## Corollary 2.11 of Jetten & Iersel 2016
isTreeBased <- function(net){
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
          curr_nd<-nd_parents[nd_parents!=curr_nd]
          prev_nd<-temp_nd
          zig<-F
        }
      }else{ ##zag
        ##The node should be an omnian. i.e. all children should be reticulate
        nd_children <- edges[edges[,1]==nd,2]
        if(!all(nd_children %in% hyb_nds)){
          path_finished <- T
        }else if(curr_nd %in% hyb_nds){ ##check if the node itself is reticulate
          return(F) ##two reticulate omnians are connected by a zigzag path
        }else{ ##move on to the next node
          temp_nd<-curr_nd
          curr_nd<-nd_children[nd_children!=curr_nd]
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

}



