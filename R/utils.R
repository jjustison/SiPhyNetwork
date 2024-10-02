

#' Determine whether a network is tree-child
#' @import ape
#' @importFrom stats rbeta rbinom rexp rnbinom runif
#'
#' @description This function determines whether a network is tree-child
#'
#' @param net A phylogenetic network of class `evonet`.

#' @return A logical that is `TRUE` if the network is tree-child
#' @details A phylogenetic network is said to be tree-child if all internal nodes have at least one tree-like or leaf node as children.
#' @export
#'
#' @examples
#' net<- read.net(text= "((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")
#' isTreeChild(net) ##returns TRUE

isTreeChild <-function(net){

  hyb_nds <- net$reticulation[,2]
  edges<-rbind(net$edge,net$reticulation)
  parent_nds <- edges[,1]


  for(nd in parent_nds){##Make sure all parent nodes have at least one tree node child
    children <- edges[nd==edges[,1],2]
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
#' @details A phylogenetic network is said to be tree-based if it can be constructed with a base tree that has additional linking arcs added. See jetten 2016 Corollary 2.11 for the algorithm used to determine whether the network is tree-based
#'
#' @export
#'
#' @examples
#' net<- read.net(text= "((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")
#' isTreeBased(net) ##returns TRUE

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
#' net<- read.net(text= "((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")
#' getNetworkLevel(net) ##returns 1
getNetworkLevel <- function(net){
  nrets<-nrow(net$reticulation)
  if(nrets==0){ ##if there are no reticulations then the network must be level 0
    return(0)
  }
  if(nrets==1){ ##having 1 reticulation always means level 1
    return(1)
  }

  rt<-as.integer(length(net$tip.labels)+1)
  edges<- rbind(net$edge,net$reticulation)
  hyb_nds <- net$reticulation[,2]
  mode(edges)<-'integer'
  nNode <- length(net$tip.label)+net$Nnode ##The total number of nodes

  blobs <- biconnectedComponents(edges,rt,nNode) ##Get the biconnected components

  net_level<- 1 ##any network with reticulations is at least level 1
  for(blob in blobs){
    blob<-blob+1
    blob_nds<-unique(as.vector(edges[blob,])) ##get all nodes from the blob
    blob_nds<-blob_nds[blob_nds %in% hyb_nds]
    net_level<-max(c(net_level,length(blob_nds)))
  }

  return(net_level)
}


#' Determine whether a phylogeny is FU-stable
#'
#' @description This function assesses whether a network is FU-stable
#'
#' @param net A phylogenetic network of class `evonet`.
#' @return A logical that is `TRUE` if the network is FU-stable
#' @export
#'
#' @examples
#' net<- read.net(text= "((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")
#' isFUstable(net) ##returns TRUE
isFUstable <- function(net){

  hyb_nds<-net$reticulation[,2]
  if(length(hyb_nds)==0){##No hybridizations
    return(TRUE)
  }
  edges<-rbind(net$edge,net$reticulation)
  ## First check to see if compressed
  children<- edges[edges[,1] %in% hyb_nds,2] ##children of the hyb_nodes
  if(any(children %in% hyb_nds)){ ##none of the children should be reticulate
    return(FALSE)
  }

  ##Now check to see if any of the tree vertices have the same children
  internal_nds <- unique(as.vector(net$edge))
  internal_nds <- internal_nds[internal_nds > length(net$tip.label)] ##internal nodes will have an index greater than the tips
  tree_nds <- internal_nds[!(internal_nds %in% hyb_nds)] ##only worry about tree nodes.
  tree_children <- list() ##store the children of all the tree nodes
  for(nd in tree_nds){
    nd_children<-edges[edges[,1]==nd,2]
    for(nd2_children in tree_children){
      if(identical(sort(nd_children),sort(nd2_children))){ ##nd and nd2 share the same children
        return(FALSE)
      }
    }
    tree_children[[nd]] <- nd_children
  }
  return(TRUE) ##Network is both compressed and has no tree vertexes with the same set of children

}


#' #' Determine whether a phylogeny is Normal
#'
#' @importFrom rstackdeque rpqueue empty insert_back peek_front without_front
#'
#' @description This function assesses whether a network is Normal
#'
#' @param net A phylogenetic network of class `evonet`.

#' @return A logical that is `TRUE` if the network is Normal
#' @export
#'
#' @examples
#' net<- read.net(text= "((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")
#' isNormal(net) ##returns TRUE
isNormal <-function(net){

  if(!isTreeChild(net)){ ##The Network can only be normal if it is tree-child
    return(FALSE)
  }

  hyb_nds<-net$reticulation[,2]
  edges<-rbind(net$edge,net$reticulation) ##put all tree and reticulate edges in one place

  tips<- 1:length(net$tip.label)

  times_visited<-rep(0,length(unique(edges))) ##the number of times we have visited an internal node
  num_children<-times_visited
  for(nd in unique(edges)){
    num_children[nd]<-sum(edges[,1]==nd)
  }

  descs<-list() ##This is where we'll store all the 'left' and 'right' descendants of each node
  for( nd in unique(edges)){
    descs[[nd]]<-list()
    descs[[nd]][[1]]<-integer(0) ##The descendants on the 'left' of the node
    descs[[nd]][[2]]<-integer(0) ##The descendants on the 'right' of the node
  }

  nd_queue<-rpqueue() ##This is will dictate how we search the phylogeny to find descendents
  for( nd in tips){ ##add the tips to the queue
    nd_queue<-insert_back(nd_queue,nd)
  }
  while(!empty(nd_queue)){
    ##Pop from the queue
    nd<-peek_front(nd_queue)
    nd_queue<-without_front(nd_queue)

    parents<-edges[edges[,2]==nd,1] ##get the parents of nd
    ##Only add the parents to the queue if they've been visited by all their children
    ## This ensures that we will be able to compute the descs of the parent as it ensures we have the descs of the children
    times_visited[parents]<-times_visited[parents]+1
    for( par_nd in parents){ ##only add the parent node if it has been visited by each child node
      if(times_visited[par_nd]== num_children[par_nd]){
        nd_queue<-insert_back(nd_queue,par_nd)
      }
    }

    children<-edges[edges[,1]==nd,2] ##The children of nd
    if(length(children) != num_children[nd]){
      stop("the number of children don't match the assumed number of children")
    }
    counter<-1 ## index for whether which child we are at
    for(child_nd in children){
      descs[[nd]][[counter]]<- unique(c(unlist(descs[[child_nd]]),child_nd))
      counter<-counter+1
    }
  }

  shortcut_edges <- edges[edges[,2] %in% hyb_nds,] ##get potential 'shortcut' edges to a hybrid node
  shortcut_edges <- shortcut_edges[!(shortcut_edges[,1] %in% hyb_nds),] ## we care about the parent nodes that aren't hybrids themselves

  for(rw  in seq_len(nrow(shortcut_edges))){ ## we now check if there are multiple paths from the parent node to the hybrid child
    par_nd <- shortcut_edges[rw,1]
    hyb_nd <- shortcut_edges[rw,2]

    if( sum(unlist(descs[[par_nd]])==hyb_nd)>1 ){
      ##There are multiple paths if the hyb_nd appears more than once in the descs list
      ##We also know that one of those paths is a 'shortcut' edge since we used the edge from par_nd to hyb_nd as reference
      return(FALSE)
    }
  }
  return(TRUE)
}


vcv_net <- function(net,tol=1e-8){
  hyb_nds<-net$reticulation[,2]
  edges<-rbind(net$edge,net$reticulation)
  rt<-as.integer(length(net$tip.label)+1)

  ##TODO add hybrid edge branch lengths to bl vector


  nnodes<-max(edges)

  vcv_mat<-matrix(data = NA,nrow=nnodes,ncol = nnodes)

  nd_queue<-rpqueue() ##This is will dictate how we search the phylogeny to build the vcv
  rt_descs<-edges[edges[,1]==rt,2] ##Get the children of rt
  for(desc in rt_descs){ ##add the descs to the nd queue
    nd_queue<-insert_back(nd_queue,desc)
  }
  while(!empty(nd_queue)){
    ##Pop from the queue
    nd<-peek_front(nd_queue)
    nd_queue<-without_front(nd_queue)

    ##Check if the nd is a hybrid and deal with it accordingly
    if(!(nd %in% hyb_nds)){ ##non-hybrid case
      ##nd varies with itself and its children during the edge that leads to nd
      children<-edges[edges[,1]==nd,2]
      bran_len <- net$edge.length[edges[,2]==nd]
      vcv_mat[c(children,nd),nd] <- vcv_mat[nd,nd] + bran_len

      ##have the children 'inherit' the same covariances as the parent
      vcv_mat[,children] <- vcv_mat[,nd]

      ##add the children to the queue
      for(child in children){
        nd_queue <-insert_back(nd_queue,child)
      }

    }else{ ##hybrid case

    }


  }

}



