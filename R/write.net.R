#' Write a Network in Parenthetic Format
#'
#' @description This function writes a network to file in the Extended Newick format.
#'
#' @details The node labels and the root edge length, if available, are written in the file.
#'
#' If inheritance probabilities are included in the network object as the 'inheritance' attribute, they are also written to file.
#'
#' If tree.names == TRUE then a variant of the Newick format is written for which the name of a tree precedes the Newick format tree (parentheses are eventually deleted beforehand). The tree names are taken from the names attribute if present (they are ignored if tree.names is a character vector).
#'
#' The tip labels (and the node labels if present) are checked before being printed: the leading and trailing spaces, and the leading left and trailing right parentheses are deleted; the other spaces are replaced by underscores; the commas, colons, semicolons, and the other parentheses are replaced with dashes
#'
#' @param net A phylogenetic network of class `evonet`. The network may include an optional attribute `inheritance`` that represents the inheritance probbilities on the edges found in the 'reticulation' attribute
#' @param file a file name specified by either a variable of mode character, or a double-quoted string; if `file = ""`` (the default) then the tree is written on the standard output connection (i.e. the console).
#' @param append a logical, if TRUE the tree is appended to the file without erasing the data possibly existing in the file, otherwise the file (if it exists) is overwritten (`FALSE`` the default).
#' @param digits a numeric giving the number of digits used for printing branch lengths.
#' @param tree.names either a logical or a vector of mode character. If TRUE then any tree names will be written prior to the tree on each line. If character, specifies the name of "phylo" objects which can be written to the file.
#'
#' @return a vector of mode character if file = "", none (invisible NULL) otherwise
#' @export
#'
#' @examples
write.net<-function(net,file="",append = FALSE, digits = 10, tree.names = FALSE){


  w.tree(phy=net,file=file,append=append,digits=digits,tree.names = tree.names)

  }




##This function returns the node timings.
##branching.times() and node.depth.edgelength() both call reorder() which causes memory errors in the RCPP calls
##This function serves the same purpose as those without making that call - makes phytools a dependency
node.times<-function(tree){
  node_times<-node.depth.edgelength(tree)
  ##Trying to remove dependency on 'phytools'
  # nd_heights<-as.numeric(nodeHeights(tree))
  # nd_nums<-as.integer(tree$edge)
  #
  # num_nds<-tree$Nnode+length(tree$tip.label)
  # node_times<-rep(NA,length(num_nds))
  # for(i in 1:length(nd_heights)){
  #   node_times[nd_nums[i]] <- nd_heights[i]
  # }
  return(node_times)
}



e2p <-  function(x)
{
  inh<-x$inheritance
  x$inheritance<-rep(NA,nrow(x$edge))
  nTips <- as.integer(length(x$tip.label))
  if (!is.null(x$edge.length)) {
    nd <- node.times(x)
    x$edge.length <- c(x$edge.length, nd[x$reticulation[, 2]] - nd[x$reticulation[, 1]])
  }
  if (!is.null(x$node.label)){
    x$tip.label <- c(x$tip.label, x$node.label[x$reticulation[, 2] - nTips])
  } else {
    newLabels <- paste0("#H", x$reticulation[, 2])
    x$tip.label <- c(x$tip.label, newLabels)
    x$node.label <- rep("", x$Nnode)
    ind <- which((x$reticulation[, 2] > nTips) & !duplicated(x$reticulation[, 2]))
    x$node.label[x$reticulation[ind, 2] - nTips] <- newLabels[ind]
  }

  if(!is.null(inh) ){
    for( i in 1:nrow(x$reticulation)){
      x$inheritance[which(x$edge[,2]==x$reticulation[i,2])]<- 1-inh[i]
    }
    x$inheritance<-c(x$inheritance,inh)
  }


  nrets <- as.integer(nrow(x$reticulation))
  x$edge[x$edge > nTips] <-  x$edge[x$edge > nTips] + nrets
  x$reticulation[, 1] <- x$reticulation[, 1] + nrets
  x$reticulation[, 2] <- nTips + (1L:nrets)
  colnames(x$edge)<-c('','')
  colnames(x$reticulation)<-c('','')

  x$edge <- rbind(x$edge, x$reticulation)
  x
}

###########################################
###functions below are modified from ape###
###########################################
#TODO properly credit ape package

## write.tree.R (2010-12-07)

##   Write Tree File in Parenthetic Format

## Copyright 2002-2010 Emmanuel Paradis, Daniel Lawson, and Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.


chckLbl <- function(phy, x, ...)
{
  ## delete all leading and trailing spaces and tabs, and
  ## the leading left and trailing right parentheses:
  ## (the syntax will work with any mix of these characters,
  ##  e.g., "    ( ( ((  " will correctly be deleted)
  x <- gsub("^[[:space:]\\(]+", "", x)
  x <- gsub("[[:space:]\\)]+$", "", x)
  ## replace all spaces and tabs by underscores:
  x <- gsub("[[:space:]]", "_", x)
  ## remove all commas, colons, and semicolons
  x <- gsub("[,:;]", "", x)
  ## replace left and right parentheses with dashes:
  x <- gsub("[\\(\\)]", "-", x)
  ## delete extra underscores and extra dashes:
  x <- gsub("_{2,}", "_", x)
  x <- gsub("-{2,}", "-", x)

  ## Output nhx annotations if they exist.
  tags <- phy$.tags
  if (!is.null(tags)) {
    for (i in 1:length(x)) {
      if (length(tags[[i]]) > 0) {
        cur.tags <- tags[[i]]
        str <- "[&&NHX:"
        tag.and.value <- paste(names(cur.tags), as.character(cur.tags), sep='=')
        tagvalues <- paste(tag.and.value, collapse=':')
        x[i] <- paste(x[i], str, tagvalues, ']', sep='')
      }
    }
  }

  x
}

w.tree <- function(phy, file = "", append = FALSE, digits = 10, tree.names = FALSE)
  {
    if (!(inherits(phy, c("phylo", "multiPhylo"))))
      stop("object \"phy\" has no trees")

    if (inherits(phy, "phylo")) phy <- c(phy)
    N <- length(phy)
    res <- character(N)

    if (is.logical(tree.names)) {
      if (tree.names) {
        tree.names <-
          if (is.null(names(phy))) character(N)
        else names(phy)
      } else tree.names <- character(N)
    }

    for (i in 1:N)
      p<-e2p(phy[[i]])

      res[i] <- .w.tree2(p, digits = digits,
                             tree.prefix = tree.names[i])

    if (file == "") return(res)
    else cat(res, file = file, append = append, sep = "\n")
  }

.w.tree2 <- function(phy, digits = 10, tree.prefix = "")
{
  brl <- !is.null(phy$edge.length)
  nodelab <- !is.null(phy$node.label)
  phy$tip.label <- chckLbl(phy, phy$tip.label)
  if (nodelab) phy$node.label <- chckLbl(phy, phy$node.label)
  f.d <- paste("%.", digits, "g", sep = "")
  cp <- function(x){
    STRING[k] <<- x
    k <<- k + 1
  }
  add.internal <- function(i) { ## Internal inheritance probs not added correctly
    cp("(")
    desc <- kids[[i]]
    for (j in desc) {
      if (j > n) add.internal(j)
      else add.terminal(ind[j])
      if (j != desc[length(desc)]) cp(",")
    }
    cp(")")
    if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
    if (brl) {
      cp(":")
      cp(sprintf(f.d, phy$edge.length[ind[i]]))
      if(!is.na(phy$inheritance[ind[i]])){
        cp("::") ##TODO make a case when we have inheritance probabilities but no edge lengths (?)
        cp( phy$inheritance[ind[i]] )
      }
    }
  }
  add.terminal <- function(i) {
    tip_lab<-phy$tip.label[phy$edge[i, 2]]
    cp(tip_lab)
    if (brl) {
      cp(":")
      cp(sprintf(f.d, phy$edge.length[i]))
      if( tip_lab %in% phy$node.label && (!is.na(phy$inheritance[i])) ){
        cp("::") ##TODO make a case when we have inheritance probabilities but no edge lengths (?)
        cp( phy$inheritance[i] )
      }

    }
  }

  n <- length(phy$tip.label)

  ## borrowed from phangorn:
  parent <- phy$edge[, 1]
  children <- phy$edge[, 2]
  kids <- vector("list", n + phy$Nnode)
  for (i in 1:length(parent))
    kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])

  ind <- match(1:max(phy$edge), phy$edge[, 2])

  LS <- 4*n + 5
  if (brl) LS <- LS + 4*n
  if (nodelab)  LS <- LS + n
  if (inherits(phy,"evonet")) LS <- LS + (2*nrow(phy$reticulation))

  STRING <- character(LS)
  k <- 1
  cp(tree.prefix)
  cp("(")
  getRoot <- function(phy)
    phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
  root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
  desc <- kids[[root]]
  for (j in desc) {
    if (j > n) add.internal(j)
    else add.terminal(ind[j])
    if (j != desc[length(desc)]) cp(",")
  }

  if (is.null(phy$root.edge)) {
    cp(")")
    if (nodelab) cp(phy$node.label[1])
    cp(";")
  }
  else {
    cp(")")
    if (nodelab) cp(phy$node.label[1])
    cp(":")
    cp(sprintf(f.d, phy$root.edge))
    cp(";")
  }
  paste(STRING, collapse = "")
}





