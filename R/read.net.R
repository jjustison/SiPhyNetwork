## read.tree.R (2010-09-27)

##   Read Tree Files in Parenthetic Format

## Copyright 2002-2010 Emmanuel Paradis, Daniel Lawson and Klaus Schliep

## This file is a modified part of the R-package 'ape'.


tree.build2 <- function(tp)
{
  add.internal <- function() {
    edge[j, 1] <<- current.node
    edge[j, 2] <<- current.node <<- node <<- node + 1L
    index[node] <<- j # set index
    j <<- j + 1L
  }
  add.terminal <- function() {
    edge[j, 1] <<- current.node
    edge[j, 2] <<- tip
    index[tip] <<- j # set index
    X <- unlist(strsplit(tpc[k], ":"))
    tip.label[tip] <<- X[1]
    edge.length[j] <<- as.numeric(X[2])
    inheritance[j] <<-as.numeric(X[4])
    k <<- k + 1L
    tip <<- tip + 1L
    j <<- j + 1L
  }
  go.down <- function() {
    l <- index[current.node]
    X <- unlist(strsplit(tpc[k], ":"))
    node.label[current.node - nb.tip] <<- X[1]
    edge.length[l] <<- as.numeric(X[2])
    inheritance[j] <<-as.numeric(X[4])
    k <<- k + 1L
    current.node <<- edge[l, 1]
  }
  if (!length(grep(",", tp))) {
    obj <- list(edge = matrix(c(2L, 1L), 1, 2))
    tp <- unlist(strsplit(tp, "[\\(\\):;]"))
    obj$edge.length <- as.numeric(tp[3])
    obj$Nnode <- 1L
    obj$tip.label <- tp[2]
    if (tp[4] != "") obj$node.label <- tp[4]
    class(obj) <- "phylo"
    return(obj)
  }

  tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
  tpc <- tpc[nzchar(tpc)]
  ## the following 2 lines are (slightly) faster than using gsub()
  tsp <- unlist(strsplit(tp, NULL))
  skeleton <- tsp[tsp %in% c("(", ")", ",", ";")]
  nsk <- length(skeleton)
  nb.node <- sum(skeleton == ")")
  nb.tip <- sum(skeleton == ",") + 1
  ## We will assume there is an edge at the root;
  ## if so, it will be removed and put into a vector
  nb.edge <- nb.node + nb.tip
  node.label <- character(nb.node)
  tip.label <- character(nb.tip)

  edge.length <- numeric(nb.edge)
  inheritance <- numeric(nb.edge) ##vectors for reading inheritance of Extended Newick format. added 2021-03-05
  edge <- matrix(0L, nb.edge, 2)
  current.node <- node <- as.integer(nb.tip + 1) # node number
  edge[nb.edge, 2] <- node
  index <- numeric(nb.edge + 1) # hash index to avoid which
  index[node] <- nb.edge

  ## j: index of the line number of edge
  ## k: index of the line number of tpc
  ## tip: tip number
  j <- k <- tip <- 1L

  for (i in 2:nsk) {
    if (skeleton[i] == "(") add.internal() # add an internal branch (on top)
    if (skeleton[i] == ",") {
      if (skeleton[i - 1] != ")") add.terminal() # add a terminal branch
    }
    if (skeleton[i] == ")") {
      if (skeleton[i - 1] != ")") { # add a terminal branch and go down one level
        add.terminal()
        go.down()
      }
      if (skeleton[i - 1] == ")") go.down() # go down one level
    }
  }
  edge <- edge[-nb.edge, ]
  obj <- list(edge = edge, Nnode = nb.node, tip.label = tip.label)
  root.edge <- edge.length[nb.edge]
  edge.length <- edge.length[-nb.edge]
  if (!all(is.na(edge.length))){ # added 2005-08-18
    obj$edge.length <- edge.length
  }
  if(any(!is.na(inheritance))){ # added 2021-03-05
    obj$prob <- inheritance
  }
  if (is.na(node.label[1])) node.label[1] <- ""
  if (any(nzchar(node.label))) obj$node.label <- node.label
  if (!is.na(root.edge)) obj$root.edge <- root.edge
  class(obj) <- "phylo"
  obj
}

read.tree2 <- function(file = "", text = NULL, tree.names = NULL, skip = 0,
                      comment.char = "", keep.multi = FALSE, ...)
{
  if (!is.null(text)) {
    if (!is.character(text))
      stop("argument 'text' must be of mode character")
    tree <- text
  } else {
    tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                 skip = skip, comment.char = comment.char, ...)
  }
  
  ## Suggestion from Eric Durand and Nicolas Bortolussi (added 2005-08-17):
  if (identical(tree, character(0))) {
    warning("empty character string.")
    return(NULL)
  }
  
  ## make a single string
  if (length(tree) > 1) tree <- paste(tree, collapse = "")
  
  single_quotes <- function(x, z) {
    x <- charToRaw(x)
    z <- which(x == as.raw(39))
    if (length(z) %% 2) stop("wrong number of single quotes around labels")
    l <- length(z) / 2
    opening <- z[c(TRUE, FALSE)]
    closing <- z[c(FALSE, TRUE)]
    from <- c(1, closing + 1L)
    to <- c(opening - 1L, length(x))
    i <- mapply(":", from = from, to = to, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    keep <- lapply(i, function(i) x[i])
    tmp_label <- paste0("IMPROBABLEPREFIX", 1:l, "IMPROBABLESUFFIX")
    tmpLabsRaw <- lapply(tmp_label, charToRaw)
    n <- 2 * l + 1L
    res <- vector("list", n)
    res[seq(1, n, 2)] <- keep
    res[seq(2, n - 1, 2)] <- tmpLabsRaw
    tree <<- rawToChar(unlist(res))
    i <- mapply(":", from = opening, to = closing, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    orig_label <- lapply(i, function(i) x[i])
    sapply(orig_label, rawToChar)
  }
  
  ## replace labels with single quotes (if needed)
  SINGLE.QUOTES.FOUND <- grepl("'", tree)
  if (SINGLE.QUOTES.FOUND) tmp_label <- single_quotes(tree)
  
  y <- unlist(gregexpr(";", tree))
  if (all(y == -1))  {
    warning("no semicolon(s) [end(s) of tree] found")
    return(NULL)
  }
  
  ## if one tree per line much faster
  if (identical(y, nchar(tree))) { # check if always one tree per line
    Ntree <- length(y)
    STRING <- character(Ntree)
    for (i in 1:Ntree) {
      STRING[i] <- gsub("\\[[^]]*\\]", "", tree[i]) # delete comments (fix 2015-01-12)
    }
  } else {
    tree <- unlist(strsplit(tree, NULL))
    y <- which(tree == ";")
    Ntree <- length(y)
    x <- c(1, y[-Ntree] + 1)
    ## Suggestion from Olivier Francois (added 2006-07-15):
    if (is.na(y[1])) return(NULL)
    STRING <- character(Ntree)
    for (i in 1:Ntree) {
      tmp <- paste0(tree[x[i]:y[i]], collapse = "")
      STRING[i] <- gsub("\\[[^]]*\\]", "", tmp) # delete comments (fix 2015-01-12)
    }
  }
  
  ## remove possible leading and trailing underscores
  STRING <- gsub("^_+|_+$", "", STRING)
  STRING <- gsub("[ \t]", "", STRING) # spaces and TABs within quoted labels are not deleted
  
  getTreeName <- function(x) {
    res <- rep("", length(x))
    i <- regexpr("\\(", x)
    s <- i > 1
    if (any(s)) res[s] <- substr(x[s], 1, i[s] - 1)
    res
  }
  
  tmpnames <- getTreeName(STRING)
  if (is.null(tree.names) && any(nzchar(tmpnames))) tree.names <- tmpnames
  
  colon <- grep(":", STRING)
  if (!length(colon)) {
    stop("we need a network with branch lengths")
    #obj <- lapply(STRING, clado.build2)
  } else if (length(colon) == Ntree) {
    obj <- lapply(STRING, tree.build2)
  } else {
    stop("we need a network with branch lengths")
    # obj <- vector("list", Ntree)
    # obj[colon] <- lapply(STRING[colon], tree.build2)
    # nocolon <- (1:Ntree)[!1:Ntree %in% colon]
    # obj[nocolon] <- lapply(STRING[nocolon], clado.build)
  }
  
  if (SINGLE.QUOTES.FOUND) {
    FOO <- function(x) {
      i <- gsub("^IMPROBABLEPREFIX|IMPROBABLESUFFIX$", "", x)
      tmp_label[as.integer(i)]
    }
    for (i in 1:Ntree) {
      lab <- obj[[i]]$tip.label
      k <- grep("IMPROBABLEPREFIX", lab)
      if (length(k)) {
        lab[k] <- FOO(lab[k])
        obj[[i]]$tip.label <- lab
      }
      lab <- obj[[i]]$node.label
      k <- grep("IMPROBABLEPREFIX", lab)
      if (length(k)) {
        lab[k] <- FOO(lab[k])
        obj[[i]]$node.label <- lab
      }
    }
  }
  if (Ntree == 1 && !keep.multi) obj <- obj[[1]] else {
    if (!is.null(tree.names)) names(obj) <- tree.names
    class(obj) <- "multiPhylo"
  }
  return(obj)
}


as.evonet.phylo2 <- function(x, ...)
{
  
  
  pos <- grep("#", x$tip.label)
  if(length(pos)==0 ){ ##if we didn't find any '#' then it is a tree and return as is
    return(x)
  }
  ind <- match(pos, x$edge[, 2])
  reticulation <- x$edge[ind, , drop = FALSE]
  inheritance<-x$prob[ind]
  edge <- x$edge[-ind, , drop = FALSE]
  nTips <- as.integer(length(x$tip.label))
  reticulation[, 2] <- as.integer(match(x$tip.label[pos], x$node.label) + nTips)
  for (i in sort(pos, TRUE)) {
    edge[edge > i ] <- edge[edge > i] - 1L
    reticulation[reticulation > i] <- reticulation[reticulation > i] - 1L
  }
  x$edge <- edge
  x$reticulation <- reticulation
  if (!is.null(x$edge.length)) x$edge.length <- x$edge.length[-ind]
  x$prob<-NULL
  x$inheritance<-inheritance
  x$tip.label <- x$tip.label[-pos]
  class(x) <- c("evonet", "phylo")
  x
}


#' Read a Network from Parenthetic Format
#'
#' @description This function reads a network from file using the Rich Newick format.
#'
#' @details
#'
#' If inheritance probabilities are included in the string, the returned `evonet` object will include an `inheritance` element. `inheritance[i]` corresponds to the inheritance probability of the hybrid edge denoted in `reticulation[i,]`
#'
#' This function also accepts the optional arguments `skip` and `tree.names`. `tree.names` is used  if there are several trees to be read and is a vector of mode character that gives names to the individual trees; if `NULL` (the default), the trees are named `"tree1"`, `"tree2"`, ...
#' The optional argument `skip` denotes the number of lines of the input file to skip before beginning to read data (this is passed directly to `scan()`).
#' @param file a file name specified by either a variable of mode character, or a double-quoted string; if `file = ""` (the default) then the tree is input on the keyboard, the entry being terminated with a blank line.
#' @param text 	alternatively, the name of a variable of mode character which contains the tree(s) in parenthetic format. By default, this is ignored (set to `NULL`, meaning that the tree is read in a file); if text is not `NULL`, then the argument file is ignored.
#' @param comment.char a single character, the remaining of the line after this character is ignored (this is passed directly to `scan()`).
#' @param ... further arguments to be passed to `scan()` and `read.tree`.
#' @return A phylogenetic network of class `evonet`.
#' @export
#'
#' @examples
#' net<-read.net(text="((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")

read.net <- function(file = "", text = NULL, comment.char = "", ...)
{
  x <- read.tree2(file = file, text = text, comment.char = comment.char, ...)
  as.evonet.phylo2(x[[1]])
  if("multiPhylo" %in% class(x)){
    lapply(x,as.evonet.phylo2)
  }else{
    as.evonet.phylo2(x)
  }
}

