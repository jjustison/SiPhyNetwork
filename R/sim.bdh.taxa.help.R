sim.bdh.taxa.help <-function(dummy,n,lambda,mu,nu,hybprops,alpha,beta,frac=1,complete=TRUE,stochsampling=FALSE){
	out<-sim.net.bdh.taxa.loop(n,1,lambda,mu,nu,hybprops,alpha,beta,frac,complete,stochsampling)


	if(class(out[[1]][[1]])=="phylo"){
	  out[[1]][[1]]<-determineNetwork(phy=out[[1]][[1]],n=n,lambda=lambda,mu=mu,nu=nu,hybprops=hybprops,alpha=alpha,beta=beta,frac=frac,complete=complete,stochsampling=stochsampling)
	}
	out
}

sim.net.bdh.taxa.loop <-function(n,numbsim,lambda,mu,nu,hybprops,alpha,beta,frac=1,complete=TRUE,stochsampling=FALSE){

  phy <- sim2.net.bdh.reverse(round(n/frac),numbsim,lambda,mu,nu,hybprops)

  phy
}


sim2.net.bdh.reverse <- function(n,numbsim,lambda,mu,nu,hybprops){	phy <- list()
time<-vector()
frac<-1
for (j in 1:numbsim){
  temp <- sim2.net.bdh.reverse.single(n,lambda,mu,nu,hybprops,frac)
  phy <- c(phy, list(temp[[1]]))
  time<-c(time,temp[[2]])
}
phy2<-list(phy,time)
phy2
}

sim2.net.bdh.reverse.single <- function(n,lambda,mu,nu,hybprops,frac){
  hyb_lambda<-hybprops[1]*nu
  hyb_mu<-hybprops[2]*nu

  maxleaf <- round(n/frac)
  leaves<- 1:maxleaf
  nodes<- leaves
  timecreation <- 0*leaves
  edge <- vector()
  edge.length <- vector()
  #extinct <- vector()
  newspecies<-0
  time=0
  extinct=0
  ##TODO currently handling which type of speciation and 'extinction' occur after generating the phylogeny. It may be more copmutationally efficient to do that here
  while (length(leaves)>0){

    lambda1<-length(leaves)*lambda #speciation
    mu1<-length(leaves)*mu #extinction
    lambda2<-choose(length(leaves),2)*hyb_lambda # speciating hybridization
    mu2<-choose(length(leaves),2)*hyb_mu # lineage degenerative hybridization

    timestep <-rexp(1,(mu1+mu2+lambda1+lambda2) )
    time = time+timestep
    timecreation <- c(timecreation,time)
    specevent <- runif(1,0,1)   #event speciation or extinction
    if ( ((lambda1+lambda2)/((lambda1+lambda2)+(mu1+mu2))) >= specevent) { #speciation
      if (length(leaves)>1) {
        newspecies<- newspecies-1
        nodes<-c(nodes,newspecies)

        species <- sample((1:(length(leaves))),2)
        spec1 <- (which(nodes==leaves[species[1]]))
        spec2 <- (which(nodes==leaves[species[2]]))

        leaves<- c(leaves,newspecies)
        edge<-rbind(c(newspecies,leaves[species[1]]),edge)
        edge.length<-c((time-timecreation[spec1]),edge.length)
        edge<-rbind(c(newspecies,leaves[species[2]]),edge)
        edge.length<-c((time-timecreation[spec2]),edge.length)
        if (species[1]>species[2]){
          leaves <- leaves[-species[1]]
          leaves <- leaves[-species[2]]
        } else {
          leaves <- leaves[-species[2]]
          leaves <- leaves[-species[1]]
        }
      } else {
        newspecies<- newspecies-1
        edge <- rbind(edge,c(newspecies,leaves[1]))
        spec <- (which(nodes==leaves[1]))
        edge.length<- c(edge.length,(time-timecreation[spec]))
        leaves <- leaves[-1]
      }
    } else { #extinction
      extinct <- extinct+1
      leaves<- c(leaves,(maxleaf+extinct))
      nodes<-c(nodes,(maxleaf+extinct))
    }
  }

  minval <- -min(edge)
  for (j in 1:length(edge.length)){
    if (edge[j,1]<0){
      edge[j,1]<- edge[j,1]+minval+maxleaf+extinct+1
    }
    if (edge[j,2]<0){
      edge[j,2]<- edge[j,2]+minval+maxleaf+extinct+1
    }
  }

  phy <- list(edge = edge)
  phy$tip.label <- paste("t", sample(maxleaf+extinct), sep = "")

  phy$edge.length <- edge.length
  phy$Nnode <- maxleaf+extinct
  class(phy) <- "phylo"
  phy2<-collapse.singles(phy)
  phy2<-reorder(phy2)
  phy2<-list(phy2,time)
  phy2
}


