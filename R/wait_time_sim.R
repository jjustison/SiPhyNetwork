#' Waiting times of the birth-death-hybridization process
#'
#' @description This function simulates a birth-death-hybridization-process and records the waiting times and population size
#'
#' @return A data frame that has three columns. The first column is the waiting time until the next event. The second column indicates the population size immediately after the event occurs. The last column is `1` when the event is a hybridizaiton event and `0` otherewise. NOTE: The first row is the root with a waiting time of 0 and the population size being the number of individuals that started the process. 
#' @export
#'


sim.times.age <- function(
    age,numbsim,
    lambda,mu,
    nu,hybprops,
    startinglineages){
  
  lapply(1:numbsim,FUN=sim.bdh.times.help,
              age=age,n=NULL,m=NULL,
              lambda=lambda,mu=mu,
              nu=nu,hybprops=hybprops,
              mrca=startinglineages)
}


sim.times.taxa.ssa <- function(
    n,numbism,
    lambda,mu,
    nu,hybprops,
    startinglineages){
  
  out<-lapply(1:numbsim,sim.bdh.times.help,
              age=NULL,n=n,m=NULL,
              lambda=lambda,mu=mu,
              nu=nu, hybprops=hybprops,
              mrca=startinglineages)
  class(out)<-c('list')
}

sim.times.taxa.gsa <- function(
    n,m,numbism,
    lambda,mu,
    nu,hybprops,
    startinglineages){
  
  out<-lapply(1:numbsim,sim.bdh.times.help,
              age=NULL,n=n,m=m,
              lambda=lambda,mu=mu,
              nu=nu, hybprops=hybprops,
              mrca=startinglineages)
  class(out)<-c('list')
}


sim.bdh.times.help <- function(dummy,
                               age,n,m,
                               lambda,mu,
                               nu,hybprops,
                               mrca){
  stop_age<-Inf
  if(!is.null(m)){ ##Stop at m for gsa
    stopping_number<-m
  }else if(!is.null(n)){ # Stop at n for ssa
    stopping_number<-n
  }else{ ##don't stop at a number for age but stop_age
    stopping_number<-Inf
    stop_age<-age
  }
  
  curr_age<-0
  nleaf<-mrca
  nevents<-2
  ltt<-list()
  ltt[[1]]<-c(curr_age,nleaf,F)
  
  
  while( !(nleaf %in% c(stopping_number,0)) ){
    
    ## make a draw for the waiting time
    event_rate<-nleaf*(lambda+mu)
    if(1<nleaf){
      event_rate<-event_rate + (choose(nleaf,2)*nu)
    }
    waiting_time<-rexp(1,rate=event_rate)
    
    if(waiting_time+curr_age>age){##we passed the end time
      ##Set waiting time st that the curr age is the max age
      waiting_time<- age - curr_age
      e_type<- F
      ltt[[nevents]]<-data.frame(waiting_time,nleaf,e_type)
      break
    }else{ ##Determine which type of event occurs 
      event_weights<-c(nleaf*lambda,             ##Speciation: sp
                 nleaf*mu,                       ##Extinction: ext 
                 choose(nleaf,2)*nu*hybprops[1], ##lin generative: hyb_pos
                 choose(nleaf,2)*nu*hybprops[2], ##lin degenerative: hyb_neg
                 choose(nleaf,2)*nu*hybprops[3]) ##lin neutral: hyb_neu
      event_weights<-event_weights/sum(event_weights) #normalize the weights
      event_type<-sample(c('sp','ext','hyb_pos','hyb_neg','hyb_neu'),1,prob = event_weights)
      
      
      if(event_type=='sp'){ #speciaition
        nleaf<-nleaf+1
        e_type<-F
      }else if(event_type=='ext'){
        nleaf<-nleaf-1
        e_type<-F
      }else if(event_type=="hyb_pos"){
        nleaf<-nleaf+1
        e_type<-T
      }else if(event_type=='hyb_neu'){
        e_type<-T
      }else if(event_type=='hyb_neg'){
        nleaf<-nleaf-1
        e_type<-T
      }
      
      curr_age<-curr_age+waiting_time
      ltt[[nevents]]<-data.frame(waiting_time,nleaf,e_type)
      nevents<-nevents+1
    }

  }
  ltt<-as.data.frame(do.call(rbind,ltt))
  colnames(ltt)<-c('waiting_time',"nspecies","hyb_event")
  return(ltt)
}


#' Waiting times to a phylogeny
#'
#' @description This function takes the waiting times and population sizes and converts it into a phylogeny
#'
#' @return a phylogenetic network of class `evonet`
#' @export
#'
waiting.to.phylo.beta<-function(ltt,beta,complete=T){
  root<-unlist(ltt[1,])
  names(root)<-NULL
  if(root[1]!=0){
    stop("Please provide a root state with 0.0 as the waiting time")
  }
  
  
  ##Set up the phylo object
  edge<-matrix(nrow=0,ncol = 2) ##TODO preallocate this by calculating how many edges we would have based on the events 
  hyb_edge<-list() ##TODO preallocate this by calculating how many edges we would have based on the events 
  ##Rt <-nleaves
  ##spec <-2
  ##ext <-0
  ##hyb_- <-1
  ##hyb_+ <- 5
  ##hyb_0 <- 3
  leaves<-c()
  leaf_weights<-c()
  extinct <- vector()
  maxspecies<- -1
  num_hybs<-0
  
  
  
  
  for(i in seq(root[2])){####Set up each edge at the root
    maxspecies<-maxspecies-1
    edge<-rbind(edge,c(-1,maxspecies))
    leaves<- c(leaves,maxspecies)
    
  }
  edge.length<-  rep(0,root[2])
  timecreation<- rep(0,root[2]+1)
  leaf_weights<- rep(1/unlist(root[2]),unlist(root[2]))
  N_curr<-root[2]
  
  time<-0
  for( i in 2:nrow(ltt)){
    time<-time+ltt[i,1]

    species <- sample(leaves,1+(ltt[i-1,2]>1),prob = leaf_weights) ## Sample species for events 
    del <- match(species,leaves) ##Get their indices in leaves
    edgespecevent <- match(species,edge) - length(edge.length) ##get their indices in edge and edge.length


    ##determine the type of event
    N_prev<-N_curr
    N_curr<- ltt[i,2]
    E_curr<- ltt[i,3]
    
    N_diff<-N_curr-N_prev
    if(N_diff==1){ ## we gained a species
      if(E_curr==0){##Speciation
        ##only use the first species
        species<-species[1]
        del<-del[1]
        edgespecevent<-edgespecevent[1]
        
        edge.length[edgespecevent] <- time-timecreation[- species]
        edge <- rbind(edge,c(species,maxspecies-1))
        edge <- rbind(edge,c(species,maxspecies-2))
        edge.length <- c(edge.length,0,0)
        leaves <- c(leaves,maxspecies-1,maxspecies-2)
        leaves <- leaves[- del]
        timecreation <- c(timecreation,time,time)
        timecreation[-species]<-time
        maxspecies <- maxspecies-2
        
        lf <- leaf_weights[del]*rbeta(1,beta,beta)
        leaf_weights <- c(leaf_weights,lf,leaf_weights[del]-lf)
        leaf_weights<-leaf_weights[-del]
      }else{##Hybridization
        num_hybs<-num_hybs+1
        ##update leaves
        leaves <- c(leaves,(maxspecies- 2:4))
        leaves <- leaves[- del]
        
        lf1 <- leaf_weights[del[1]]*rbeta(1,beta,beta) ## how weight is split for parent 1
        lf2 <- leaf_weights[del[2]]*rbeta(1,beta,beta) ## how weight is split for parent 2
        leaf_weights <- c(leaf_weights,lf1,lf2, ((leaf_weights[del[1]]-lf1)+(leaf_weights[del[2]]-lf2)) ) ##Hybrid gets the other weights from both parents: 1-lf1+1-lf2
        leaf_weights<-leaf_weights[-del]
        edge <- rbind(edge,c(species[1],maxspecies-2))  ##Create new leaf for a parent lineage
        edge <- rbind(edge,c(species[2],maxspecies-3)) ##Create a new leaf for the other parent lineage
        edge <- rbind(edge,c(species[1],maxspecies-1)) ##link the parent lineage to the hybrid node
        edge <- rbind(edge,c(maxspecies-1,maxspecies-4)) ##create a leaf for the hybrid lineage
        hyb_edge[[num_hybs]]<-c(species[2],maxspecies-1) ##link the other parent lineage to the hybrid node
        ##Update edge lengths
        edge.length[edgespecevent] <- time-timecreation[- species] ##update the edge length of the lineages with the event
        edge.length <- c(edge.length,0,0,0,0)
        
        timecreation <- c(timecreation,time,time,time,time)
        maxspecies <- maxspecies-4
      }
    }else if(N_diff==-1){ ## we lost a species
      if(E_curr==0){##Extinction
        species<-species[1]
        del<-del[1]
        edgespecevent<-edgespecevent[1]
        
        ##update leaves and edges
        extinct <- c(extinct,leaves[del])
        leaves <- leaves[- del]
        ##update edge lengths
        edge.length[edgespecevent] <- time-timecreation[- species]
        
        timecreation[-species]<-time
        leaf_weights<-leaf_weights[-del]
      }else{##Hybridization
        num_hybs<-num_hybs+1
        extinct <- c(extinct,species[2])
        leaves<-c(leaves,maxspecies-1)
        leaves <- leaves[- del]
        
        hyb_weight <- sum(leaf_weights[del]) # hybrid gets the weight of both parents 
        leaf_weights <- c(leaf_weights,hyb_weight)
        leaf_weights<-leaf_weights[-del]
        
        ##update edge matrices
        edge <- rbind(edge,c(species[1],maxspecies-1))
        hyb_edge[[num_hybs]]<-c(species[2],species[1])
        ##update edge lengths
        edge.length[edgespecevent] <- time-timecreation[- species]
        edge.length<- c(edge.length,0) ##new edge for hybrid lineage
        
        timecreation <- c(timecreation,time)
        maxspecies <- maxspecies-1
      }
    }else if(N_diff==0){##species number stayed the same
      if(E_curr==0){##Nothing/End
        
        
      }else{##Hybridization
        num_hybs<-num_hybs+1
        ##update leaves
        leaves<-c(leaves,maxspecies-1:2) ##maxspecies-1 will be the hybrid and 'recipient' lineage while maxspecies-2 will be the 'donor' lineage
        leaves <- leaves[- del]
        
        
        lf <- leaf_weights[del[2]] * rbeta(1,beta,beta) ##Leaf of the donor loses some of its weight...
        hyb_weight <- leaf_weights[del[1]] + (leaf_weights[del[2]]-lf) ##... and gives it to the hybrid
        leaf_weights <- c(leaf_weights,hyb_weight,lf)
        leaf_weights<-leaf_weights[-del]
        
        ##update edge matrices
        edge <- rbind(edge,c(species[1],maxspecies-1))
        edge <- rbind(edge,c(species[2],maxspecies-2))
        hyb_edge[[num_hybs]]<-c(species[2],species[1])
        ##update edge lengths
        edge.length[edgespecevent] <- time-timecreation[- species]
        edge.length<- c(edge.length,0,0)
        
        timecreation <- c(timecreation,time,time)
        maxspecies <- maxspecies-2
      }
      timecreation[-species]<-time
    }
  }
        
        
  #assign pendant edge length
  for (j in (1:length(leaves))){
    k = which( edge == leaves[j]  ) - length(edge.length)
    edge.length[ k ] <- time - timecreation[- leaves[j]]
  }
  timecreation[-leaves]<-time

  
  ##convert hyb_edge to a data frame
  if(num_hybs==0){
    hyb_edge<-matrix(nrow=0,ncol=2)
  }else{
    hyb_edge<-matrix(unlist(hyb_edge),nrow = num_hybs,ncol = 2,byrow = T)
  }
  colnames(hyb_edge)<-c('from','to')
  num_nodes <- abs(maxspecies)
  ## we want to create things into a 'phylo' to report what happened
  ##Convert 'maxspecies' node numberings to ape friendly numberings
  leaf=1 ##Start numbering for leaves
  interior=nleaves+length(extinct)+1 ##Start numberings. start at n+1 because 1:n is reserved for tips
  extinct_tips<-rep(NA,length(extinct)) ##record all extinct tips
  timecreation_order<-rep(NA,length(timecreation)) ##reorder timecreation so node numbers and creation match
  Nextinct<-1
  
  
  for (j in (1:num_nodes)){
    inleaves <- (leaves  %in% (-j))
    inextinct<- (extinct %in% (-j))
    if (!(any(inleaves) | any(inextinct))){ ##we are at an internal node
      replaced_value <- interior
      interior <- interior +1
    } else { ## We are at a tip of some sort (extant,extinct, or hybrid)
      

      
      replaced_value <- leaf
      leaf <- leaf +1
      if(sum(match(extinct,-j,0))==1){
        extinct_tips[Nextinct]<-replaced_value
        Nextinct<-Nextinct+1
      }
    }
    timecreation_order[replaced_value] <- timecreation[j]
    hyb_edge[hyb_edge== - j]<-replaced_value
    edge[edge == -j]<-replaced_value
  }
  
  ##create a list with all the 'phylo' bits
  phy <- list(edge = edge)
  phy$tip.label <- paste("t", sample(nleaves+length(extinct)), sep = "")
  phy$edge.length <- edge.length
  phy$Nnode <- num_nodes-(nleaves+length(extinct)) ##Nnode reports the number of internal nodes
  class(phy) <- c("evonet","phylo")
  phy$reticulation<-hyb_edge
  #phy$dists<-genetic_dists
  #phy$inheritance<-inheritance
  
  hyb_tip_indices<-match(extinct_tips,phy$reticulation[,1],0)>0
  phy$hyb_tips<-extinct_tips[hyb_tip_indices] ## Keep track of tips that have an edge going into a hybrid one. These will be fused later
  phy$extinct<-phy$tip.label[extinct_tips[!hyb_tip_indices]] ##Keep track of extinct tips. Will be used if using the reconstructed tree
  
  
  SiPhyNetwork:::handleTipsTaxa(phy=phy,complete=complete,target_ntaxa = nleaves, current_n = nleaves)
}




