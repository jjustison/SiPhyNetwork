sim.bdh.taxa.ssa.help <-function(dummy,n,numbsim,
                             lambda,mu,
                             nu, hybprops, hyb.rate.fxn=NULL,
                             hyb.inher.fxn,
                             frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE,
                             trait.model,
                             beta){

  effective_n<-getEffectiveN(n=n,frac = frac,stochsampling = stochsampling)
	phy<-sim2.bdh.origin(m=effective_n,n=0,age=Inf,
	                     lambda=lambda,mu=mu,
	                     nu=nu,
	                     hyb.inher.fxn = hyb.inher.fxn,
	                     hybprops=hybprops,hyb.rate.fxn=hyb.rate.fxn,
	                     mrca=mrca,
	                     trait.model=trait.model,
	                     beta=beta)

	##handle incomplete sampling and reconstruction
	if( "phylo" %in% class(phy)){
    phy<- handleTipsTaxa(phy=phy,complete=complete,target_ntaxa = n,current_n =effective_n)
	}

  phy
}



sim.bdh.taxa.gsa.help <-function(dummy,m,n,numbsim,
                                 lambda,mu,
                                 nu, hybprops, hyb.rate.fxn=NULL,
                                 hyb.inher.fxn,
                                 frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE,
                                 trait.model,
                                 beta){

  effective_n<-getEffectiveN(n=n,frac = frac,stochsampling = stochsampling)

  phy<-sim2.bdh.origin(m=m,n=effective_n,age=Inf,
                       lambda=lambda,mu=mu,
                       nu=nu,
                       hyb.inher.fxn = hyb.inher.fxn,
                       hybprops=hybprops,hyb.rate.fxn=hyb.rate.fxn,
                       mrca=mrca,
                       trait.model,
                       beta=beta)


  ##only do things if we have a phylogeny
  if( "phylo" %in% class(phy) ){
    phy<-internal.network.gsa(net = phy,times=phy$time_in_n,timecreation = phy$timecreation,extinct_labels = phy$extinct)
    ##Get rid of these attributes as we don't need them anymore
    phy$time_in_n<-NULL
    phy$timecreation<-NULL

    ##recompute phy$hyb_tips since we deleted tips and changed node numberings
    phy$hyb_tips <- which(1:length(phy$tip.label) %in% phy$reticulation[,1])

    phy<-handleTipsTaxa(phy=phy,complete=complete,target_ntaxa=n,current_n=effective_n)

  }

  phy
}








# sim2.net.bdh.reverse.single <- function(n,lambda,mu,nu,hybprops,frac){
#
#   maxleaf <- round(n/frac)
#   leaves<- 1:maxleaf
#   nleaves<-length(leaves)
#   nodes<- leaves
#   timecreation <- 0*leaves
#   edge <- matrix(nrow=0,ncol=2)
#   hyb_edge <- matrix(nrow=0,ncol=2)
#   num_hybs<-0
#   inheritance<-c()
#   edge.length <- vector()
#   #extinct <- vector()
#   newspecies<-0
#   time=0
#   extinct=0
#   while (nleaves>0){
#     print('new round in the while loop')
#     print(paste("theese are the leaves",paste(leaves,collapse = ' ')))
#
#     spec_rate<-nleaves*lambda
#     ext_rate<-nleaves*mu
#     hyb_rate<-choose(nleaves,2)*nu
#     total_rate<-spec_rate+ext_rate+hyb_rate
#
#     timestep <-rexp(1, total_rate)
#     time = time+timestep
#     timecreation <- c(timecreation,time)
#     randevent <- runif(1,0,1)   #event speciation or extinction
#     if (randevent <= (spec_rate/total_rate) ){ #speciation
#       print("speciation event")
#       if (length(leaves)>1) {
#         newspecies<- newspecies-1
#         nodes<-c(nodes,newspecies)
#
#         species <- sample((1:nleaves),2)
#         spec1 <- (which(nodes==leaves[species[1]]))
#         spec2 <- (which(nodes==leaves[species[2]]))
#
#         leaves<- c(leaves,newspecies)
#         edge<-rbind(c(newspecies,leaves[species[1]]),edge)
#         edge.length<-c((time-timecreation[spec1]),edge.length)
#         edge<-rbind(c(newspecies,leaves[species[2]]),edge)
#         edge.length<-c((time-timecreation[spec2]),edge.length)
#         if (species[1]>species[2]){
#           leaves <- leaves[-species[1]]
#           leaves <- leaves[-species[2]]
#         } else {
#           leaves <- leaves[-species[2]]
#           leaves <- leaves[-species[1]]
#         }
#       } else {
#         print('special case speciation')
#         newspecies<- newspecies-1
#         edge <- rbind(edge,c(newspecies,leaves[1]))
#         spec <- (which(nodes==leaves[1]))
#         edge.length<- c(edge.length,(time-timecreation[spec]))
#         leaves <- leaves[-1]
#       }
#       nleaves <- nleaves-1
#     } else if (randevent <= ((spec_rate+ext_rate)/total_rate) ){ #extinction
#       print("extinction")
#       extinct <- extinct+1
#       leaves<- c(leaves,(maxleaf+extinct))
#       nodes<-c(nodes,(maxleaf+extinct))
#       nleaves<- nleaves+1
#     } else{  ##hybridization event
#       print("hybridization")
#       ##Update things that are constant between all hybridization types
#       num_hybs<-num_hybs+1
#       inheritance[num_hybs]<- rbeta(n=1,shape1 = alpha,shape2 = beta)
#
#       randevent <- runif(1,0,1)
#       if( randevent<=(hybprops[1]/sum(hybprops)) ){ #Lineage Generating
#         if(nleaves<3){ ##we can't have lineage generating hybridization if there are fewer than three lineages
#           inheritance<-inheritance[-num_hybs]
#           num_hybs<-num_hybs-1
#           print('special case hybridization')
#         }else{
#           species <- sample((1:nleaves),3)
#           spec1 <- (which(nodes==leaves[species[1]]))
#           spec2 <- (which(nodes==leaves[species[2]]))
#           spec3 <- (which(nodes==leaves[species[3]]))
#
#           leaves<- c(leaves,newspecies-1,newspecies-2)
#
#           edge<-rbind(c(newspecies-1,leaves[species[1]]),edge)
#           edge.length<-c((time-timecreation[spec1]),edge.length)
#
#           edge<-rbind(c(newspecies-2,leaves[species[2]]),edge)
#           edge.length<-c((time-timecreation[spec2]),edge.length)
#
#           edge<-rbind(c(newspecies-3,leaves[species[3]]),edge)
#           edge.length<-c((time-timecreation[spec3]),edge.length)
#
#           edge<-rbind(c(newspecies-1,newspecies-3),edge)
#           edge.length<-c(0,edge.length)
#
#           hyb_edge<-rbind(c(newspecies-2,newspecies-3),hyb_edge)
#
#           leaves<-leaves[-species]
#           nodes<-c(nodes,newspecies-(1:3))
#           nleaves<-nleaves-1
#           newspecies<-newspecies-3
#
#           timecreation <- c(timecreation,time,time)
#         }
#       }else if( randevent<=(sum(hybprops[1:2])/sum(hybprops)) ){ #Lineage Degenerative
#
#
#       }else{ ##Lineage Neutral
#         species <- sample((1:nleaves),2)
#         spec1 <- (which(nodes==leaves[species[1]]))
#         spec2 <- (which(nodes==leaves[species[2]]))
#
#         leaves<- c(leaves,newspecies-1,newspecies-2)
#
#         edge<-rbind(c(newspecies-1,leaves[species[1]]),edge)
#         edge.length<-c((time-timecreation[spec1]),edge.length)
#
#         edge<-rbind(c(newspecies-2,leaves[species[2]]),edge)
#         edge.length<-c((time-timecreation[spec2]),edge.length)
#
#         hyb_edge<-rbind(c(newspecies-1,newspecies-2),hyb_edge)
#
#         leaves<-leaves[-species]
#         nodes<-c(nodes,newspecies-(1:2))
#         newspecies<-newspecies-2
#       }
#     }
#     print(nrow(edge))
#     print(length(edge.length))
#   }
#
#   print(nrow(edge))
#   print(length(edge.length))
#
#   minval <- -min(edge)
#   for (j in 1:length(edge.length)){
#     if (edge[j,1]<0){
#       edge[j,1]<- edge[j,1]+minval+maxleaf+extinct+1
#     }
#     if (edge[j,2]<0){
#       edge[j,2]<- edge[j,2]+minval+maxleaf+extinct+1
#     }
#   }
#
#   phy <- list(edge = edge)
#   phy$tip.label <- paste("t", sample(maxleaf+extinct), sep = "")
#
#   phy$edge.length <- edge.length
#   phy$Nnode <- maxleaf+extinct
#   class(phy) <- "phylo"
#   phy<-reorder(phy)
#   phy<-list(phy,time)
#   phy
# }


