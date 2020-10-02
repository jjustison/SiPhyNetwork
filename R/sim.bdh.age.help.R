sim.bdh.age.help <-function(dummy,age,lambda,mu,
                           nu, hybprops,
                           alpha,beta,
                           frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE){
	out<-sim.bdh.age.loop(age=age,numbsim=1,lambda=lambda,mu=mu,nu=nu, hybprops=hybprops,alpha=alpha,beta=beta,frac=frac,mrca=mrca,complete=complete,stochsampling=stochsampling)
	out
}


sim.bdh.age.loop <-
  function(age,numbsim,lambda,mu,nu,hybprops,alpha, beta, frac=1,mrca=FALSE,complete=TRUE,stochsampling=stochsampling) {
    phy <- sim2.bdh.age(age,numbsim,lambda,mu,nu,hybprops,mrca)[[1]]


    if( class(phy)=="phylo"){
      if(!mrca){
        phy$root.edge<-age-max(getx(phy,sersampling=1)[,1])
      }
      Nextant<-phy$Nextant
      phy$Nextant<-NULL ##Not a real part of a phylo object so we remove it
      out<-determineNetwork(phy=phy,n=Nextant,lambda=lambda,mu=mu,nu=nu,hybprops=hybprops,alpha=alpha,beta=beta,frac=frac,complete=complete,stochsampling = stochsampling)

      ###TODO After sampling the tree may have no extant tips. return 0 in this case
      ###TODO  Sometimes reconstructing/sampling will create an edge on the root. I'm not handling this and putting it in the 'root.edge' at all

      phy<-out
    }

    phy
  }



sim2.bdh.age <-
  function(age,numbsim,lambda,mu,nu, hybprops,mrca){
    phy <- list()
    for (j in 1:numbsim){
      temp <- sim2.bdh(0,age,lambda,mu,nu, hybprops,mrca)
      phy <- c(phy, list(temp))
    }
    phy
  }



sim2.bdh <-
  function(n,age, lambda,mu,nu,hybprops,mrca){
    phy2 <- sim2.bdh.origin(n=n,age=age,lambda=lambda,mu=mu,nu=nu,hybprops=hybprops,mrca=mrca)
    if (class(phy2)=="phylo") {
      #Tanja 9.7.: include root age!
      phy2<-collapse.singles(phy2)
    }
    phy2
  }



sim2.bdh.origin <-
  function(n,age, lambda,mu,nu,hybprops,mrca){

    hyb_lambda<-hybprops[1]*nu
    hyb_mu<-hybprops[2]*nu
    lambda0<-lambda
    edge <- c(-1,-2)		#matrix of edges
    leaves <- c(-2)			#list of extant leaves
    timecreation <-c(0,0)		#time when species -2, -3, ... -n was created after origin
    extinct <- vector()		#list of extinct leaves
    time <-0				#time after origin
    maxspecies <- -2		#smallest species
    edge.length <- c(0)		#edge length. if 0: leaf which didn't speciate /extinct yet
    extincttree = 0
    stop = 0

    if(mrca){##fake a speciation event to start with two lineges
      timestep<-pi ##Clearly fake value for a edgelength
      age<-timestep+age ##make the simulation longer by pi since we will ignore this edgelength

      time = time+timestep			#time after origin
      species <- sample(leaves,1)  	#the leaf undergoing the next event
      del <- which(leaves == species)
      specevent <- runif(1,0,1)		#event speciation or extinction
      edgespecevent <- which(edge == species) - length(edge.length)

      edge.length[edgespecevent] <- time-timecreation[- species]
      edge <- rbind(edge,c(species,maxspecies-1))
      edge <- rbind(edge,c(species,maxspecies-2))
      edge.length <- c(edge.length,0,0)
      leaves <- c(leaves,maxspecies-1,maxspecies-2)
      maxspecies <- maxspecies-2
      leaves <- leaves[- del]
      timecreation <- c(timecreation,time,time)
      if (length(leaves) == n){
        stop = 1
      }

    }


    while (stop == 0 ){
      if (length(leaves) == 0){
        #phy2 = 0
        if (age >0) {
          phy2=0
        }
        extincttree=1
        stop = 1
      } else {
        timestep <- rexp(1,((length(leaves)*(lambda+mu))+((choose(length(leaves),2)*(hyb_lambda+hyb_mu))) ) )    #time since last event
        if ((age > 0 && (time+timestep) < age) || age == 0) {
          time = time+timestep			#time after origin
          species <- sample(leaves,1)  	#the leaf undergoing the next event
          del <- which(leaves == species)
          specevent <- runif(1,0,1)		#event speciation or extinction
          edgespecevent <- which(edge == species) - length(edge.length)
          if (( ((length(leaves)*lambda)+(choose(length(leaves),2)*hyb_lambda)) /((length(leaves)*(lambda+mu))+((choose(length(leaves),2)*(hyb_lambda+hyb_mu)))) ) > specevent) {
            edge.length[edgespecevent] <- time-timecreation[- species]
            edge <- rbind(edge,c(species,maxspecies-1))
            edge <- rbind(edge,c(species,maxspecies-2))
            edge.length <- c(edge.length,0,0)
            leaves <- c(leaves,maxspecies-1,maxspecies-2)
            maxspecies <- maxspecies-2
            leaves <- leaves[- del]
            timecreation <- c(timecreation,time,time)
            if (length(leaves) == n){
              stop = 1
            }
          } else {
            extinct <- c(extinct,leaves[del])
            leaves <- leaves[- del]
            edge.length[edgespecevent] <- time-timecreation[- species]
          }
        } else {
          stop = 1
        }
      }
    }

    if (extincttree==0){	#length pendant edge
      if (age==0){
        timestep <- rexp(1,((length(leaves)*(lambda+mu))+((choose(length(leaves),2)*(hyb_lambda+hyb_mu))) ) )      #time since last event
        time = time+timestep			#time after origin
      } else {
        time = age
      }
    }


    if (extincttree==0 || age ==0){
      #assign pendant edge length
      for (j in (1:length(leaves))){
        k = which( edge == leaves[j]  ) - length(edge.length)
        edge.length[ k ] <- time - timecreation[- leaves[j]]
      }

      nodes <- (length(leaves)+length(extinct))*2
      leaf=1
      interior=length(leaves)+length(extinct)+1
      edgetemp<-edge

      if (nodes == 2) {
        #edge = c(2,1)
        #leaves = c(1,2)
        phy2 <-1
      } else {


        for (j in (1:nodes)){
          if (sum(match(leaves,- j,0)) + sum(match(extinct,- j,0)) == 0) {
            # replace all -j values in edge by interior
            posvalue <- interior
            interior <- interior +1
          } else {
            posvalue <- leaf
            leaf <- leaf +1
          }
          replacel <- which(edge == - j)
          if (length(replacel)>0)  {
            for (k in 1:length(replacel)) {
              if ((replacel[k]-1) < length(edge.length)) {
                edge[replacel[k],1] <- posvalue
              } else {
                edge[(replacel[k]-length(edge.length)),2] <- posvalue
              }
            }
          }
        }


        phy <- list(edge = edge)
        phy$tip.label <- paste("t", sample(length(leaves)+length(extinct)), sep = "")
        phy$edge.length <- edge.length
        phy$Nnode <- length(leaves)+length(extinct)
        class(phy) <- "phylo"
        #phy2<-collapse.singles(phy)
        phy2 <- phy

        phy2$Nextant<-length(leaves) ##record the number of extant species
        if(mrca){
          ##Remove the fake root edge with length pi
          rt_num<-phy2$edge[1,1]
          phy2$edge<-phy2$edge[-1,]
          phy2$edge.length<-phy2$edge.length[-1]
          phy2$edge[phy2$edge>rt_num]<-phy2$edge[phy2$edge>rt_num]-1
          phy2$Nnode<-phy2$Nnode-1

        }

      }
    }

    phy2
  }

