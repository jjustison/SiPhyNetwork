test_that("tests that write.net creates the correct newick", {
  library(ape)
  n<-read.evonet(text = "((A:7,((B:2,C:2):3)#H1:2):3,(#H1:1,D:6):4);")
  n$inheritance<-0.4
  expect_equal(write.net(n),"((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")
  n$inheritance<-0.6
  expect_equal(write.net(n),"((A:7,#H1:2::0.4):3,(((B:2,C:2):3)#H1:1::0.6,D:6):4);") ## write.net makes the minor edge the hybrid leaf

  ##Check Classes of networks

  expect_true(isTreeChild(n))
  expect_true(isTreeBased(n))
  expect_true(isFUstable(n))
  expect_equal(getNetworkLevel(n),1)

  n2<-read.net(text="(((((t2:1,t3:1):1)#H10:1::0.5)#H9:1.5::0.5,#H7:0.5::0.5):0.5,((((t5:1,t1:1):2,#H10:1::0.5):0.5,#H9:0.5::0.5):0.5)#H7:1::0.5);")
  expect_true(isTreeBased(n2))
  expect_false((isFUstable(n2)))
  expect_false(isTreeChild(n2))
  expect_equal(getNetworkLevel(n2),3)


  ##Check that swap minor edges in write.net works properly

  set.seed(1347) # current time
  numbsim = 100
  n = 7 # number of taxa to simulate to

  ###

  inheritance.fxn <- make.beta.draw(1,1) # beta(1,1) distribution (i.e. unif(0,1))
  hybrid_proportions <-c(0.5,  ##Lineage Generative
                         0.25, ##Lineage Degenerative
                         0.25) ##Lineage Neutral
  ssa_nets <- sim.bdh.taxa.ssa(
    #age = age,
    n = 7,
    numbsim = numbsim,
    lambda=1, # speciation rate
    mu=0.2, # extinction rate
    nu=0.25, # hybridization rate,
    hybprops = hybrid_proportions,
    hyb.inher.fxn = inheritance.fxn,
    complete=FALSE
  )
  expect_true(is.ultrametric(read.net(text=write.net(ssa_nets[[100]],swap.minor=T))))
  expect_true(is.ultrametric(read.net(text=write.net(ssa_nets[[100]],swap.minor=F))))


  ##Test whether network classes work properly

  set.seed(1234)
  nreps<-1000
  count<-1
  phys<-list()

  while(count<=nreps){
    phy<-sim.bdh.age(age=1,numbsim=1,
                     lambda=3,mu=2,
                     nu=0.25,hybprops=c(0,0.5,0.5),
                     hyb.inher.fxn = make.beta.draw(10,10))[[1]]
    if('phylo' %in% class(phy)){
      phys[[count]]<-phy
      count<-count+1
    }

  }

  tc_nets<-unlist(lapply(phys, isTreeChild)) ##Phylogenies that are Tree Child
  tb_nets<-unlist(lapply(phys, isTreeBased)) ##Phylogenies that are Tree Based
  fu_nets<-unlist(lapply(phys, isFUstable))  ##Phylogenies that are FU Stable

  expect_true(all( !tc_nets | tb_nets)) ##Tree Child networks should also be Tree Based
  expect_true(all( !tc_nets | fu_nets)) ##Tree Child networks should also be FU Stable

  expect_true(all(tb_nets | !tc_nets)) ##Not Tree Based implies not Tree Child i.e. !(!tb & tc) with De Morgans law
  expect_true(all(tb_nets | !tc_nets)) ##Not Tree Based implies not FU Stable

  expect_true(all( !fu_nets | tb_nets)) ##FU Stable implies Tree Based

  expect_true(all(fu_nets | !tc_nets)) ## Not FU Stable implies not Tree Child







})
