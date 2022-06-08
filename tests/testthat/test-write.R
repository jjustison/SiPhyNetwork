test_that("tests that write.net creates the correct newick", {
  library(ape)
  n<-read.evonet(text = "((A:7,((B:2,C:2):3)#H1:2):3,(#H1:1,D:6):4);")
  n$inheritance<-0.4
  expect_equal(write.net(n),"((A:7,((B:2,C:2):3)#H1:2::0.6):3,(D:6,#H1:1::0.4):4);")
  n$inheritance<-0.6
  expect_equal(write.net(n),"((A:7,#H1:2::0.4):3,(((B:2,C:2):3)#H1:1::0.6,D:6):4);") ## write.net makes the minor edge the hybrid leaf


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

})
