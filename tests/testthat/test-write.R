test_that("tests that write.net creates the correct newick", {
  library(ape)
  library(phytools)
  n<-read.evonet(text = "((A:7,((B:2,C:2):3)#H1:2):3,(#H1:1,D:6):4);")
  expect_equal(write.net(n),"((A:7,((B:2,C:2):3)#H1:2):3,(D:6,#H1:1):4);")
  n$inheritance<-0.6
  expect_equal(write.net(n),"((A:7,((B:2,C:2):3)#H1:2::0.4):3,(D:6,#H1:1::0.6):4);")


})
