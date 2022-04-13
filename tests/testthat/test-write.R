test_that("tests that write.net creates the correct newick", {
  library(ape)
  n<-read.evonet(text = "((A:7,((B:2,C:2):3)#H1:2):3,(#H1:1,D:6):4);")
  expect_equal(write.net(n),"((A:7,((B:2,C:2):3)#H1:2):3,(D:6,#H1:1):4);")
  n$inheritance<-0.6
  expect_equal(write.net(n),"((A:7,#H1:3::0.4):3,(((B:2,C:2):3)#H1:2::0.6,D:6):4);") ## write.net makes the minor edge the hybrid leaf


})
