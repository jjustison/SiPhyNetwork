
<!-- README.md is generated from README.Rmd. Please edit that file -->
NetSim
======

<!-- badges: start -->
<!-- badges: end -->
The goal of NetSim is to generate Phylogenetic networks in a process-based manner. NetSim considers a birth-death-hybridization proccess for simulating Networks. We consider allow for hybridization to be lineage generative, lineage degenerative, and lineage neutral.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jjustison/NetSim")
```

Example
-------

This is a basic example of simulating networks up to a certain age or number of taxa:

``` r
library(NetSim)
library(phytools)
#> Warning: package 'phytools' was built under R version 3.5.3
#> Loading required package: ape
#> Loading required package: maps
#> Warning: package 'maps' was built under R version 3.5.3
library(ape)
net_age <- sim.bdh.age(age=2,numbsim=1,
            lambda=0.4, mu=0.1,
            nu=0.05,hybprops=c(0.5,0.3,0.2),
            alpha=1,beta=1,
            frac=1,mrca = F,complete=TRUE,stochsampling = F)[[1]]

##Currently not working
# net_taxa <- sim.bdh.taxa(n=5, numbsim=1,
#                          lambda=0.4, mu=0.3,
#                          nu=0.05,hybprops = c(0.5,0.3,0.2),
#                          alpha=1,beta=1,
#                          frac=1, complete=TRUE,stochsampling = F)
```

We can also create newick strings for the networks we generated

``` r
write.net(net_age)
#> [1] "adding inheritance for internal"
#> [1] "(t2:0.8791724257,(((t5:0.2574956422,#H15:0::0.653201974229887):0.3143894536,t3:0.5718850958):0.2277070172,((t6:0.01753684847,t4:0.01753684847):0.7155743147,((t1:0.2574956422)#H15:0::0.346798025770113,t7:0.1397551533):0.475615521):0.06648094984):0.07958031268):1.120827574;"
```
