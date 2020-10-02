
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
set.seed(10)
net_age <- sim.bdh.age(age=2,numbsim=1,
            lambda=0.4, mu=0.1,
            nu=0.1,hybprops=c(0.5,0.3,0.2),
            alpha=1,beta=1,
            frac=1,mrca = F,complete=TRUE,stochsampling = F)[[1]]


net_taxa <- sim.bdh.taxa(n=5, numbsim=1,
                         lambda=0.4, mu=0.3,
                         nu=0.05,hybprops = c(0.5,0.3,0.2),
                         alpha=1,beta=1,
                         frac=1, complete=TRUE,stochsampling = F)
```

We can also create newick strings for the networks we generated. It is worth noting that the newick here includes the inheritance probabilities, unlike the write.evonet() function found in *ape*

``` r
write.net(net_age)
#> [1] "((((t2:0.5195518503,t4:0.7195636879):0.8929270535,(t6:0.03552307698,#H16:0::0.531064524082467):1.576967664):0.3072734784,#H10:0::0.985566094052047):0.05032296758,(((t1:0.07076038526,(t3:0.03552307698)#H16:0.03523730828::0.468935475917533):1.086977901,(t7:0.378549421,t5:0.378549421):0.7791888649):0.7620259339)#H10:0.05032296758::0.0144339059479534):0.02991281264;"
```

We can try plotting this network in R

``` r
plot(net_age)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

What a beaut. However, we got lucky with this tree, some trees won't plot in R. Instead we can try plotting the network in [IcyTree](https://icytree.org/), this often works when plotting in R fails. Although, oddly enough, if we try plotting the network we just made in IcyTree it doesn't work. Plotting networks in Julia is the only plotting function I've seen that works with any network but even that plotting function has its quirks. Unfortunately, it might be worth considering a network plotting function at some point.
