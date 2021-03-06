Moving window analysis of cell number density of lymph node images
================

## Introduction

This repository contains source codes of moving-window cell density
analysis presented in the paper “A local regulatory T cell feedback
circuit maintains immunological homeostasis by pruning self-activated T
cells” and examples of the analysis.

## Prerequisit

The actual analyses are done in R. However, core functions of obtaining
3D geometric intersections and volumes were imported from the geometric
libraries, libigl and CGAL written in c++ through the Rcpp package.
Thus, users are required to install libigl (<https://libigl.github.io>)
and CGAL (<https://doc.cgal.org/latest/Manual/usage.html>).
Additionally, Eigen
(<https://eigen.tuxfamily.org/dox/GettingStarted.html>), Boost
(<https://www.boost.org/>), GMP (<https://gmplib.org/>) and MPFR
(<https://www.mpfr.org/>) are required to run libigl and CGAL libraries.

## How to use

Download and copy the files functions.R, mesh\_boolean\_R.cpp, and
mesh\_volume\_R.cpp to a folder where you perform your analysis.

### Load R packages and source functions.R.

In functions.R, there are various functions defined for use in the
following. Sourcing functions.R includes compilation of CPP codes
through the Rcpp package. Users are required to specify relevant
directiries. Please see comments in functions.R.

``` r
library(rgl)
library(Rcpp)
library(nat)
library(alphashape3d)
library(Rvcg)
source("functions.R")
```

### Load data.

The data should be in the form of a matrix, with rows for cells and
columns for coordinates (X, Y, Z) (first three columns) and other
informations (marker intensities or cell id.). The example is for
comparing regulatory T cell (Treg) densities surounding wildtype (WT) or
TXA23 CD4 T cells at 24 hours after transfer in gastric lymph nodes for
the isotype control and IL2 blocking condition. Data is loaded in the
list form with each element from different samples. Tregs, WT T cells,
and TXA23 T cells are loaded in seperate list objects.

#### Regulatory T cells

``` r
tregs <- list()
tregs[["ab_1"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 1/Tregs XYZ.csv")
tregs[["ab_2"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 2/Tregs XYZ.csv")
tregs[["ab_3"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 3/Tregs XYZ.csv")
tregs[["ab_4"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=2 - 2/Tregs XYZ.csv")
tregs[["ab_5"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=2 - 3/Tregs XYZ.csv")
tregs[["ab_6"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=3 - 1/Tregs XYZ.csv")
tregs[["ab_7"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=4 - 1/Tregs XYZ.csv")

tregs[["iso_1"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 1/Tregs XYZ.csv")
tregs[["iso_2"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 2/Tregs XYZ.csv")
tregs[["iso_3"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 3/Tregs XYZ.csv")
tregs[["iso_4"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 1/Tregs XYZ.csv")
tregs[["iso_5"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 2/Tregs XYZ.csv")
tregs[["iso_6"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 3/Tregs XYZ.csv")
tregs[["iso_7"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 4/Tregs XYZ.csv")
tregs[["iso_8"]] <- read.csv("data/IL2 block/isotype/isotype - n=3 - 1/Tregs XYZ.csv")
```

#### Antigen specific T cells (TXA23)

``` r
TXA23 <- list()
TXA23[["ab_1"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 1/TXA23 XYZ.csv")
TXA23[["ab_2"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 2/TXA23 XYZ.csv")
TXA23[["ab_3"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 3/TXA23 XYZ.csv")
TXA23[["ab_4"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=2 - 2/TXA23 XYZ.csv")
TXA23[["ab_5"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=2 - 3/TXA23 XYZ.csv")
TXA23[["ab_6"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=3 - 1/TXA23 XYZ.csv")
TXA23[["ab_7"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=4 - 1/TXA23 XYZ.csv")

TXA23[["iso_1"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 1/TXA23 XYZ.csv")
TXA23[["iso_2"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 2/TXA23 XYZ.csv")
TXA23[["iso_3"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 3/TXA23 XYZ.csv")
TXA23[["iso_4"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 1/TXA23 XYZ.csv")
TXA23[["iso_5"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 2/TXA23 XYZ.csv")
TXA23[["iso_6"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 3/TXA23 XYZ.csv")
TXA23[["iso_7"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 4/TXA23 XYZ.csv")
TXA23[["iso_8"]] <- read.csv("data/IL2 block/isotype/isotype - n=3 - 1/TXA23 XYZ.csv")
```

#### WT (Control) T cells

``` r
WT <- list()
WT[["ab_1"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 1/WT XYZ.csv")
WT[["ab_2"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 2/WT XYZ.csv")
WT[["ab_3"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=1 - 3/WT XYZ.csv")
WT[["ab_4"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=2 - 2/WT XYZ.csv")
WT[["ab_5"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=2 - 3/WT XYZ.csv")
WT[["ab_6"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=3 - 1/WT XYZ.csv")
WT[["ab_7"]] <- read.csv("data/IL2 block/anti-IL2/anti IL-2 - n=4 - 1/WT XYZ.csv")

WT[["iso_1"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 1/WT XYZ.csv")
WT[["iso_2"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 2/WT XYZ.csv")
WT[["iso_3"]] <- read.csv("data/IL2 block/isotype/isotype - n=1 - 3/WT XYZ.csv")
WT[["iso_4"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 1/WT XYZ.csv")
WT[["iso_5"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 2/WT XYZ.csv")
WT[["iso_6"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 3/WT XYZ.csv")
WT[["iso_7"]] <- read.csv("data/IL2 block/isotype/isotype - n=2 - 4/WT XYZ.csv")
WT[["iso_8"]] <- read.csv("data/IL2 block/isotype/isotype - n=3 - 1/WT XYZ.csv")
```

### Treg extraction centered at CD4 T cells at distances from (0,1) to (99,100) um.

This step extracts Tregs’ coordinates and marker intensities at varying
distance from (0,1) to (99,100) um centered at each CD4 T cell
(reference cells). Users can specify the distance ranges upon their
interests. We call each 1um interval as a
shell.

``` r
tregs.TXA23.extr <- extract.cells.increasing.shells( cells.list = tregs, 
                                                       refs.list = TXA23, 
                                                       shells = c(0:99) )

tregs.WT.extr <- extract.cells.increasing.shells( cells.list = tregs, 
                                                       refs.list = WT, 
                                                       shells = c(0:99) )
```

### Count the number of Tregs in each shells.

This step gives cumulative counts or shell counts of Tregs.

``` r
#TXA23
tregs.TXA23.num <- count.number.per.shell(extracted.list = tregs.TXA23.extr, shells = seq(0,99,1))
```

    ## [1] "ab_1"
    ## [1] "ab_2"
    ## [1] "ab_3"
    ## [1] "ab_4"
    ## [1] "ab_5"
    ## [1] "ab_6"
    ## [1] "ab_7"
    ## [1] "iso_1"
    ## [1] "iso_2"
    ## [1] "iso_3"
    ## [1] "iso_4"
    ## [1] "iso_5"
    ## [1] "iso_6"
    ## [1] "iso_7"
    ## [1] "iso_8"

``` r
tregs.TXA23.num.cumul <- tregs.TXA23.num$num.cumul
tregs.TXA23.num.1um <- tregs.TXA23.num$num.per.shell

#WT
tregs.WT.num <- count.number.per.shell(extracted.list = tregs.WT.extr, shells = seq(0,99,1))
```

    ## [1] "ab_1"
    ## [1] "ab_2"
    ## [1] "ab_3"
    ## [1] "ab_4"
    ## [1] "ab_5"
    ## [1] "ab_6"
    ## [1] "ab_7"
    ## [1] "iso_1"
    ## [1] "iso_2"
    ## [1] "iso_3"
    ## [1] "iso_4"
    ## [1] "iso_5"
    ## [1] "iso_6"
    ## [1] "iso_7"
    ## [1] "iso_8"

``` r
tregs.WT.num.cumul <- tregs.WT.num$num.cumul
tregs.WT.num.1um <- tregs.WT.num$num.per.shell
```

### Define 3D mesh surface objects from Treg coordinate data for each sample.

This step defines 3D mesh surface objects from section or whole lymph
node images from coordinates of constituent cells. Here, Tregs are most
abundant, thus mainly used. The value of the parameter alpha in the
function, def.surface.whole.LN() is determined by trial-and-errors to
make the mesh surface object smooth enough locally but preserving the
distinct shape
globally.

``` r
tregs.TXA23.whole.LN.surf <- def.surface.whole.LN( tregs, alpha = 500) #500~1500 is reasonable
```

    ## [1] "ab_1"
    ## [1] "ab_2"
    ## [1] "ab_3"
    ## [1] "ab_4"
    ## [1] "ab_5"
    ## [1] "ab_6"
    ## [1] "ab_7"
    ## [1] "iso_1"
    ## [1] "iso_2"
    ## [1] "iso_3"
    ## [1] "iso_4"
    ## [1] "iso_5"
    ## [1] "iso_6"
    ## [1] "iso_7"
    ## [1] "iso_8"

``` r
plot3d(tregs.TXA23.whole.LN.surf$ab_2, aspect = F) # check the shape of the surface object
```

### Measuring edge-corrected volumes of shells (Takes long).

This is the most important yet most time consuming step. This obtains
edge-corrected volumes of shells. Depending on the number of reference
cells, this may take more than a day. If users are accessible to high
performance computing facilities (HPC), it is recommended to split this
step into multiple jobs and run those in parallel (HPC) if the number of
reference cells exceed several hundreds.

``` r
#Edge correected volume of sphere (Takes long)
if(0){ 
vol.surf.WT <- vol.edge.corrected.whole.LN(surf.list = tregs.TXA23.whole.LN.surf, refs.list = WT, shells = seq(0,99,1) )
vol.surf.TXA23 <- vol.edge.corrected.whole.LN(surf.list = tregs.TXA23.whole.LN.surf, refs.list = TXA23, shells = seq(0,99,1) )
}
#for the purpose of generating this markdown
vol.surf.WT <- readRDS(file = "data/vol.surf.TXA23.IL2block.WT.Rds")
vol.surf.TXA23 <- readRDS(file = "data/vol.surf.TXA23.IL2block.TXA23.Rds")

#Edge corrected volume of each shell
vol.surf.TXA23.delt <- vol.surf.TXA23
vol.surf.TXA23.delt[,2:ncol(vol.surf.TXA23)] <- vol.surf.TXA23[,2:ncol(vol.surf.TXA23)]  - vol.surf.TXA23[,1:(ncol(vol.surf.TXA23)-1)]
vol.surf.WT.delt <- vol.surf.WT
vol.surf.WT.delt[,2:ncol(vol.surf.WT)] <- vol.surf.WT[,2:ncol(vol.surf.WT)]  - vol.surf.WT[,1:(ncol(vol.surf.WT)-1)]

#Plot of edge-corrected densities
# plot(1:100, vol.surf.TXA23[1,],type = 'l')
# for(i in 2:nrow(vol.surf.TXA23)){
#   points(1:100, vol.surf.TXA23[i,],type = 'l')
# }
# 
# plot(1:100, vol.surf.WT[1,],type = 'l')
# for(i in 2:nrow(vol.surf.WT)){
#   points(1:100, vol.surf.WT[i,],type = 'l')
# }
```

### Compute shell density.

Having the number of Tregs and edge-corrected volume for each shell,
obtaining the number density per shell is straightforwad: number/volume.

``` r
##concatenate lists to matrices
tregs.TXA23.num.1um.all <- NULL
tregs.WT.num.1um.all <- NULL
idx.list <-names(tregs)
for(ln in idx.list){
  tregs.TXA23.num.1um.all <- rbind(tregs.TXA23.num.1um.all, tregs.TXA23.num.1um[[ln]] )
  tregs.WT.num.1um.all <- rbind(tregs.WT.num.1um.all, tregs.WT.num.1um[[ln]] )
}

##Shell density
shell.density.TXA23.1um <- tregs.TXA23.num.1um.all[,-1]/vol.surf.TXA23.delt
#Remove infinit values
for(i in 1:100){
  shell.density.TXA23.1um[is.nan(shell.density.TXA23.1um[,i]),i] <- 0
}

shell.density.WT.1um <- tregs.WT.num.1um.all[,-1]/vol.surf.WT.delt
#Remove infinit values
for(i in 1:100){
  shell.density.WT.1um[is.nan(shell.density.WT.1um[,i]),i] <- 0
}


sum(unlist(lapply(tregs.TXA23.num.1um, nrow)))
```

    ## [1] 33

``` r
idx.ar.ab <- 1:14
idx.ar.iso <- 15:33

sum(unlist(lapply(tregs.WT.num.1um, nrow)))
```

    ## [1] 39

``` r
idx.wt.ab <- 1:22
#remove #28 outlier
idx.wt.iso <- c(23:27,29:39)

##plot for mean
##TXA23
plot(1:100,colMeans(shell.density.TXA23.1um[idx.ar.ab,1:100]), col = "red", type = "l", lwd = 5, ylim = c(0, 0.003), xlim = c(0,90), xlab = "r, distance", ylab = "rho(r), shell density", main = "Shell density of Tregs (TXA23) at 24h: 1um, mean" )
points(1:100,colMeans(shell.density.WT.1um[idx.ar.iso,1:100]),col = "blue", type = "l", lwd = 5)
legend("topright", legend=c("ab", "iso"), pch=16, col=c("red", "blue"))
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
##WT
plot(1:100,colMeans(shell.density.WT.1um[idx.ar.ab,1:100]), col = "red", type = "l", lwd = 5, ylim = c(0, 0.003), xlim = c(0,90), xlab = "r, distance", ylab = "rho(r), shell density", main = "Shell density of Tregs (WT) at 24h: 1um, mean" )
points(1:100,colMeans(shell.density.TXA23.1um[idx.ar.iso,1:100]),col = "blue", type = "l", lwd = 5)
legend("topright", legend=c("ab", "iso"), pch=16, col=c("red", "blue"))
```

![](README_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

Densities of 1um shells are noisy due to the finite size of Tregs (~ 5um
of radius) and point estimations of coordinates.

### Obtain density with moving windows with a specified width.

Finally to mitigate noise in density estimation, numbers and volumes of
multiple adjacent shells are pooled. Considering the average radius of
Tregs being 5um, we chose the width = 5 (5 adjascent shells
pooled).

``` r
shell.density.TXA23.mov.5um <- shell.density.moving.width(cell.num.per.shell = 
                                                              tregs.TXA23.num.1um.all[,-1],
                                                           volume.per.shell = vol.surf.TXA23.delt,
                                                           shell.range = c(1,100),
                                                           width = 5
                                                           )

shell.density.WT.mov.5um <- shell.density.moving.width(cell.num.per.shell = tregs.WT.num.1um.all[,-1],
                                                           volume.per.shell = vol.surf.WT.delt,
                                                           shell.range = c(1,100),
                                                           width = 5
                                                           )
#Remove infinit values
for(i in 1:ncol(shell.density.TXA23.mov.5um)){
  shell.density.TXA23.mov.5um[is.nan(shell.density.TXA23.mov.5um[,i]),i] <- 0
  shell.density.WT.mov.5um[is.nan(shell.density.WT.mov.5um[,i]),i] <- 0
}
```

Plots

``` r
###TXA23
r <- as.numeric(colnames(shell.density.TXA23.mov.5um)) 
plot(r, shell.density.TXA23.mov.5um[1, ], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.2), ylim = c(0,0.006),ylab = "rho(r), shell density", xlab = "r, distance", main = "Shell density (TXA23) at 24h: Tregs, moving width = 5")
for(i in 2:length(idx.ar.ab)){
  points(r, shell.density.TXA23.mov.5um[idx.ar.ab[i],], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.3))
}

for(i in 1:length(idx.ar.iso)){
  points(r, shell.density.TXA23.mov.5um[idx.ar.iso[i], ], type = "l",  col =rgb(red=0, green=0, blue=1, alpha=0.3))
}
legend("topright", legend=c("ab", "iso"), pch=16, col=c("red", "blue"))
points(r, colMeans( shell.density.TXA23.mov.5um[idx.ar.ab,]), type = "l", col = "red", lwd = 5) 
points(r, colMeans(shell.density.TXA23.mov.5um[idx.ar.iso,]),  type = "l", col = "blue", lwd = 5) 
abline(v = 7.5)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
###WT
r <- as.numeric(colnames(shell.density.WT.mov.5um)) 
plot(r, shell.density.WT.mov.5um[1, ], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.2), ylim = c(0,0.006),ylab = "rho(r), shell density", xlab = "r, distance", main = "Shell density (WT) at 24h: Tregs, moving width = 5")
for(i in 2:length(idx.wt.ab)){
  points(r, shell.density.WT.mov.5um[idx.wt.ab[i],], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.3))
}

##exclude #28 from WT
for(i in 1:(length(idx.wt.iso)-1)){
  points(r, shell.density.WT.mov.5um[idx.wt.iso[i], ], type = "l",  col =rgb(red=0, green=0, blue=1, alpha=0.3))
}
legend("topright", legend=c("ab", "iso"), pch=16, col=c("red", "blue"))
points(r, colMeans( shell.density.WT.mov.5um[idx.wt.ab,]), type = "l", col = "red", lwd = 5) 
points(r, colMeans(shell.density.WT.mov.5um[idx.wt.iso,]),  type = "l", col = "blue", lwd = 5) 
abline(v = 7.5)
```

![](README_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

We also tried the width
10.

``` r
shell.density.TXA23.mov.10um <- shell.density.moving.width(cell.num.per.shell = 
                                                              tregs.TXA23.num.1um.all[,-1],
                                                           volume.per.shell = vol.surf.TXA23.delt,
                                                           shell.range = c(1,90),
                                                           width = 10
                                                           )

shell.density.WT.mov.10um <- shell.density.moving.width(cell.num.per.shell = tregs.WT.num.1um.all[,-1],
                                                           volume.per.shell = vol.surf.WT.delt,
                                                           shell.range = c(1,90),
                                                           width = 10
                                                           )
#Remove infinit values
for(i in 1:ncol(shell.density.TXA23.mov.10um)){
  shell.density.TXA23.mov.10um[is.nan(shell.density.TXA23.mov.10um[,i]),i] <- 0
  shell.density.WT.mov.10um[is.nan(shell.density.WT.mov.10um[,i]),i] <- 0
}
```

Plot

``` r
###TXA23
r <- as.numeric(colnames(shell.density.TXA23.mov.10um)) 
plot(r, shell.density.TXA23.mov.10um[1, ], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.2), ylim = c(0,0.006),ylab = "rho(r), shell density", xlab = "r, distance", main = "Shell density (TXA23) at 24h: Tregs, moving width = 5")
for(i in 2:length(idx.ar.ab)){
  points(r, shell.density.TXA23.mov.10um[idx.ar.ab[i],], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.3))
}

for(i in 1:length(idx.ar.iso)){
  points(r, shell.density.TXA23.mov.10um[idx.ar.iso[i], ], type = "l",  col =rgb(red=0, green=0, blue=1, alpha=0.3))
}
legend("topright", legend=c("ab", "iso"), pch=16, col=c("red", "blue"))
points(r, colMeans( shell.density.TXA23.mov.10um[idx.ar.ab,]), type = "l", col = "red", lwd = 5) 
points(r, colMeans(shell.density.TXA23.mov.10um[idx.ar.iso,]),  type = "l", col = "blue", lwd = 5) 
abline(v = 7.5)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
###WT
r <- as.numeric(colnames(shell.density.WT.mov.10um)) 
plot(r, shell.density.WT.mov.10um[1, ], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.2), ylim = c(0,0.006),ylab = "rho(r), shell density", xlab = "r, distance", main = "Shell density (WT) at 24h: Tregs, moving width = 5")
for(i in 2:length(idx.wt.ab)){
  points(r, shell.density.WT.mov.10um[idx.wt.ab[i],], type = "l",  col =rgb(red=1, green=0, blue=0, alpha=0.3))
}

##exclude #28 from WT
for(i in 1:(length(idx.wt.iso)-1)){
  points(r, shell.density.WT.mov.10um[idx.wt.iso[i], ], type = "l",  col =rgb(red=0, green=0, blue=1, alpha=0.3))
}
legend("topright", legend=c("ab", "iso"), pch=16, col=c("red", "blue"))
points(r, colMeans( shell.density.WT.mov.10um[idx.wt.ab,]), type = "l", col = "red", lwd = 5) 
points(r, colMeans(shell.density.WT.mov.10um[idx.wt.iso,]),  type = "l", col = "blue", lwd = 5) 
abline(v = 7.5)
```

![](README_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->
