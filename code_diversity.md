Libraries
=========

Load the libraries required.

``` r
library(vegan)
library(ade4)
library(adiv)
library(adegraphics)
library(RRPP)
```

Data
====

Set a working directory
-----------------------

This is an optional step. Researchers could be interested in saving data and results in the same folder, the working directory. If data is stored in the working directory they will call data everytime with the name of the file and it will not be necessary to write the path to a file anymore.

``` r
setwd("path_to_your_directory")
```

Databases Case 1
----------------

### Host-parasite abundance data

We upload the host-parasite abundance matrix. We subtset this matrix to perform the analysis of the Case study 1 (see main text). Case 1 includes parasite species of *Chelon auratus*, *Chelon ramada* and *Mugil cephalus* from Santa Pola Sea (SPS) autumn 2004 (2) and autumn 2005 (4).

``` r
host<- read.delim("diversity_data.txt", row.names = 1)
host1<- host[host$case<= 1,]
h1data<-host1[,1:5] # host info
h1data<- cbind(h1data, struc=as.factor(rep(c("2CSPS", "4CSPS", "2ASPS", "4ASPS",
                        "2RSPS", "4RSPS"), c(20, 30, 12, 30, 30, 15))))

# host1[remove host without parasites, remove host information]
host1<-host1[apply(host1[,6:ncol(host1)], 1, sum)>0,-(1:5)]

# remove parasite species that are not present in this Case study
host1<-host1[,apply(host1, 2, sum)>0]
```

### Parasite functional trait data

We upload the functional trait information of all the parasite species in the sample. Then, we log-transformed numerical traits.

``` r
traits<- read.delim("traits_total.txt", row.names = 1)
traits$Esiz<-log(traits$Esiz+1)
traits$Enum<-log(traits$Enum+1)
traits$Bmas<-log(traits$Bmas+1)
```

We subset the parasite species present in Case 1.

``` r
# functional traits of parasite species present in Case 1
traits1<- traits[intersect(rownames(traits), colnames(host1)),]
```

To give the same weight to all traits, we transform each column into a single matrix. Otherwise, a category of a categorical functional trait (e.g. attachment sucker) would have the same weight as a numerical functional trait (e.g. egg size).

``` r
# Give the same weight to all traits
n_atta1<-data.frame(atta=traits1[,1], row.names = row.names(traits1))
q_esiz1<-data.frame(esiz=traits1[,2], row.names = row.names(traits1))
q_enum1<-data.frame(enum=traits1[,3], row.names = row.names(traits1))
q_bmas1<-data.frame(bmas=traits1[,4], row.names = row.names(traits1))
n_lcyc1<-data.frame(lcyc=traits1[,5], row.names = row.names(traits1))
```

Then, we have to create a ktab object. It is a list of objects, each belonging to the class data.frame (see Supplementary Material in Pavoine et al., 2009, <https://doi.org/10.1111/j.1600-0706.2008.16668.x>). We calculate Gower distances between parasite species in terms of functional traits. Notice that the type of functional trait must be especified (i.e. Q: quantitative; N: nominal). Finally, following Pavoine et al. (2009), we transform the Gower distances into Euclidean distances and rescale them between 0 and 1.

``` r
ktab1<-ktab.list.df(list(q_esiz1, q_enum1, q_bmas1, n_atta1, n_lcyc1))
# Gower distance -> matrix of functional distances:
traits1dist<- dist.ktab(ktab1, c("Q", "Q", "Q", "N", "N"))
# transform Gower distances into Euclidean ones:
traits1dist<-lingoes(traits1dist)
traits1dist<-traits1dist/max(traits1dist) # [0,1]
#is.euclid(traits1dist) # double-check if the distance matrix is Euclidean
```

### Parasite proxy of phylogeny

We upload the phylogenetic-like information of all the parasite species in the sample.

``` r
phylo<- read.delim("phylo_total.txt", row.names = 1)
```

We subset the parasite species present in Case 1.

``` r
# phylogenetic-like information of parasite species present in Case 1
phylo1<- phylo[intersect(rownames(phylo), colnames(host1)),]
phylo1dist<-taxa2dist(phylo1, varstep=TRUE) # phylogenetic-like distances
phylo1dist<-phylo1dist/max(phylo1dist) #[0,1]
#is.euclid(phylo1dist)
```

Databases Case 1 sample autumn 2004
-----------------------------------

### Host-parasite abundance data

``` r
# Subset host information from autumn 2004
h1data_s2<- as.data.frame.matrix(h1data[c(1:20, 51:62, 93:122),])

# Subset community from autumn 2004
host1_s2<-host1[c(1:20, 51:62, 93:122),]

# remove parasite species that are not present in Case 1 sample autumn 2004
host1_s2<-host1_s2[,apply(host1_s2, 2, sum)>0]
```

### Parasite functional trait data

``` r
# Functional traits of species present in Case 1 autumn 2004:
traits1_s2<- traits1[intersect(rownames(traits1), colnames(host1_s2)),]

# Give the same weight to all traits:
n_atta1s2<-data.frame(atta=traits1_s2[,1], row.names = row.names(traits1_s2))
q_esiz1s2<-data.frame(esiz=traits1_s2[,2], row.names = row.names(traits1_s2))
q_enum1s2<-data.frame(enum=traits1_s2[,3], row.names = row.names(traits1_s2))
q_bmas1s2<-data.frame(bmas=traits1_s2[,4], row.names = row.names(traits1_s2))
n_lcyc1s2<-data.frame(lcyc=traits1_s2[,5], row.names = row.names(traits1_s2))
ktab1s2<-ktab.list.df(list(q_esiz1s2,q_enum1s2,q_bmas1s2,n_atta1s2,n_lcyc1s2))
traits1dist_s2<- dist.ktab(ktab1s2, c("Q", "Q", "Q", "N", "N")) # Gower distance
traits1dist_s2<-lingoes(traits1dist_s2) # Euclidean distance
traits1dist_s2<-traits1dist_s2/max(traits1dist_s2) # [0,1]
#is.euclid(phylo1dist_s2)
```

### Parasite proxy of phylogeny

``` r
# Phylogenetic-like information of species present in Case 1 autumn 2004:
phylo1_s2<-phylo1[intersect(row.names(phylo1), colnames(host1_s2)),]
phylo1dist_s2<-taxa2dist(phylo1_s2, varstep=TRUE)
phylo1dist_s2<-phylo1dist_s2/max(phylo1dist_s2)
#is.euclid(phylo1dist_s2)
```

Databases Case 1 sample autumn 2005
-----------------------------------

### Host-parasite abundance data

``` r
# Subset host information from autumn 2005:
h1data_s4<- as.data.frame.matrix(h1data[c(21:50, 63:92, 123:137),])
# Subset community from autumn 2004:
host1_s4<-host1[c(21:50, 63:92, 123:137),]
# Remove parasite species that are not present in Case 1 autumn 2005:
host1_s4<-host1_s4[,apply(host1_s4, 2, sum)>0]
```

### Parasite functional trait data

``` r
# Functional traits of species present in Case 1 autumn 2005:
traits1_s4<- traits1[intersect(rownames(traits1), colnames(host1_s4)),]
# Give the same weight to all traits:
n_atta1s4<-data.frame(atta=traits1_s4[,1], row.names = row.names(traits1_s4))
q_esiz1s4<-data.frame(esiz=traits1_s4[,2], row.names = row.names(traits1_s4))
q_enum1s4<-data.frame(enum=traits1_s4[,3], row.names = row.names(traits1_s4))
q_bmas1s4<-data.frame(bmas=traits1_s4[,4], row.names = row.names(traits1_s4))
n_lcyc1s4<-data.frame(lcyc=traits1_s4[,5], row.names = row.names(traits1_s4))
ktab1s4<-ktab.list.df(list(q_esiz1s4,q_enum1s4,q_bmas1s4,n_atta1s4,n_lcyc1s4))
traits1dist_s4<- dist.ktab(ktab1s4, c("Q", "Q", "Q", "N", "N")) # Gower distance
traits1dist_s4<-lingoes(traits1dist_s4) # Euclidean distance
traits1dist_s4<-traits1dist_s4/max(traits1dist_s4) # [0,1]
#is.euclid(phylo1dist_s4)
```

### Parasite proxy of phylogeny

``` r
# Proxy of phylogeney Case 1 sample 4:
# Phylogenetic-like information of species present in Case 1 autumn 2005:
phylo1_s4<-phylo1[intersect(row.names(phylo1), colnames(host1_s4)),]
phylo1dist_s4<-taxa2dist(phylo1_s4, varstep=TRUE)
#is.euclid(phylo1dist_s4)
phylo1dist_s4<-phylo1dist_s4/max(phylo1dist_s4)
#is.euclid(phylo1dist_s4)
```

Diversity analyses
==================

Influence of one factor on diversity
------------------------------------

### α diversity analyses

The dpcoa function is a combined version of the Rao index of diversity and the Weighted Principal Coordinate Analysis. DPCoA allows comparing the partitioning of diversity at different levels of an organisational scale and the different facets of diversity (Pavoine et al., 2004, <https://doi.org/10.1016/j.jtbi.2004.02.014>). When the argument dis is set as NULL, the funcion calculates de Taxonomic Diversity (TD) of the sample. In contrast, when a matrix of euclidean distances is given, it calculates either the Functional Diversity (FD) or the Phylogenetic Diversity (PD, or a Proxy of the Phylogenetic Diversity - PPD) of the parasite community.

#### Case 1 sample autumn 2004

``` r
dpcoa1TD_s2<- dpcoa(df = host1_s2, dis = NULL, scannf=T, RaoDecomp = T, full=T)
pco1TD_s2<- dudi.pco(dpcoa1TD_s2$RaoDis, scannf = F, full = T)
dpcoa1FD_s2<- dpcoa(host1_s2, traits1dist_s2, scannf=T, RaoDecomp = T, full=T)
pco1FD_s2<- dudi.pco(dpcoa1FD_s2$RaoDis, scannf = F, full = T)
dpcoa1PD_s2<- dpcoa(host1_s2, phylo1dist_s2, scannf=T, RaoDecomp = T, full=T)
pco1PD_s2<- dudi.pco(dpcoa1PD_s2$RaoDis, scannf = F, full = T)
```

We transform the α diversity values into their Equivalent Numbers (Ricotta and Szeidl, 2009, <https://doi.org/10.1016/j.tpb.2009.10.001>).

``` r
alphaTD1_s2<-1/(1-dpcoa1TD_s2$RaoDiv) # Ricotta and Szeild 2009
alphaFD1_s2<-1/(1-dpcoa1FD_s2$RaoDiv) # Ricotta and Szeild 2009
alphaPD1_s2<-1/(1-dpcoa1PD_s2$RaoDiv) # Ricotta and Szeild 2009
```

##### Statistical analyses

Function lm.rrpp of package RRPP (Collyer and Adams 2018a, <https://cran.r-project.org/web/packages/RRPP/index.html>) performs a linear model by residual randomisation and provides empirical sampling distributions for further ANOVAs. Following Collyer and Adams (2018b, <https://doi.org/10.1111/2041-210X.13029>), univariate α values were log-transformed. Then, we performed ANOVAs (type I of sums of squares) using random distributions of the F-statistics (Collyer and Adams, 2018b) for TD, FD and PPD, independently. When differences between samples from different host species or localities were significant, we ran a posteriori pairwise comparisons of α TD, FD and PPD between host species or localities using function pairwise in RRPP. 


-   α TD

``` r
alphaTD1_s2_st<-lm.rrpp(log(alphaTD1_s2)~h1data_s2$host.species, SS.type = "I",
                        print.progress = F)
anova(alphaTD1_s2_st, effect.type = "F")
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 1000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##                        Df      SS      MS     Rsq      F       Z Pr(>F)
    ## h1data_s2$host.species  2  0.3109 0.15547 0.02625 0.7952 0.27411  0.447
    ## Residuals              59 11.5357 0.19552 0.97375                      
    ## Total                  61 11.8466                                      
    ## 
    ## Call: lm.rrpp(f1 = log(alphaTD1_s2) ~ h1data_s2$host.species, SS.type = "I",  
    ##     print.progress = F)

-   α FD

``` r
alphaFD1_s2_st<-lm.rrpp(log(alphaFD1_s2)~h1data_s2$host.species, SS.type = "I",
                        print.progress = F)
anova(alphaFD1_s2_st, effect.type = "F")
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 1000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##                        Df      SS        MS     Rsq      F        Z Pr(>F)
    ## h1data_s2$host.species  2 0.00601 0.0030060 0.01154 0.3443 -0.39509  0.693
    ## Residuals              59 0.51518 0.0087318 0.98846                       
    ## Total                  61 0.52119                                         
    ## 
    ## Call: lm.rrpp(f1 = log(alphaFD1_s2) ~ h1data_s2$host.species, SS.type = "I",  
    ##     print.progress = F)

-   α PPD

``` r
alphaPD1_s2_st<-lm.rrpp(log(alphaPD1_s2)~h1data_s2$host.species, SS.type = "I",
                        print.progress = F)
anova(alphaPD1_s2_st, effect.type = "F")
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 1000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##                        Df      SS       MS     Rsq      F     Z Pr(>F)  
    ## h1data_s2$host.species  2 0.12865 0.064325 0.12298 4.1368 1.495  0.023 *
    ## Residuals              59 0.91743 0.015550 0.87702                      
    ## Total                  61 1.04608                                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Call: lm.rrpp(f1 = log(alphaPD1_s2) ~ h1data_s2$host.species, SS.type = "I",  
    ##     print.progress = F)

``` r
alphaPD1_s2_stpairwise<- pairwise(alphaPD1_s2_st,
                                  groups = h1data_s2$host.species)
summary(alphaPD1_s2_stpairwise, confidence = 0.95, test.type = "dist",
        stat.table = FALSE) # distances between means
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: Ca Cr Mc 
    ## 
    ## RRPP: 1000 permutations
    ## 
    ## LS means:
    ## Vectors hidden (use show.vectors = TRUE to view)
    ## 
    ## Pairwise distances between means
    ##            Ca         Cr         Mc
    ## Ca 0.00000000 0.05904259 0.04320946
    ## Cr 0.05904259 0.00000000 0.10225205
    ## Mc 0.04320946 0.10225205 0.00000000
    ## 
    ## Pairwise 95% upper confidence limits between means
    ##            Ca         Cr         Mc
    ## Ca 0.00000000 0.08391996 0.09312965
    ## Cr 0.08391996 0.00000000 0.07647935
    ## Mc 0.09312965 0.07647935 0.00000000
    ## 
    ## Pairwise effect sizes (Z) between means
    ##           Ca        Cr        Mc
    ## Ca 0.0000000 0.9282631 0.1645996
    ## Cr 0.9282631 0.0000000 3.0304877
    ## Mc 0.1645996 3.0304877 0.0000000
    ## 
    ## Pairwise P-values between means
    ##       Ca    Cr    Mc
    ## Ca 1.000 0.164 0.373
    ## Cr 0.164 1.000 0.009
    ## Mc 0.373 0.009 1.000

#### Case 1 sample autumn 2005

``` r
dpcoa1TD_s4<- dpcoa(host1_s4, dis = NULL, scannf=T, RaoDecomp = T, full=T)
pco1TD_s4<- dudi.pco(dpcoa1TD_s4$RaoDis, scannf = F, full = T)
dpcoa1FD_s4<- dpcoa(host1_s4, traits1dist_s4, scannf=T, RaoDecomp = T, full=T)
pco1FD_s4<- dudi.pco(dpcoa1FD_s4$RaoDis, scannf = F, full = T)
dpcoa1PD_s4<- dpcoa(host1_s4, phylo1dist_s4, scannf=T, RaoDecomp = T, full=T)
pco1PD_s4<- dudi.pco(dpcoa1PD_s4$RaoDis, scannf = F, full = T)
```

We transform the α diversity values into their Equivalent Numbers (Ricotta and Szeidl, 2009).

``` r
alphaTD1_s4<-1/(1-dpcoa1TD_s4$RaoDiv) # Ricotta & Szeild 2009
alphaFD1_s4<-1/(1-dpcoa1FD_s4$RaoDiv) # Ricotta & Szeild 2009
alphaPD1_s4<-1/(1-dpcoa1PD_s4$RaoDiv) # Ricotta & Szeild 2009
```

##### Statistical analyses

-   α TD

``` r
alphaTD1_s4_st<-lm.rrpp(log(alphaTD1_s4)~h1data_s4$host.species,
                        SS.type = "I", print.progress = F)
anova(alphaTD1_s4_st, effect.type = "F")
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 1000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##                        Df      SS      MS     Rsq      F     Z Pr(>F)   
    ## h1data_s4$host.species  2  4.1981 2.09907 0.25023 12.015 2.356  0.001 **
    ## Residuals              72 12.5789 0.17471 0.74977                       
    ## Total                  74 16.7771                                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Call: lm.rrpp(f1 = log(alphaTD1_s4) ~ h1data_s4$host.species, SS.type = "I",  
    ##     print.progress = F)

``` r
# Pairwise distances between the parasite community of the three host species:
alphaTD1_s4_stpairwise<- pairwise(alphaTD1_s4_st,
                                  groups = h1data_s4$host.species)
summary(alphaTD1_s4_stpairwise, confidence = 0.95, test.type = "dist",
        stat.table = FALSE) 
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: Ca Cr Mc 
    ## 
    ## RRPP: 1000 permutations
    ## 
    ## LS means:
    ## Vectors hidden (use show.vectors = TRUE to view)
    ## 
    ## Pairwise distances between means
    ##            Ca        Cr         Mc
    ## Ca 0.00000000 0.6046250 0.02794706
    ## Cr 0.60462501 0.0000000 0.57667794
    ## Mc 0.02794706 0.5766779 0.00000000
    ## 
    ## Pairwise 95% upper confidence limits between means
    ##           Ca        Cr        Mc
    ## Ca 0.0000000 0.2914843 0.2337740
    ## Cr 0.2914843 0.0000000 0.2917055
    ## Mc 0.2337740 0.2917055 0.0000000
    ## 
    ## Pairwise effect sizes (Z) between means
    ##            Ca       Cr         Mc
    ## Ca  0.0000000 5.380400 -0.9167596
    ## Cr  5.3803998 0.000000  4.9479336
    ## Mc -0.9167596 4.947934  0.0000000
    ## 
    ## Pairwise P-values between means
    ##       Ca    Cr    Mc
    ## Ca 1.000 0.001 0.809
    ## Cr 0.001 1.000 0.001
    ## Mc 0.809 0.001 1.000

-   α FD

``` r
alphaFD1_s4_st<-lm.rrpp(log(alphaFD1_s4)~h1data_s4$host.species, SS.type = "I",
                        print.progress = F)
anova(alphaFD1_s4_st, effect.type = "F")
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 1000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##                        Df     SS        MS     Rsq      F      Z Pr(>F)   
    ## h1data_s4$host.species  2 0.0622 0.0311023 0.14373 6.0426 1.7815  0.002 **
    ## Residuals              72 0.3706 0.0051472 0.85627                        
    ## Total                  74 0.4328                                          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Call: lm.rrpp(f1 = log(alphaFD1_s4) ~ h1data_s4$host.species, SS.type = "I",  
    ##     print.progress = F)

``` r
alphaFD1_s4_stpairwise<- pairwise(alphaFD1_s4_st,
                                  groups = h1data_s4$host.species)
summary(alphaFD1_s4_stpairwise, confidence = 0.95, test.type = "dist",
        stat.table = FALSE)
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: Ca Cr Mc 
    ## 
    ## RRPP: 1000 permutations
    ## 
    ## LS means:
    ## Vectors hidden (use show.vectors = TRUE to view)
    ## 
    ## Pairwise distances between means
    ##             Ca        Cr          Mc
    ## Ca 0.00000e+00 0.0719501 9.58008e-05
    ## Cr 7.19501e-02 0.0000000 7.20459e-02
    ## Mc 9.58008e-05 0.0720459 0.00000e+00
    ## 
    ## Pairwise 95% upper confidence limits between means
    ##            Ca         Cr         Mc
    ## Ca 0.00000000 0.04778893 0.03879255
    ## Cr 0.04778893 0.00000000 0.04769027
    ## Mc 0.03879255 0.04769027 0.00000000
    ## 
    ## Pairwise effect sizes (Z) between means
    ##           Ca       Cr        Mc
    ## Ca  0.000000 3.602653 -1.325722
    ## Cr  3.602653 0.000000  3.607696
    ## Mc -1.325722 3.607696  0.000000
    ## 
    ## Pairwise P-values between means
    ##       Ca    Cr    Mc
    ## Ca 1.000 0.002 0.998
    ## Cr 0.002 1.000 0.004
    ## Mc 0.998 0.004 1.000

-   α PPD

``` r
alphaPD1_s4_st<-lm.rrpp(log(alphaPD1_s4)~h1data_s4$host.species, SS.type = "I",
                        print.progress = F)
anova(alphaPD1_s4_st, effect.type = "F")
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 1000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##                        Df      SS       MS     Rsq      F      Z Pr(>F)   
    ## h1data_s4$host.species  2 0.39847 0.199236 0.42441 26.545 2.9674  0.001 **
    ## Residuals              72 0.54041 0.007506 0.57559                        
    ## Total                  74 0.93888                                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Call: lm.rrpp(f1 = log(alphaPD1_s4) ~ h1data_s4$host.species, SS.type = "I",  
    ##     print.progress = F)

``` r
alphaPD1_s4_stpairwise<- pairwise(alphaPD1_s4_st,
                                  groups = h1data_s4$host.species)
summary(alphaPD1_s4_stpairwise, confidence = 0.95, test.type = "dist",
        stat.table = FALSE)
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: Ca Cr Mc 
    ## 
    ## RRPP: 1000 permutations
    ## 
    ## LS means:
    ## Vectors hidden (use show.vectors = TRUE to view)
    ## 
    ## Pairwise distances between means
    ##            Ca        Cr         Mc
    ## Ca 0.00000000 0.1022249 0.09372099
    ## Cr 0.10222486 0.0000000 0.19594585
    ## Mc 0.09372099 0.1959459 0.00000000
    ## 
    ## Pairwise 95% upper confidence limits between means
    ##            Ca         Cr         Mc
    ## Ca 0.00000000 0.07087518 0.05601063
    ## Cr 0.07087518 0.00000000 0.06796588
    ## Mc 0.05601063 0.06796588 0.00000000
    ## 
    ## Pairwise effect sizes (Z) between means
    ##          Ca       Cr       Mc
    ## Ca 0.000000 3.444161 4.067709
    ## Cr 3.444161 0.000000 7.650752
    ## Mc 4.067709 7.650752 0.000000
    ## 
    ## Pairwise P-values between means
    ##       Ca    Cr    Mc
    ## Ca 1.000 0.005 0.001
    ## Cr 0.005 1.000 0.001
    ## Mc 0.001 0.001 1.000



```{r}
par (mfrow = c(2,3), cex.axis=0.8)
boxplot(alphaTD1_s2~h1data_s2$struc, col="grey", xlab=NULL,
        ylab = expression(paste(alpha, " diversity autumn 2004")),
        main="Taxonomic", frame.plot = FALSE, ylim=c(min(alphaTD1_s2),
        max(alphaTD1_s2)+0.3), xaxt='n')
text(1, 4.3, "a")
text(2, 4.3, "a")
text(3, 4.3, "a")

boxplot(alphaFD1_s2~h1data_s2$struc, col="grey", xlab=NULL,
        ylab = NULL,
        main="Functional", frame.plot = FALSE, ylim=c(min(alphaFD1_s2),
        max(alphaFD1_s2)+0.1), xaxt='n')
text(1, 1.4, "a")
text(2, 1.4, "a")
text(3, 1.4, "a")

boxplot(alphaPD1_s2~h1data_s2$struc, col="grey", xlab=NULL,
        ylab = NULL,
        main="Phylogenetic Proxy", frame.plot = FALSE, ylim=c(min(alphaPD1_s2),
        max(alphaPD1_s2)+0.1), xaxt='n')
text(1, 1.6, "a b")
text(2, 1.6, "b")
text(3, 1.6, "a")

splabels <- c(expression(italic("Chelon auratus"), italic("Mugil cephalus"), italic("Chelon ramada")))

boxplot(alphaTD1_s4~h1data_s4$struc, col="grey", xlab="host species",
        ylab = expression(paste(alpha, " diversity autumn 2005")),
        main=NULL, frame.plot = FALSE, ylim=c(min(alphaTD1_s4),
        max(alphaTD1_s4)+0.5), names= splabels)
text(1, 6.6, "a")
text(2, 6.6, "a")
text(3, 6.6, "b")

boxplot(alphaFD1_s4~h1data_s4$struc, col="grey", xlab="host species",
        ylab = NULL,
        main = NULL, frame.plot = FALSE, ylim=c(min(alphaFD1_s4),
        max(alphaFD1_s4)+0.15), names= splabels)
text(1, 1.4, "a")
text(2, 1.4, "a")
text(3, 1.4, "b")

boxplot(alphaPD1_s4~h1data_s4$struc, col="grey", xlab="host species",
        ylab = NULL,
        main =NULL, frame.plot = FALSE, ylim=c(min(alphaPD1_s4),
        max(alphaPD1_s4)+0.15), names= splabels)
text(1, 1.6, "a")
text(2, 1.6, "b")
text(3, 1.6, "c")
```


<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/Fig1.png">
<figcaption>
  Figure 1. Parasite a diversity in terms of taxonomy (TD), functional traits (FD) and a proxy of the phylogeny (PPD) for each host individual of each fish species. Different lowercase letters indicate significant differences between host species.
  </figcaption>
</figure>


### β diversity analyses

#### Case 1 sample autumn 2004

We calculate β1 and β2 diversities under the equivalent number approach (Ricotta and Szeidl, 2009), using the third proposition of the Rao index of diversity in Pavoine et al. (2016, <https://doi.org/10.1111/2041-210X.12591>) specifically developed for unbalanced samplings.

-   β TD

``` r
# We obtain β1 and β2 diversities. normed according to de Bello et al. 2010,
# https://doi.org/10.1111/j.1654-1103.2010.01195.x
eqrao1TDbetas2<-EqRao(host1_s2, dis=NULL,
                      structure = as.data.frame(h1data_s2$host.species,
                                  row.names= row.names(host1_s2)),
                      formula="QE", option= "normed1", wopt = "speciesab")

# β2: Inter-h1data_s2$host.species; β1: Inter-sites Intra-h1data_s2$host.species
eqrao1TDbetas2
```

    ##                                          Normed contributions to diversity
    ## Inter-h1data_s2$host.species                                     0.2721615
    ## Inter-sites Intra-h1data_s2$host.species                         0.1047394

-   β FD

``` r
eqrao1FDbetas2<-EqRao(host1_s2, dis=traits1dist_s2,
                      structure = as.data.frame(h1data_s2$host.species,
                                  row.names=row.names(host1_s2)),
                      formula="QE", option= "normed1", wopt = "speciesab")
eqrao1FDbetas2
```

    ##                                          Normed contributions to diversity
    ## Inter-h1data_s2$host.species                                    0.06919576
    ## Inter-sites Intra-h1data_s2$host.species                        0.06511553

-   β PPD

``` r
eqrao1PDbetas2<-EqRao(host1_s2, dis=phylo1dist_s2,
                      structure = as.data.frame(h1data_s2$host.species,
                                  row.names=row.names(host1_s2)),
                      formula="QE", option= "normed1", wopt = "speciesab")
eqrao1FDbetas2
```

    ##                                          Normed contributions to diversity
    ## Inter-h1data_s2$host.species                                    0.06919576
    ## Inter-sites Intra-h1data_s2$host.species                        0.06511553

##### Statistical analyses

We compared each of the TD, FD and PPD β1 and β2 diversities with 999 randomly simulated β1 and β2 values in order to establish whether the observed values significantly differ from those randomly simulated (p &lt; 0.05).

When significant, we compared observed and simulated results to determine whether the observed β1 or β2 were greater or lower than expected at random. This allows determining whether parasite communities from fish of the same species (β1) or of different fish species (β2) are more similar (the observed value is lower than simulated values) or more dissimilar (the observed value is greater than simulated values) to each other than expected by chance. Finally, we used the standardised β1 and β2 given by EqRao function (observed β – mean of randomly simulated βs/ standard deviation of randomly simulated βs) to infer if the parasite species, traits or the phylogenetic proxy are overdispersed (negative standardised β) or clustered (positive standardised β) (Head et al., 2018, <https://doi.org/10.1002/ece3.3969>) within a level of a factor (β1) or between levels of a factor (β2).

-   β TD

``` r
# This is the permutation test of β1. We obtain the stadandardised observed
# value of β1 and a significance p.value in comparison to 999 randomly generated
# β1 values. alter = "two-sided" because we do not have a priori idea;
# wopt = "speciesab" to weight the sampling units by their sum of
# parasite species’ abundances
eqrao1TDbeta1rt_s2 <- rtestEqRao(host1_s2, NULL,
                                 structure=as.data.frame(h1data_s2$host.species,
                                            row.names=row.names(host1_s2)), 
                                 formula="QE", level=1, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")

# This is the permutation test of β2. We obtain the stadandardised observed
# value of β2 and a significance p.value in comparison to 999 randomly generated
# β2 values. alter = "two-sided" because we do not have a priori idea;
# wopt = "speciesab" to weight the sampling units by their sum of
# parasite species’ abundances
eqrao1TDbeta2rt_s2 <- rtestEqRao(host1_s2, NULL,
                                 structure=as.data.frame(h1data_s2$host.species,
                                          row.names=row.names(host1_s2)), 
                                 formula="QE", level=2, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab") 
eqrao1TDbeta1rt_s2
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s2, dis = NULL, structures = as.data.frame(h1data_s2$host.species, 
    ##     row.names = row.names(host1_s2)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 1, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.1047394 
    ## 
    ## Based on 811 replicates
    ## Simulated p-value: 0.007389163 
    ## Alternative hypothesis: two-sided 
    ## 
    ##       Std.Obs   Expectation      Variance 
    ## -2.6518695324  0.1642304469  0.0005032678

``` r
eqrao1TDbeta2rt_s2
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s2, dis = NULL, structures = as.data.frame(h1data_s2$host.species, 
    ##     row.names = row.names(host1_s2)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 2, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.2721615 
    ## 
    ## Based on 999 replicates
    ## Simulated p-value: 0.001 
    ## Alternative hypothesis: two-sided 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 1.468904e+01 2.237710e-02 2.891636e-04

-   β FD

``` r
eqrao1FDbeta1rt_s2 <- rtestEqRao(host1_s2, traits1dist_s2,
                                 structure=as.data.frame(h1data_s2$host.species,
                                            row.names=row.names(host1_s2)), 
                                 formula="QE", level=1, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt = "speciesab")
eqrao1FDbeta2rt_s2 <- rtestEqRao(host1_s2, traits1dist_s2,
                                 structure=as.data.frame(h1data_s2$host.species,
                                            row.names=row.names(host1_s2)),  
                                 formula="QE", level=2, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")
eqrao1FDbeta1rt_s2
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s2, dis = traits1dist_s2, structures = as.data.frame(h1data_s2$host.species, 
    ##     row.names = row.names(host1_s2)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 1, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.06511553 
    ## 
    ## Based on 863 replicates
    ## Simulated p-value: 0.08333333 
    ## Alternative hypothesis: two-sided 
    ## 
    ##       Std.Obs   Expectation      Variance 
    ## -1.6902099767  0.0890186752  0.0001999994

``` r
eqrao1FDbeta2rt_s2
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s2, dis = traits1dist_s2, structures = as.data.frame(h1data_s2$host.species, 
    ##     row.names = row.names(host1_s2)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 2, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.06919576 
    ## 
    ## Based on 999 replicates
    ## Simulated p-value: 0.001 
    ## Alternative hypothesis: two-sided 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 1.411588e+01 6.987469e-03 1.942139e-05

-   β PPD

``` r
eqrao1PDbeta1rt_s2 <- rtestEqRao(host1_s2, phylo1dist_s2,
                                 structure=as.data.frame(h1data_s2$host.species,
                                          row.names=row.names(host1_s2)),  
                                 formula="QE", level=1, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")
eqrao1PDbeta2rt_s2 <- rtestEqRao(host1_s2, phylo1dist_s2,
                                 structure=as.data.frame(h1data_s2$host.species,
                                          row.names=row.names(host1_s2)), 
                                 formula="QE", level=2, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")
eqrao1PDbeta1rt_s2
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s2, dis = phylo1dist_s2, structures = as.data.frame(h1data_s2$host.species, 
    ##     row.names = row.names(host1_s2)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 1, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.0570059 
    ## 
    ## Based on 854 replicates
    ## Simulated p-value: 0.1508772 
    ## Alternative hypothesis: two-sided 
    ## 
    ##       Std.Obs   Expectation      Variance 
    ## -1.4364567781  0.0750786040  0.0001582927

``` r
eqrao1PDbeta2rt_s2
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s2, dis = phylo1dist_s2, structures = as.data.frame(h1data_s2$host.species, 
    ##     row.names = row.names(host1_s2)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 2, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.04934041 
    ## 
    ## Based on 999 replicates
    ## Simulated p-value: 0.001 
    ## Alternative hypothesis: two-sided 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 1.256073e+01 5.314619e-03 1.228527e-05

``` r
ADEgS(c(plot(eqrao1TDbeta1rt_s2,
        main=expression(paste("(a) Case 1, autumn 2004, TD, ",beta,"1")),
        plot=F),
        plot(eqrao1TDbeta2rt_s2,
        main=expression(paste("(d) Case 1, autumn 2004, TD, ",beta,"2")),
        plot=F),
        plot(eqrao1FDbeta1rt_s2,
        main=expression(paste("(b) Case 1, autumn 2004, FD, ",beta,"1")),
        plot=F),
        plot(eqrao1FDbeta2rt_s2,
        main=expression(paste("(e) Case 1, autumn 2004, FD, ",beta,"2")),
        plot=F),
        plot(eqrao1PDbeta1rt_s2,
        main=expression(paste("(c) Case 1, autumn 2004, PPD, ",beta,"1")),
        plot=F),
        plot(eqrao1PDbeta2rt_s2,
        main=expression(paste("(f) Case 1, autumn 2004, PPD, ",beta,"2")),
        plot=F)),
      layout = c(3, 2))
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/Fig2.png">
<figcaption>
  Figure 2. Observed and simulated β diversity values (Case 1: autumn 2004). (a, b, c) β1 diversity or extent of dissimilarity in the diversity of parasite communities among host individuals within each host species (Chelon auratus, Mugil cephalus and Chelon ramada). (d, e, f) β2 diversity or extent of dissimilarity in the diversity of parasite communities between host species. Diversity was measured in terms of (a, d) Taxonomic Diversity (TD), (b, e) Functional Diversity (FD) and (c, f) the Proxy of Phylogenetic Diversity (PPD). Samples are from Santa Pola Lagoon and autumn 2004 (Case 1). Observed β values (black diamond on the top of the black vertical line) and distribution of the simulated (x-axis: sim) β values (grey bars).
</figcaption>
</figure>


#### Case 1 sample autumn 2005

We calculate β1 and β2 diversities under the equivalent number approach (Ricotta and Szeidl, 2009) using the third proposition of the Rao index of diversity in Pavoine et al. (2016, <https://doi.org/10.1111/2041-210X.12591>) since it is specifically developed for unbalanced samplings.

-   β TD

``` r
eqrao1TDbetas4<-EqRao(host1_s4, dis=NULL,
                      structure=as.data.frame(h1data_s4$host.species,
                                              row.names=row.names(host1_s4)),
                      formula="QE", option= "normed1", wopt = "speciesab")
eqrao1TDbetas4
```

    ##                                          Normed contributions to diversity
    ## Inter-h1data_s4$host.species                                     0.7817061
    ## Inter-sites Intra-h1data_s4$host.species                         0.1646932

-   β FD

``` r
eqrao1FDbetas4<-EqRao(host1_s4, dis=traits1dist_s4,
                      structure=as.data.frame(h1data_s4$host.species,
                                              row.names=row.names(host1_s4)),
                      formula="QE", option= "normed1", wopt = "speciesab")
eqrao1FDbetas4
```

    ##                                          Normed contributions to diversity
    ## Inter-h1data_s4$host.species                                     0.1626783
    ## Inter-sites Intra-h1data_s4$host.species                         0.0559430

-   β PPD

``` r
eqrao1PDbetas4<-EqRao(host1_s4, dis=phylo1dist_s4,
                      structure=as.data.frame(h1data_s4$host.species,
                                              row.names=row.names(host1_s4)),
                      formula="QE", option= "normed1", wopt = "speciesab")
eqrao1PDbetas4
```

    ##                                          Normed contributions to diversity
    ## Inter-h1data_s4$host.species                                    0.15635556
    ## Inter-sites Intra-h1data_s4$host.species                        0.04457452

##### Statistical analyses

-   β TD

``` r
eqrao1TDbeta1rt_s4 <- rtestEqRao(host1_s4, NULL,
                                 structure=as.data.frame(h1data_s4$host.species,
                                          row.names=row.names(host1_s4)),
                                 formula = "QE", level=1, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")
eqrao1TDbeta2rt_s4 <- rtestEqRao(host1_s4, NULL,
                                 structure=as.data.frame(h1data_s4$host.species,
                                          row.names=row.names(host1_s4)), 
                                 formula="QE", level=2, nrep=999,
                                 option="normed1", alter = "two-sided",
                                 wopt = "speciesab")
eqrao1TDbeta1rt_s4
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s4, dis = NULL, structures = as.data.frame(h1data_s4$host.species, 
    ##     row.names = row.names(host1_s4)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 1, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.1646932 
    ## 
    ## Based on 978 replicates
    ## Simulated p-value: 0.8283963 
    ## Alternative hypothesis: two-sided 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## -0.370478690  0.183476858  0.002570586

``` r
eqrao1TDbeta2rt_s4
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s4, dis = NULL, structures = as.data.frame(h1data_s4$host.species, 
    ##     row.names = row.names(host1_s4)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 2, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.7817061 
    ## 
    ## Based on 999 replicates
    ## Simulated p-value: 0.001 
    ## Alternative hypothesis: two-sided 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 14.154195687  0.074920066  0.002493478

-   β FD

``` r
eqrao1FDbeta1rt_s4 <- rtestEqRao(host1_s4, traits1dist_s4,
                                 structure=as.data.frame(h1data_s4$host.species,
                                          row.names=row.names(host1_s4)), 
                                 formula = "QE", level=1, nrep=999,
                                 option="normed1", alter = "two-sided",
                                 wopt = "speciesab")
eqrao1FDbeta2rt_s4 <- rtestEqRao(host1_s4, traits1dist_s4,
                                 structure=as.data.frame(h1data_s4$host.species,
                                          row.names=row.names(host1_s4)), 
                                 formula="QE", level=2, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")
eqrao1FDbeta1rt_s4
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s4, dis = traits1dist_s4, structures = as.data.frame(h1data_s4$host.species, 
    ##     row.names = row.names(host1_s4)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 1, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.055943 
    ## 
    ## Based on 968 replicates
    ## Simulated p-value: 0.5923633 
    ## Alternative hypothesis: two-sided 
    ## 
    ##       Std.Obs   Expectation      Variance 
    ## -0.6419653973  0.0659364088  0.0002423284

``` r
eqrao1FDbeta2rt_s4
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s4, dis = traits1dist_s4, structures = as.data.frame(h1data_s4$host.species, 
    ##     row.names = row.names(host1_s4)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 2, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.1626783 
    ## 
    ## Based on 999 replicates
    ## Simulated p-value: 0.001 
    ## Alternative hypothesis: two-sided 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 1.001596e+01 2.019781e-02 2.023607e-04

-   β PPD

``` r
eqrao1PDbeta1rt_s4 <- rtestEqRao(host1_s4, phylo1dist_s4,
                                 structure=as.data.frame(h1data_s4$host.species,
                                          row.names=row.names(host1_s4)),
                                 formula="QE", level=1, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")
eqrao1PDbeta2rt_s4 <- rtestEqRao(host1_s4, phylo1dist_s4,
                                 structure=as.data.frame(h1data_s4$host.species,
                                          row.names=row.names(host1_s4)),
                                 formula="QE", level=2, nrep=999,
                                 option="normed1", alter="two-sided",
                                 wopt="speciesab")
eqrao1PDbeta1rt_s4
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s4, dis = phylo1dist_s4, structures = as.data.frame(h1data_s4$host.species, 
    ##     row.names = row.names(host1_s4)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 1, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.04457452 
    ## 
    ## Based on 975 replicates
    ## Simulated p-value: 0.4815574 
    ## Alternative hypothesis: two-sided 
    ## 
    ##       Std.Obs   Expectation      Variance 
    ## -0.7611700450  0.0540092932  0.0001536385

``` r
eqrao1PDbeta2rt_s4
```

    ## Monte-Carlo test
    ## Call: rtestEqRao(comm = host1_s4, dis = phylo1dist_s4, structures = as.data.frame(h1data_s4$host.species, 
    ##     row.names = row.names(host1_s4)), formula = "QE", option = "normed1", 
    ##     wopt = "speciesab", level = 2, nrep = 999, alter = "two-sided")
    ## 
    ## Observation: 0.1563556 
    ## 
    ## Based on 999 replicates
    ## Simulated p-value: 0.001 
    ## Alternative hypothesis: two-sided 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 1.106049e+01 1.826557e-02 1.558748e-04

``` r
ADEgS(c(plot(eqrao1TDbeta1rt_s4,
             main=expression(paste("(a) Case 1, autumn 2005, TD, ",beta,"1")),
             plot=F),
        plot(eqrao1TDbeta2rt_s4,
             main=expression(paste("(d) Case 1, autumn 2005, TD, ",beta,"2")),
             plot=F),
        plot(eqrao1FDbeta1rt_s4,
             main=expression(paste("(b) Case 1, autumn 2005, FD, ",beta,"1")),
             plot=F),
        plot(eqrao1FDbeta2rt_s4,
             main=expression(paste("(e) Case 1, autumn 2005, FD, ",beta,"2")),
             plot=F),
        plot(eqrao1PDbeta1rt_s4,
             main=expression(paste("(c) Case 1, autumn 2005, PPD, ",beta,"1")),
             plot=F),
        plot(eqrao1PDbeta2rt_s4,
             main=expression(paste("(f) Case 1, autumn 2005, PPD, ",beta,"2")),
             plot=F)),
      layout = c(3, 2))
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/Fig3.png">
<figcaption>
  Figure 3. Observed and simulated β diversity values (Case 1: autumn 2005). (a, b, c) β1 diversity or extent of dissimilarity in the diversity of parasite communities among host individuals within each host species (Chelon auratus, Mugil cephalus and Chelon ramada). (d, e, f) β2 diversity or extent of dissimilarity in the diversity of parasite communities between host species. Diversity was measured in terms of (a, d) Taxonomic Diversity (TD), (b, e) Functional Diversity (FD) and (c, f) the Proxy of Phylogenetic Diversity (PPD). Samples are from Santa Pola Lagoon and autumn 2005 (Case 1). Observed β values (black diamond on the top of the black vertical line) and distribution of the simulated (x-axis: sim) β values (grey bars).
</figcaption>
</figure>

Influence of two-crossed factors on diversity
---------------------------------------------

To evaluate and disentangle the effect of crossed factors on diversity, we use the crossed-DPCoA (Pavoine et al., 2013, <https://doi.org/10.1371/journal.pone.0054530>). We analyse the effect of two-crossed factors simultaneously: host species and season. The crossed-DPCoA is grounded on the DPCoA but it analyses the effect of two-crossed factors at the same time. Thus, it distinguishes the proportional contribution of the sampling unit, each factor individually and the effect of the interaction of both factors on the diversity of the community.

### Main effect A

The crossed-DPCoA consists of three consecutive analyses. Following Pavoine et al.’s (2013) terminology, each parasite community is associated with a component of the factor A (hosts species) and a component of a factor B (season). The main version of the crossed-DPCoA plots in a DPCoA space the parasite species, the sampling units and the variables of the main factor A, without taking into account seasonal differences (factor B).

-   crossed-DPCoA TD

``` r
ab1main<- crossdpcoa_maineffect(host1, h1data$host.species,
                                as.factor(h1data$year), dis =NULL, scannf = F)

# Positions of the levels of factor A - hosts species - on its principal axes
# in the main TD space:
s.label(ab1main$l2)
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoA_TD.png">
<figcaption>
  Figure 4. To view the position of host species in the TD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figure>


``` r
# Percentange of variation explained by each axes: most of the variation is
# reflected by the first axis:
ab1main$eig[1:2]/sum(ab1main$eig)
```

    ## [1] 0.8508505 0.1491495

``` r
# Percentage of diversity associated with each factor:
(ab1main$div*100)/ab1main$div[5]

```

    ##        SSW        SSA        SSB       SSAB        SST 
    ##  50.330382  31.235617   1.622172  16.811829 100.000000

``` r
# ssw: diversity within communities or percentaje of diversity due to
# host individual factor
# ssa: diversity related to factor A or host species;
# ssb: diversity related to factor B or season;
# ssab: diversity due to the interaction of factors or
# host species x season
# sst: the total diversity over all communities
```

-   crossed-DPCoA FD

``` r
ab1mainFD<-crossdpcoa_maineffect(host1, h1data$host.species,
                                as.factor(h1data$year),
                                dis=traits1dist, scannf=F)

# Positions of the levels of factor A - hosts species -
# on its principal axes in the FD space:
s.label(ab1mainFD$l2)
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoA_FD.png">
<figcaption>
  Figure 5. To view the position of host species in the FD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figcaption>
</figure>

``` r
# Percentange of variation explained by each axes: most of the variation
# is reflected by the first axis:
ab1mainFD$eig[1:2]/sum(ab1mainFD$eig)
```

    ## [1] 0.93234594 0.06765406

``` r
# Percentage of diversity associated with each factor:
(ab1mainFD$div*100)/ab1mainFD$div[5]
```

    ##        SSW        SSA        SSB       SSAB        SST 
    ##  50.197098  28.021099   2.058065  19.723738 100.000000

-   crossed-DPCoA PPD

``` r
ab1mainPD<- crossdpcoa_maineffect(host1, h1data$host.species,
                                  as.factor(h1data$year),
                                  dis=phylo1dist, scannf=F)

# Positions of the levels of factor A - hosts species -
# on its principal axes in the FD space
s.label(ab1mainFD$l2)
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoA_PPD.png">
<figcaption>
  Figure 6. To view the position of host species in the PPD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figure>

``` r
# Percentange of variation explained by each axes:
ab1mainPD$eig[1:2]/sum(ab1mainPD$eig)
```

    ## [1] 0.8598161 0.1401839

``` r
# Percentage of diversity associated with each factor:
(ab1mainPD$div*100)/ab1mainPD$div[5]
```

    ##        SSW        SSA        SSB       SSAB        SST 
    ##  56.345906  20.649526   1.447008  21.557561 100.000000


### crossed-DPCoA version 1

The first version of the crossed-DPCoA removes the amount of diversity among sampling units due to the sole effect of factor B, but retains combined effects of factors A and B (i.e. the interaction between factors A and B). This first version should always be performed to confirm that the two factors are not independent, otherwise versions 1 and 2 of the crossed-DPCoA will produce the same results. This can be done by visually comparing ellipse plots generated from version 1 and version 2.

-   crossed-DPCoA TD

``` r
a1v1 <- crossdpcoa_version1(host1, h1data$host.species, as.factor(h1data$year),
                            dis =NULL, w = "classic", scannf = F)

# Positions of the levels of factor A - hosts species - on its principal axes
# in the Version 1 TD space
s.class(a1v1$l3, h1data$host.species, 
        col=c("indianred", "darkgoldenrod1", "cadetblue3"))
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoAv1_TD.png">
<figcaption>
  Figure 7. Position of all communities (host individuals) in the crossed-DPCoA version 1 TD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figure>


-   crossed-DPCoA FD

``` r
a1v1FD <- crossdpcoa_version1(host1, h1data$host.species, as.factor(h1data$year),
                              dis =traits1dist, w = "classic", scannf = F)

# Positions of the levels of factor A - hosts species - on its principal axes
# in the Version 1 FD space
s.class(a1v1FD$l3, h1data$host.species, 
        col=c("indianred", "darkgoldenrod1", "cadetblue3"))
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoAv1_FD.png">
<figcaption>
  Figure 8. Position of all communities (host individuals) in the crossed-DPCoA version 1 FD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figure>


-   crossed-DPCoA PPD

``` r
a1v1PD <- crossdpcoa_version1(host1, h1data$host.species, as.factor(h1data$year),
                              dis =phylo1dist, w = "classic", scannf = F)

# Positions of the levels of factor A - hosts species -
# on its principal axes in the Version 1 PPD space
s.class(a1v1PD$l3, h1data$host.species, 
        col=c("indianred", "darkgoldenrod1", "cadetblue3"))
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoAv1_PPD.png">
<figcaption>
  Figure 9. Position of all communities (host individuals) in the crossed-DPCoA version 1 PPD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figure>


### crossed-DPCoA version 2

Finally, the second version of the crossed-DPCoA eliminates any influence of the factor B on the factor A (including the interaction term). Thus, it provides diversity exclusively under the light of the main crossed-factor, factor A (host species).

-   crossed-DPCoA TD

``` r
a1v2 <- crossdpcoa_version2(host1, h1data$host.species, as.factor(h1data$year),
                            dis =NULL, w = "classic", scannf = F)

# Positions of the levels of factor A - hosts species -
# on its principal axes in the Version 2 TD space
s.class(a1v2$l3, h1data$host.species, 
        col=c("indianred", "darkgoldenrod1", "cadetblue3"))
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoAv2_TD.png">
<figcaption>
  Figure 10. Position of all communities (host individuals) in the crossed-DPCoA version 2 TD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figure>

-   crossed-DPCoA FD

``` r
a1v2FD <- crossdpcoa_version2(host1, h1data$host.species, as.factor(h1data$year),
                              dis =traits1dist, w = "classic", scannf = F)

# Positions of the levels of factor A - hosts species -
# on its principal axes in the Version 2 FD space
s.class(a1v2FD$l3, h1data$host.species, 
        col=c("indianred", "darkgoldenrod1", "cadetblue3")) 
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoAv2_FD.png">
<figcaption>
  Figure 11. Position of all communities (host individuals) in the crossed-DPCoA version 2 FD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.
</figcaption>
</figure>

-   crossed-DPCoA PPD

``` r
a1v2PD <- crossdpcoa_version2(host1, h1data$host.species, as.factor(h1data$year),
                              dis =phylo1dist, w = "classic", scannf = F)

# Positions of the levels of factor A - hosts species -
# on its principal axes in the Version 2 PPD space
s.class(a1v2PD$l3, h1data$host.species,
        col=c("indianred", "darkgoldenrod1", "cadetblue3"))
```

<figure>
<img src="https://github.com/crisLB/diversity/blob/master/figures/crossedDPCoAv2_PPD.png">
<figcaption>
  Figure 12. Position of all communities (host individuals) in the crossed-DPCoA version 2 PPD space. Ca: Chelon auratus; Mc: Mugil cephalus; Cr: Chelon ramada. d (top-right) provides the scale.</figcaption>
</figure>
