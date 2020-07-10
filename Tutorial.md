Description
================
As described in the [`README.md`](https://github.com/danieledurante/ESBM/blob/master/README.md) file, this tutorial contains guidelines and code to perform the analyses for the first scenario [`network 1`] considered in the simulation study of the article [**extended stochastic block models**](https://github.com/danieledurante/ESBM). In particular, you will find a detailed step-by-step guide and `R` code to **implement the collapsed Gibbs-sampler presented in the article** and to **perform posterior inference** under **ESBM**. For implementation purposes, **execute the code below considering the same order in which is presented**.

Upload the data and check performance of algorithmic methods
================
To start the analysis, **set the working directory** where the simulated network [`network_1.RData`](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/network_1.RData), and the source codes [`esbm.R`](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/esbm.R) and [`stirling.cpp`](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/stirling.cpp) are placed. Once this has been done, **clean the workspace, and load the data along with useful** `R` **packages**.

``` r
rm(list=ls())
source("esbm.R")
Rcpp::sourceCpp('stirling.cpp')

library(reshape)
library(gdata)
library(igraph)
library(mcclust.ext)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(coda)
library(dummies)

load("network_1.RData")
V <- dim(Y)[1]

# Note that Y must have diagonal equal to 0
diag(Y)
```

As discussed in the article, the network under analysis has *V=100* and *5* equally–sized clusters of *20* nodes each, displaying not only classical community structures, but also core-periphery patterns [see Figure 2 in the article]. Hence, before implementing the proposed **ESBM**, it is useful to **check the performance of state-of-the-art algorithmic approaches for community detection** in this illustrative network. To accomplish this goal, let us focus on the popular [`Louvain algorithm`](https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/meta).

``` r
net <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)
class(net)
cluster_louvain(net)
```

According to the above results, the [`Louvain algorithm`](https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/meta) detects only *3* communities, wrongly collapsing clusters *2* and *3* as well as *4* and *5*. In fact, as shown Figure 2 in the article, such communities display core-periphery patterns that cannot be detected by classical algorithmic approaches. This motivates our focus on the **ESBM** class.


Setting the hyperparameters
================
For the hyperparameters of the `Beta(a,b)` **priors on the block probabilities** we follow common implementations of stochastic block models and consider the default values `a=1` and `b=1` to induce a **uniform** prior. Less straightforward is instead the choice of the hyperparameters for the **Gibb-type priors on the random partition**. A possibility to address this issue is to specify such quantities in order to obtain a value for the expectation of the non-empty number of communities `H` that matches some prior knowledge. Below, we provide the **code to obtain such a quantity as a function of pre-specified hyperparameters for the four relevant examples of Gibbs-type priors** discussed in the article. 

``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
sigma_dm <- 0   
H_dm <- 50 # Conservative upper bound 
beta_dm <- 3/H_dm 
expected_cl_py(V, sigma = sigma_dm, theta = beta_dm*H_dm, H = H_dm)

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
sigma_dp <- 0   
H_dp <- Inf 
alpha_dp <- 2.55
expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp)

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
sigma_py <- 0.575
H_py <- Inf 
alpha_py <- -0.325
expected_cl_py(V, sigma = sigma_py, theta = alpha_py, H = H_py)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------
gamma <- 0.475
probs_gnedin <- HGnedin(V, 1:V, gamma = gamma)
sum(1:V*probs_gnedin)
```

Here we **set the hyperparameters so that the expectation of `H` is close to *10* under all the four priors**. This is twice as many as the true number of communities to check whether our results are robust to hyperparameters settings.

Posterior computation via collapsed Gibbs sampler
================
This section contains the code to **implement collapsed Gibbs sampler for ESBM** [function `esbm()`] and to **evaluate marginal likelihoods** [function `log_pY_z()`] for model selection. Such a code is applied to the four relevant examples of Gibbs-type priors discussed in the article, both without and with node-specific attributes. See the source code [`esbm.R`]() for a detailed description of the inputs and the outputs of the two functions `esbm()` and `log_pY_z()`.

Implementation without node-specific attributes
------------------
Let us first **perform posterior computation for the model without node-attributes**. To accomplish this goal, execute the code below.

``` r
N_iter <- 20000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

my_prior <- "DM"
Z_DM <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 3/50, H_DM = 50)

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

my_prior <- "DP"
Z_DP <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = 2.55, sigma_PY = 0)

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

my_prior <- "PY"
Z_PY <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = -0.325, sigma_PY = 0.575)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

my_prior  <- "GN"
Z_GN <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.475)
```

Once the above steps have been done, **compute the log-marginal likelihoods** that will be used for comparing the performance of the different prior specifications, and **save the output** in the file `Posterior_No_Attributes.RData`.

``` r
# Compute the log-marginal likelihoods under the different priors

l_y_DM <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_DM[t] <- log_pY_z(Y,Z_DM[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

l_y_DP <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_DP[t] <- log_pY_z(Y,Z_DP[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

l_y_PY <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_PY[t] <- log_pY_z(Y,Z_PY[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

l_y_GN <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_GN[t] <- log_pY_z(Y,Z_GN[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

# Save the output
save(Z_DP,l_y_DP,Z_PY,l_y_PY,Z_GN,l_y_GN,Z_DM,l_y_DM,file="Posterior_No_Attributes.RData")
rm(Z_DP,l_y_DP,Z_PY,l_y_PY,Z_GN,l_y_GN,Z_DM,l_y_DM)
```

Implementation with node-specific attributes
------------------
Let us now **perform posterior computation for the model with node-specific attributes** coinciding, in this case, with the true membership labels. To accomplish this goal, execute the code below.

``` r
N_iter  <- 20000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

my_x <- c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))
my_alpha_xi <- rep(1,5)

# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

my_prior <- "DM"
Z_DM_x <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 3/50, H_DM = 50, x = my_x, alpha_xi = my_alpha_xi)

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

my_prior <- "DP"
Z_DP_x <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = 2.55, sigma_PY = 0, x = my_x, alpha_xi = my_alpha_xi)

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

my_prior <- "PY"
Z_PY_x <- esbm(Y,my_seed,N_iter, my_prior, my_z,a=1,b=1,alpha_PY=-0.325,sigma_PY=0.575,x=my_x,alpha_xi=my_alpha_xi)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

my_prior <- "GN"
Z_GN_x <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.475, x = my_x, alpha_xi = my_alpha_xi)
```

Also in this case we **compute the log-marginal likelihoods** that will be used for comparing the performance of the different prior specifications, and **save the output** in the file `Posterior_Attributes.RData`.

``` r
# Compute the log-marginal likelihoods under the different priors

l_y_DM_x <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_DM_x[t] <- log_pY_z(Y,Z_DM_x[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

l_y_DP_x <- rep(0,N_iter)
for (t in 1:N_iter){
l_y_DP_x[t] <-	log_pY_z(Y,Z_DP_x[,t],1,1)
if (t%%1000 == 0){print(paste("Iteration:", t))}
}

l_y_PY_x <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_PY_x[t] <- log_pY_z(Y,Z_PY_x[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

l_y_GN_x <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_GN_x[t] <- log_pY_z(Y,Z_GN_x[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

# Save the output

save(Z_DP_x,l_y_DP_x,Z_PY_x,l_y_PY_x,Z_GN_x,l_y_GN_x,Z_DM_x,l_y_DM_x,file="Posterior_Attributes.RData")
rm(Z_DP_x,l_y_DP_x,Z_PY_x,l_y_PY_x,Z_GN_x,l_y_GN_x,Z_DM_x,l_y_DM_x)
```

Posterior inference
================
This section contains the **code to perform estimation, uncertainty quantification and model selection for ESBM** leveraging the samples from the collapsed Gibbs sampler. In particular, we **reproduce the analyses in Table 2 and Figure 2 in the article**. To accomplish this goal let us first **upload the MCMC samples**, and define the `burn_in` along with the vector `z_0` containing the true community labels. 

``` r
burn_in <- 5000
z_0 <- c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20))
load("Posterior_No_Attributes.RData")
load("Posterior_Attributes.RData")
```
Before performing posterior inference, let us **visualize the traceplots for the logarithm of the likelihood in Eq. [3]**, evaluated at the MCMC samples of `z` under the different priors, both with and without nodal attributes.

``` r
traceplot <- melt(cbind(l_y_DM,l_y_DM_x,l_y_DP,l_y_DP_x,l_y_PY,l_y_PY_x,l_y_GN,l_y_GN_x))
traceplot <- traceplot[,-2]

traceplot$Group <- c(rep("DM",N_iter*2),rep("DP",N_iter*2),rep("PY",N_iter*2),rep("GN",N_iter*2))
traceplot$Group <- factor(traceplot$Group,levels=c("DM","DP","PY","GN"))

traceplot$Attr<- rep(c(rep("Without Attributes",20000),rep("With Attributes",20000)),4)
traceplot$Attr <- factor(traceplot$Attr,levels=c("Without Attributes","With Attributes"))

Trace <- ggplot(traceplot,aes(y=value,x=X1))+geom_line()+facet_grid(Attr~Group)+ theme_bw()+labs(y="",x="")
ggsave("Trace.png",width=10,height=4.5)
```
![](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/Trace.png)

The above traceplots confirm that our Gibbs sampler has **satisfactory mixing and rapid convergence**. Due to the high stability of the chains for the quantity in Eq. [3], we can reliably compute the **logarithm of the marginal likelihoods** for the different priors and models [without and with attributes] via the harmonic mean approach in Eq. [17].

``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

l_y_DM <- l_y_DM[(burn_in+1):N_iter]
neg_l_y_DM <- -c(l_y_DM)
l_y_post_DM <- log(length(l_y_DM))-max(neg_l_y_DM)-log(sum(exp(neg_l_y_DM-max(neg_l_y_DM))))
l_y_post_DM

l_y_DM_x <- l_y_DM_x[(burn_in+1):N_iter]
neg_l_y_DM_x <- -c(l_y_DM_x)
l_y_post_DM_x <- log(length(l_y_DM_x))-max(neg_l_y_DM_x)-log(sum(exp(neg_l_y_DM_x-max(neg_l_y_DM_x))))
l_y_post_DM_x

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

l_y_DP <- l_y_DP[(burn_in+1):N_iter]
neg_l_y_DP <- -c(l_y_DP)
l_y_post_DP <- log(length(l_y_DP))-max(neg_l_y_DP)-log(sum(exp(neg_l_y_DP-max(neg_l_y_DP))))
l_y_post_DP

l_y_DP_x <- l_y_DP_x[(burn_in+1):N_iter]
neg_l_y_DP_x <- -c(l_y_DP_x)
l_y_post_DP_x <- log(length(l_y_DP_x))-max(neg_l_y_DP_x)-log(sum(exp(neg_l_y_DP_x-max(neg_l_y_DP_x))))
l_y_post_DP_x

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

l_y_PY <- l_y_PY[(burn_in+1):N_iter]
neg_l_y_PY <- -c(l_y_PY)
l_y_post_PY <- log(length(l_y_PY))-max(neg_l_y_PY)-log(sum(exp(neg_l_y_PY-max(neg_l_y_PY))))
l_y_post_PY

l_y_PY_x <- l_y_PY_x[(burn_in+1):N_iter]
neg_l_y_PY_x <- -c(l_y_PY_x)
l_y_post_PY_x <- log(length(l_y_PY_x))-max(neg_l_y_PY_x)-log(sum(exp(neg_l_y_PY_x-max(neg_l_y_PY_x))))
l_y_post_PY_x

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

l_y_GN <- l_y_GN[(burn_in+1):N_iter]
neg_l_y_GN <- -c(l_y_GN)
l_y_post_GN <- log(length(l_y_GN))-max(neg_l_y_GN)-log(sum(exp(neg_l_y_GN-max(neg_l_y_GN))))
l_y_post_GN

l_y_GN_x <- l_y_GN_x[(burn_in+1):N_iter]
neg_l_y_GN_x <- -c(l_y_GN_x)
l_y_post_GN_x <- log(length(l_y_GN_x))-max(neg_l_y_GN_x)-log(sum(exp(neg_l_y_GN_x-max(neg_l_y_GN_x))))
l_y_post_GN_x
```

As it can be noticed, the **Gnedin process tends to perform slightly better** in both scenarios relative to the other priors. Moreover, as expected, the overall **learning process benefits from informative node-specific attributes**.   

The **posterior mean of the VI distance from the true partition `z_0`** can be instead obtained using the `VI()` function within the `mcclust.ext` package ([Wade and Ghahramani, 2018](https://projecteuclid.org/euclid.ba/1508378464)) as follow.


``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
VI(z_0,t(Z_DM[,(burn_in+1):N_iter]))
VI(z_0,t(Z_DM_x[,(burn_in+1):N_iter]))

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
VI(z_0,t(Z_DP[,(burn_in+1):N_iter]))
VI(z_0,t(Z_DP_x[,(burn_in+1):N_iter]))

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
VI(z_0,t(Z_PY[,(burn_in+1):N_iter]))
VI(z_0,t(Z_PY_x[,(burn_in+1):N_iter]))

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------
VI(z_0,t(Z_GN[,(burn_in+1):N_iter]))
VI(z_0,t(Z_GN_x[,(burn_in+1):N_iter]))
```

The above results **confirm the results and rankings** obtained from the analysis of the marginal likelihoods.

As discussed in the article, accurate learning of the underlying number of communities is a fundamental goal. Hence, let us study the **quantiles of the posterior distribution for the number of non-empty communities** under the different priors and models. 

``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
quantile(apply(Z_DM[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_DM_x[,(burn_in+1):N_iter],2,max))[c(2:4)]

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
quantile(apply(Z_DP[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_DP_x[,(burn_in+1):N_iter],2,max))[c(2:4)]

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
quantile(apply(Z_PY[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_PY_x[,(burn_in+1):N_iter],2,max))[c(2:4)]

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------
quantile(apply(Z_GN[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_GN_x[,(burn_in+1):N_iter],2,max))[c(2:4)]
```

Also for this measure, the inclusion of **informative node-specific attributes provides improved performance**. It also interesting to notice how, unlike DM, DP and PY, the **Gnedin process can learn accurately the true number of underlying communities even without the additional information** provided by the node-specific attributes.

To conclude our analysis, let us obtain a **point estimates** and **credible balls** for the community assignments of the different nodes. This is done by adapting the methods presented in [Wade and Ghahramani (2018)](https://projecteuclid.org/euclid.ba/1508378464) and implemented in the `R` package `mcclust.ext`. To apply these strategies we also require an estimate of the **co-clustering matrix**, whose generic element `c[v,u]` encodes the relative frequency of MCMC samples in which nodes `v` and `u` are in the same community. Such an estimate can be obtained via the function `pr_cc()` in the source code `esbm.R`. We also study the **misclassification error** using the function `misclass()` in the source code `esbm.R`.


``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

c_Z_DM <- pr_cc(Z_DM[,(burn_in+1):N_iter])
# Point estimate
memb_Z_DM_VI <- minVI(c_Z_DM,method="avg",max.k=20)
memb_Z_DM <- memb_Z_DM_VI$cl
VI(z_0,t(memb_Z_DM))
# Horizontal bound of the credible ball
credibleball(memb_Z_DM_VI$cl,t(Z_DM[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
misclass(memb_Z_DM,Y,a=1,b=1)


c_Z_DM <- pr_cc(Z_DM_x[,(burn_in+1):N_iter])
# Point estimate
memb_Z_DM_VI <- minVI(c_Z_DM,method="avg",max.k=20)
memb_Z_DM  <- memb_Z_DM_VI$cl
VI(z_0,t(memb_Z_DM))
# Horizontal bound of the credible ball
credibleball(memb_Z_DM_VI$cl,t(Z_DM_x[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
misclass(memb_Z_DM,Y,a=1,b=1)


# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

c_Z_DP <- pr_cc(Z_DP[,(burn_in+1):N_iter])
# Point estimate
memb_Z_DP_VI <- minVI(c_Z_DP,method="avg",max.k=20)
memb_Z_DP <- memb_Z_DP_VI$cl
VI(z_0,t(memb_Z_DP))
# Horizontal bound of the credible ball
credibleball(memb_Z_DP_VI$cl,t(Z_DP[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
misclass(memb_Z_DP,Y,a=1,b=1)


c_Z_DP <- pr_cc(Z_DP_x[,(burn_in+1):N_iter])
# Point estimate
memb_Z_DP_VI <- minVI(c_Z_DP,method="avg",max.k=20)
memb_Z_DP <- memb_Z_DP_VI$cl
VI(z_0,t(memb_Z_DP))
# Horizontal bound of the credible ball
credibleball(memb_Z_DP_VI$cl,t(Z_DP_x[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
misclass(memb_Z_DP,Y,a=1,b=1)


# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

c_Z_PY <- pr_cc(Z_PY[,(burn_in+1):N_iter])
# Point estimate
memb_Z_PY_VI <- minVI(c_Z_PY,method="avg",max.k=20)
memb_Z_PY <- memb_Z_PY_VI$cl
VI(z_0,t(memb_Z_PY))
# Horizontal bound of the credible ball
credibleball(memb_Z_PY_VI$cl,t(Z_PY[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
missclass(memb_Z_PY,Y,a=1,b=1)


c_Z_PY <- pr_cc(Z_PY_x[,(burn_in+1):N_iter])
# Point estimate
memb_Z_PY_VI <- minVI(c_Z_PY,method="avg",max.k=20)
memb_Z_PY <- memb_Z_PY_VI$cl
VI(z_0,t(memb_Z_PY))
# Horizontal bound of the credible ball
credibleball(memb_Z_PY_VI$cl,t(Z_PY_x[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
misclass(memb_Z_PY,Y,a=1,b=1)


# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])
# Point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl
VI(z_0,t(memb_Z_GN))
# Horizontal bound of the credible ball
credibleball(memb_Z_GN_VI$cl,t(Z_GN[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
misclass(memb_Z_GN,Y,a=1,b=1)


c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
# Point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl
VI(z_0,t(memb_Z_GN))
# Horizontal bound of the credible ball
credibleball(memb_Z_GN_VI$cl,t(Z_GN_x[,(burn_in+1):N_iter]))[[5]]
# Misclassification error
misclass(memb_Z_GN,Y,a=1,b=1)
```

Also these **results are in line with our previous discussion**. Finally, the code to **reproduce the upper panel of Figure 2** and save it, is provided below.

``` r
c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

c_Z_GN_x <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
memb_Z_GN_VI_x <- minVI(c_Z_GN_x,method="avg",max.k=20)
memb_Z_GN_x <- memb_Z_GN_VI_x$cl

row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN,V,1)))
names(row_plot_GN) <- "memb_Z_GN"
rownames(c_Z_GN) <- rownames(row_plot_GN)
mycolors <- (brewer.pal(length(unique(row_plot_GN$memb_Z_GN)),"Set3"))
names(mycolors) <- unique(row_plot_GN$memb_Z_GN)
mycolors <- list(memb_Z_GN = mycolors)

Marg <- pheatmap(c_Z_GN,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN_x,V,1)))
names(row_plot_GN) <- "memb_Z_GN_x"
rownames(c_Z_GN_x) <- rownames(row_plot_GN)
mycolors <- (brewer.pal(length(unique(row_plot_GN$memb_Z_GN_x)),"Set3"))
names(mycolors) <- unique(row_plot_GN$memb_Z_GN_x)
mycolors <- list(memb_Z_GN_x = mycolors)

Cov <- pheatmap(c_Z_GN_x,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

diag(Y) <- 1
row_plot_Y <- as.data.frame(as.factor(matrix(z_0,V,1)))
names(row_plot_Y) <-"z_0"
rownames(Y) <- rownames(row_plot_Y)
mycolors <- (brewer.pal(length(unique(row_plot_Y$z_0)),"Set3"))
names(mycolors) <- unique(row_plot_Y$z_0)
mycolors <- list(z_0 = mycolors)

Network <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

g <- grid.arrange(Network[[4]],Marg[[4]],Cov[[4]],nrow=1,vp=viewport(width=1, height=1))
g2 <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

png("scenario1.png",width=1500,height=501)
print(g2)
dev.off()
```
![](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/scenario1.png)

Refer to the simulation section in the article [**Extended stochastic block models**]()  for detailed comments on the above Figures.
