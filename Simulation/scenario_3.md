Simulation: Scenario 3
================
This tutorial contains guidelines and code to perform the analyses for the third scenario `network_3.RData` considered in the simulation study of the article **Extended stochastic block models with application to criminal networks**. In particular, you will find a detailed step-by-step guide and `R` code to **implement the collapsed Gibbs sampler presented in the article** and to **fully reproduce the results for the third scenario** presented in Section 4 of the article. For implementation purposes, please **execute the code below considering the same order in which is presented**.

Import the data
================
To start the analysis, **set the working directory** where the `README` file and the subfolders `Source`, `Simulation` and `Application` are located. Once this has been done, **clean the workspace, and load the data along with useful** `R` **packages**.

``` r
rm(list=ls())
source("Source/esbm.R")
Rcpp::sourceCpp('Source/stirling.cpp')

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
library(randnet)

load("Simulation/network_3.RData")
V <- dim(Y)[1]

# note that Y must have diagonal equal to 0
diag(Y)
```

As discussed in the article, the network under analysis has *V=80* nodes and *5* groups displaying **community and core-periphery structures** [see third row of Figure 3 in the article]. 

Setting the hyperparameters
================
For the hyperparameters of the `Beta(a,b)` **priors on the block probabilities** we follow common implementations of stochastic block models and consider the default values `a=1` and `b=1` to induce a **uniform** prior. Less straightforward is instead the choice of the hyperparameters for the **Gibbs-type priors on the random partition**. A possibility to address this issue is to specify such quantities in order to obtain a value for the expectation of the non-empty number of groups `H` that matches some prior knowledge. Below, we provide the **code to obtain such a quantity as a function of pre-specified hyperparameters for the four relevant examples of Gibbs-type priors** discussed in the article. 

``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
sigma_dm <- 0   
H_dm <- 50 # Conservative upper bound 
beta_dm <- 3.5/H_dm 
round(expected_cl_py(V, sigma = sigma_dm, theta = beta_dm*H_dm, H = H_dm))

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
sigma_dp <- 0   
H_dp <- Inf 
alpha_dp <- 3
round(expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp))

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
sigma_py <- 0.6
H_py <- Inf 
alpha_py <- -0.3
round(expected_cl_py(V, sigma = sigma_py, theta = alpha_py, H = H_py))

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------
gamma <- 0.45
probs_gnedin <- HGnedin(V, 1:V, gamma = gamma)
round(sum(1:V*probs_gnedin))
```

Here we **set the hyperparameters so that the expectation of `H` is close to *10* under all the four priors**. This is twice as many as the true number of groups to check whether our results are robust to hyperparameters settings.

Posterior computation via collapsed Gibbs sampler
================
This section contains the code to **implement the collapsed Gibbs sampler for ESBM** [function `esbm()`] and to **evaluate marginal likelihoods** [function `log_pY_z()`] for model selection. Such a code is applied to [i] first select among the four relevant examples of unsupervised Gibbs-type priors discussed in the article, and then [ii] check whether introducing informative node attributes further improves the performance of the selected best unsupervised prior. See the source code `esbm.R` for a detailed description of the inputs and the outputs of the two functions `esbm()` and `log_pY_z()`.

Implementation without node-specific attributes
------------------
Let us first **perform posterior computation for the model without node attributes**. To do this, execute the code below.

``` r
N_iter <- 50000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

my_prior <- "DM"
Z_DM <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 3.5/50, H_DM = 50)

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

my_prior <- "DP"
Z_DP <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = 3, sigma_PY = 0)

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

my_prior <- "PY"
Z_PY <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = -0.3, sigma_PY = 0.6)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

my_prior  <- "GN"
Z_GN <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.45)
```

Once the above steps have been done, **compute the logarithm of the marginal likelihoods** that will be used for comparing the performance of the different prior specifications, and **save the output** in the file `Posterior_No_Attributes3.RData`.

``` r
# compute the logarithm of the marginal likelihoods under the different priors

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

# save the output
save(Z_DP,l_y_DP,Z_PY,l_y_PY,Z_GN,l_y_GN,Z_DM,l_y_DM,file="Simulation/Posterior_No_Attributes3.RData")
rm(Z_DP,l_y_DP,Z_PY,l_y_PY,Z_GN,l_y_GN,Z_DM,l_y_DM)
```

Implementation with node-specific attributes
------------------
As shown in Table 2 in the article, the **GN process** yields the best performance among the four relevant examples of unsupervised Gibbs-type priors discussed in the article. Hence, let us now **perform posterior computation for the supervised GN prior with node-specific attributes** coinciding, in this case, with the true membership labels. To accomplish this goal, execute the code below.

``` r
N_iter  <- 50000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# define the vector with node attributes
my_x <- c(rep(1,35),rep(2,20),rep(3,15),rep(4,5),rep(5,5))

# set hyperparameters for the Dirichlet-Multinomial cohesion function (see the article)
my_alpha_xi <- rep(1,5)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

my_prior <- "GN"
Z_GN_x <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.45, x = my_x, alpha_xi = my_alpha_xi)
```

Also in this case we **compute the logarithm of the marginal likelihood** that will be used for assessing performance, and **save the output** in the file `Posterior_Attributes3.RData`.

``` r
# compute the logarithm of the marginal likelihood under supervised GN prior

l_y_GN_x <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_GN_x[t] <- log_pY_z(Y,Z_GN_x[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

# save the output
save(Z_GN_x,l_y_GN_x,file="Simulation/Posterior_Attributes3.RData")
rm(Z_GN_x,l_y_GN_x)
```

Posterior inference under ESBM [Table 2: Scenario 3]
================
This section contains the **code to perform estimation, uncertainty quantification and model selection for ESBM** leveraging the samples from the collapsed Gibbs sampler. In particular, we **reproduce the analyses in Table 2 of the article**, for **scenario 3**. To accomplish this goal let us first **import the MCMC samples**, and define the `burn_in` along with the vector `z_0` containing the true group labels. 

``` r
burn_in <- 10000
z_0 <- c(rep(1,35),rep(2,20),rep(3,15),rep(4,5),rep(5,5))
load("Simulation/Posterior_No_Attributes3.RData")
load("Simulation/Posterior_Attributes3.RData")
```
Before performing posterior inference, let us **visualize the traceplots for the logarithm of the likelihood in Eq. [1]**, evaluated at the MCMC samples of `z` under the different priors, both with and without nodal attributes.

``` r
traceplot <- melt(cbind(l_y_DM,l_y_DP,l_y_PY,l_y_GN,l_y_GN_x))
traceplot <- traceplot[,-2]

traceplot$Group <- c(rep("DM [unsup]",N_iter),rep("DP [unsup]",N_iter),rep("PY [unsup]",N_iter),rep("GN [unsup]",N_iter),rep("GN [sup]",N_iter))
traceplot$Group <- factor(traceplot$Group,levels=c("DM [unsup]","DP [unsup]","PY [unsup]","GN [unsup]","GN [sup]"))

Trace <- ggplot(traceplot,aes(y=value,x=X1)) + geom_line() + facet_grid(.~Group) + theme_bw() + labs(y="",x="")
Trace
```

The above traceplots confirm that our Gibbs sampler has **satisfactory mixing and rapid convergence**. Due to the stability of the chains for the quantity in Eq. [1], we can reliably compute the **logarithm of the marginal likelihoods** for the different priors and models via the harmonic mean in Eq. [14] (*see the third column in Table 2*).

``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

l_y_DM <- l_y_DM[(burn_in+1):N_iter]
neg_l_y_DM <- -c(l_y_DM)
l_y_post_DM <- log(length(l_y_DM))-max(neg_l_y_DM)-log(sum(exp(neg_l_y_DM-max(neg_l_y_DM))))
l_y_post_DM

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

l_y_DP <- l_y_DP[(burn_in+1):N_iter]
neg_l_y_DP <- -c(l_y_DP)
l_y_post_DP <- log(length(l_y_DP))-max(neg_l_y_DP)-log(sum(exp(neg_l_y_DP-max(neg_l_y_DP))))
l_y_post_DP

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

l_y_PY <- l_y_PY[(burn_in+1):N_iter]
neg_l_y_PY <- -c(l_y_PY)
l_y_post_PY <- log(length(l_y_PY))-max(neg_l_y_PY)-log(sum(exp(neg_l_y_PY-max(neg_l_y_PY))))
l_y_post_PY

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

As it can be noticed, the **Gnedin process performs slightly better** relative to the other priors. Moreover, as expected, the overall **learning process benefits from informative node-specific attributes**.  

The **posterior mean of the variation of information (VI) distance from the true partition** `z_0` can be instead obtained using the `VI()` function within the `mcclust.ext` package [Wade and Ghahramani, 2018] as follows (*see the sixth column in Table 2*).


``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
VI(z_0,t(Z_DM[,(burn_in+1):N_iter]))

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
VI(z_0,t(Z_DP[,(burn_in+1):N_iter]))

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
VI(z_0,t(Z_PY[,(burn_in+1):N_iter]))

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------
VI(z_0,t(Z_GN[,(burn_in+1):N_iter]))
VI(z_0,t(Z_GN_x[,(burn_in+1):N_iter]))
```

The above results **confirm the rankings** obtained from the analysis of the marginal likelihoods.

As discussed in the article, accurate learning of the underlying number of groups is a fundamental goal. Hence, let us study the **quartiles of the posterior distribution for the number of non-empty groups** under the different priors and models (*see the ninth column in Table 2*). 

``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
quantile(apply(Z_DM[,(burn_in+1):N_iter],2,max))[c(2:4)]

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
quantile(apply(Z_DP[,(burn_in+1):N_iter],2,max))[c(2:4)]

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
quantile(apply(Z_PY[,(burn_in+1):N_iter],2,max))[c(2:4)]

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------
quantile(apply(Z_GN[,(burn_in+1):N_iter],2,max))[c(2:4)]
quantile(apply(Z_GN_x[,(burn_in+1):N_iter],2,max))[c(2:4)]
```

Also for this measure, the inclusion of **informative node-specific attributes provides improved performance**. 

To complete Table 2, let us obtain **point estimates** and **credible balls** for the group assignments of the different nodes. This is done by adapting the methods presented in Wade and Ghahramani (2018) and implemented in the `R` package `mcclust.ext`. To apply these strategies we also require an estimate of the **co-clustering matrix**, whose generic element `c[v,u]` encodes the relative frequency of MCMC samples in which nodes `v` and `u` are in the same cluster. Such an estimate can be obtained via the function `pr_cc()` in the source code `esbm.R`. We also study the **misclassification error** using the function `misclass()` in the source code `esbm.R` (*see the twelfth column in Table 2 for the `VI` distance between the estimated partition and the 95% credible bound*). 


``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

c_Z_DM <- pr_cc(Z_DM[,(burn_in+1):N_iter])

# point estimate
memb_Z_DM_VI <- minVI(c_Z_DM,method="avg",max.k=20)
memb_Z_DM <- memb_Z_DM_VI$cl

# horizontal bound of the credible ball
credibleball(memb_Z_DM_VI$cl,t(Z_DM[,(burn_in+1):N_iter]))[[5]]

# misclassification error [in-sample for edges]
misclass(memb_Z_DM,Y,a=1,b=1)

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

c_Z_DP <- pr_cc(Z_DP[,(burn_in+1):N_iter])

# point estimate
memb_Z_DP_VI <- minVI(c_Z_DP,method="avg",max.k=20)
memb_Z_DP <- memb_Z_DP_VI$cl

# horizontal bound of the credible ball
credibleball(memb_Z_DP_VI$cl,t(Z_DP[,(burn_in+1):N_iter]))[[5]]

# misclassification error [in-sample for edges]
misclass(memb_Z_DP,Y,a=1,b=1)

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

c_Z_PY <- pr_cc(Z_PY[,(burn_in+1):N_iter])

# point estimate
memb_Z_PY_VI <- minVI(c_Z_PY,method="avg",max.k=20)
memb_Z_PY <- memb_Z_PY_VI$cl

# horizontal bound of the credible ball
credibleball(memb_Z_PY_VI$cl,t(Z_PY[,(burn_in+1):N_iter]))[[5]]

# misclassification error [in-sample for edges]
misclass(memb_Z_PY,Y,a=1,b=1)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# horizontal bound of the credible ball
credibleball(memb_Z_GN_VI$cl,t(Z_GN[,(burn_in+1):N_iter]))[[5]]

# misclassification error [in-sample for edges]
misclass(memb_Z_GN,Y,a=1,b=1)

# ------------------------------------

c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# horizontal bound of the credible ball
credibleball(memb_Z_GN_VI$cl,t(Z_GN_x[,(burn_in+1):N_iter]))[[5]]

# misclassification error [in-sample for edges]
misclass(memb_Z_GN,Y,a=1,b=1)
```

Also these **results are in line with our previous discussion**. 

Comparison with state-of-the-art competitors [Table 3: Scenario 3]
================
Let us now reproduce the results for **scenario 3** in **Table 3**. Here the focus is on comparing the performance of **ESBM with GN prior** and **state–of–the–art competitors** in the `R` libraries `igraph` and `randnet`. Such alternative strategies include the **Louvain algorithm** [Blondel et al., 2008], **spectral clustering** [Von Luxburg, 2007] and **regularized spectral clustering** [Amini et al., 2013]. 

To compute the errors in the column `ERROR [EST]`, we require the **matrix of true edge probabilities** which have been used to simulate `network_3.RData`.

``` r
pi_true <- matrix(0,V,V)

block_matrix <- matrix(0.25,5,5)
block_matrix[1,1] <- 0.75
block_matrix[2,2] <- 0.75
block_matrix[3,3] <- 0.75
block_matrix[4,4] <- 0.75
block_matrix[5,5] <- 0.75
block_matrix[4,1] <- block_matrix[1,4] <- 0.75
block_matrix[4,5] <- block_matrix[5,4] <- 0.75
block_matrix[2,5] <- block_matrix[5,2] <- 0.75
block_matrix[3,5] <- block_matrix[5,3] <- 0.75

for (v in 2:V){
	for (u in 1:(v-1)){
pi_true[v,u] <- block_matrix[z_0[v],z_0[u]]
	}
}
pi_true <- pi_true+t(pi_true)
```

The following code provides the performance measures in *columns 3, 6 and 9 of Table 3* for the **unsupervised GN prior**.

``` r
# ------------------------------------
# GNEDIN PROCESS UNSUPERVISED
# ------------------------------------

c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# estimated H
length(table(memb_Z_GN))

# VI distance between estimated and true partition
VI(z_0,t(memb_Z_GN))

# estimation error for edge probabilities 
mean(abs(lowerTriangle(edge_est(memb_Z_GN,Y,a=1,b=1))-lowerTriangle(pi_true)))
```

Similarly, the performance measures in *columns 3, 6 and 9 of Table 3* for the **supervised GN prior** are obtained as follows.

``` r
# ------------------------------------
# GNEDIN PROCESS SUPERVISED
# ------------------------------------

c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# estimated H
length(table(memb_Z_GN))

# VI distance between estimated and true partition
VI(z_0,t(memb_Z_GN))

# estimation error for edge probabilities 
mean(abs(lowerTriangle(edge_est(memb_Z_GN,Y,a=1,b=1))-lowerTriangle(pi_true)))
```

Let us now focus on the performance of the **Louvain algorithm** in *columns 3, 6 and 9 of Table 3*. To implement this strategy we rely on the function `cluster_louvain()` within the `R` library `igraph`.

``` r
# ------------------------------------
# LOUVAIN ALGORITHM
# ------------------------------------

# transform the adjacency matrix into an igraph object
net <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# point estimate
Louv <- cluster_louvain(net)$membership

# estimated H
length(table(Louv))

# VI distance between estimated and true partition
VI(z_0,t(Louv))

# estimation error for edge probabilities 
mean(abs(lowerTriangle(edge_est(Louv,Y,a=1,b=1))-lowerTriangle(pi_true)))
```

In implementing **spectral clustering**, we first need to **specify the number of groups** `sel_H`. To do this, we consider a variety of model selection criteria available in the `R` library `randnet`, and set the number of groups equal to the median of the values of `H` estimated under the different strategies. 

``` r
set.seed(1)
H_select <- rep(0,8)

# Le and Levina (2015)
bhmc <- BHMC.estimate(Y,K.max=10)
H_select[1] <- bhmc$K

# Wang and Bickel (2017)
lrbic <- LRBIC(Y,Kmax=10)
H_select[2] <- lrbic$SBM.K

# Chen and Lei (2018)
ncv <- NCV.select(Y,max.K=10)
H_select[3] <- which.min(ncv$l2)
H_select[4] <- which.min(ncv$dev)

# Li et al. (2020)
ecv <- ECV.block(Y,max.K=10)
H_select[5] <- which.min(ecv$l2)
H_select[6] <- which.min(ecv$dev)

# Li et al. (2020)
ecv.R <- ECV.Rank(Y,10,weighted=FALSE,mode="undirected")
H_select[7] <- ecv.R$sse.rank
H_select[8] <- ecv.R$auc.rank

sel_H <- round(median(H_select))
```

Once `sel_H` is available, we can obtain the performance measures in *columns 3, 6 and 9 of Table 3* under **spectral clustering** and **regularized spectral clustering** as follows.

``` r
# ------------------------------------
# SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=0)$cluster

# estimated H
length(table(sc))

# VI distance between estimated and true partition
VI(z_0,t(sc))

# estimation error for edge probabilities 
mean(abs(lowerTriangle(edge_est(sc,Y,a=1,b=1))-lowerTriangle(pi_true)))


# ------------------------------------
# REGULARIZED SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
r_sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=1)$cluster

# estimated H
length(table(r_sc))

# VI distance between estimated and true partition
VI(z_0,t(r_sc))

# estimation error for edge probabilities 
mean(abs(lowerTriangle(edge_est(r_sc,Y,a=1,b=1))-lowerTriangle(pi_true)))
```

Predictive performance for the group membership of new nodes
================
We study the performance of the **supervised GN process prior** (which yields the most accurate inference within the **ESBM** class for **scenario 3**) in **predicting the group membership of new incoming nodes**. To accomplish this goal, we first simulate the edges between `300` new nodes and those included in the original network `Y`. Among these incoming nodes, `50` belong to a new group not yet observed in `Y` (which is characterized by low connection probability with nodes in the original `5` clusters).

``` r
set.seed(1)
V_new <- 300

# membership indicators for the incoming nodes
memb_new <- c(rep(1,50),rep(2,50),rep(3,50),rep(4,50),rep(5,50),rep(6,50))

# create empty matrix of edges between the 300 new nodes and those in the original network
Y_new <- matrix(0,300,V)

# create true block probability matrix (including the new group not yet observed in Y)
block_matrix_pred <- matrix(0,6,6)
block_matrix_pred[1:5,1:5] <- block_matrix
block_matrix_pred[6,] <- block_matrix_pred[,6] <- 0.05

# simulate the new edges
for (v in 1:V_new){
	for (u in 1:V){
Y_new[v,u]<-rbinom(1,1,prob=block_matrix_pred[memb_new[v],z_0[u]])
	}
}

# create an augmented (V+1)x(V+1) adjacency matrix which adds to Y a last row and column that will
# contain the edges between a new incoming node and the V=80 original ones. 
# Note: prediction is done for one node at-a-time. 

Y_augmented <- matrix(0,V+1,V+1)
Y_augmented[1:V,1:V] <- Y
```

Once the edges for the new incoming nodes have been simulated, we can compute a **plug–in estimate for the predictive probabilities of the cluster allocations** by applying Eq. [15] to each new node (one at-a-time). This can be done via the function `pred_esbm()` in the source code `esbm.R`. The final prediction of the group allocation for each node is that label with the highest predicted probability. Once these predictions have been obtained, we also compute the **misclassification error** by comparing such quantities with the true labels in `memb_new`.


``` r
# create matrix which will be populated with the predictive probabilities 
Post_Prob <- matrix(0,V_new,6)

# compute predictive probabilities for the 300 nodes (one at-a-time).
for (v in 1:V_new){
Y_augmented[V+1,1:V] <- Y_augmented[1:V,V+1] <- Y_new[v,]
Post_Prob[v,] <- pred_esbm(Y_augmented, prior="GN", z_hat=memb_Z_GN,a = 1, b = 1, gamma_GN = 0.45)
}

# compute misclassification error
(V_new-sum(diag(table(memb_new,apply(Post_Prob,1,which.max)))))/V_new
```


Graphical representation [Figure 3: Row 3]
================

Finally, the code to **reproduce the third row of Figure 3** is provided below.

``` r
diag(Y) <- 0
row_plot_Y <- as.data.frame(as.factor(matrix(z_0,V,1)))
names(row_plot_Y) <-"z_0"
rownames(Y) <- rownames(row_plot_Y)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y$z_0)
mycolors <- list(z_0 = mycolors)

Network <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

# ------------------------------------

c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN,V,1)))
names(row_plot_GN) <- "memb_Z_GN"
rownames(c_Z_GN) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[4],brewer.pal(9,"YlOrBr")[3],brewer.pal(10,"RdBu")[7],brewer.pal(10,"PRGn")[7])
names(mycolors) <- unique(row_plot_GN$memb_Z_GN)
mycolors <- list(memb_Z_GN = mycolors)

Marg <- pheatmap(c_Z_GN,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

# ------------------------------------

c_Z_GN_x <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
memb_Z_GN_VI_x <- minVI(c_Z_GN_x,method="avg",max.k=20)
memb_Z_GN_x <- memb_Z_GN_VI_x$cl

row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN_x,V,1)))
names(row_plot_GN) <- "memb_Z_GN_x"
rownames(c_Z_GN_x) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_GN$memb_Z_GN_x)
mycolors <- list(memb_Z_GN_x = mycolors)

Cov <- pheatmap(c_Z_GN_x,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

# ------------------------------------

g <- grid.arrange(Network[[4]],Marg[[4]],Cov[[4]],nrow=1,vp=viewport(width=1, height=1))
g2 <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(g2)
```

Refer to Section 4 in the article for detailed comments on the above Figures.
