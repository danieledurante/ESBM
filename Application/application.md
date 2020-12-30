Application: The Infinito Network
================
This tutorial contains guidelines and code to perform the analyses in the application to the *Infinito network* [`crime_net.RData`] illustrated in the article **Extended stochastic block models with application to criminal networks**. In particular, you will find a detailed step-by-step guide and `R` code to **pre-process the original data**, **implement the collapsed Gibbs sampler developed in the article** and **fully reproduce the results** discussed in Sections 1 and 5 of the article. For implementation purposes, please **execute the code below considering the same order in which is presented**.

Data pre-processing
================
The original data [`NDRANGHETAMAFIA_2M.csv`] are available in the zip directory `Ndrangheta Mafia 2 CSV.zip` at the link https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndranghetamafia2, and comprise information on the co–participation of 156 suspects at 47 monitored summits of the criminal organization *La Lombardia*, as reported in the **judicial acts** which can be accessed at, e.g., https://liberavco.liberapiemonte.it/wp-content/uploads/sites/13/2012/04/Operazione-Infinito-Ordinanza-di-Custodia-Cautelare.pdf.  

Please download the original [`NDRANGHETAMAFIA_2M.csv`](https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndranghetamafia2/Ndrangheta%20Mafia%202%20CSV.zip?attredirects=0&d=1) file and save it in the `Application` folder. To load these data in `R`, first **set the working directory** where the `README` file and the subfolders `Source`, `Simulation` and `Application` are located. Once this has been done, **clean the workspace, and load the data in** `R`.
 
``` r
rm(list=ls())

A <- read.csv(file="Application/NDRANGHETAMAFIA_2M.csv",header=TRUE, stringsAsFactors = TRUE)
```

A first careful double check between the information in the judicial acts and the data in `A` shows that **in few cases the attendance of a suspect to a summit was not reported in** `A`. Hence, let us impute this information. 

``` r
A[23,20] <- 1
A[48,23] <- 1
A[78,35] <- 1
A[115,22] <- 1
```

As discussed in Section 1.1, our overarching goal in this application is to **shed light on the internal structure of** *La Lombardia* via inference on the block connectivity patterns among its affiliates. Hence, as a first pre-processing step, let us **exclude from the analysis those individuals who never attended a summit as well as those who have not been recognized during the investigation process**. 

``` r
# suspects who never attended a summit
sel_empty <- which(apply(as.matrix(A[,-1]),1,sum)==0)

# suspects who have not been recognized during the investigation process
A[c(38,105,106,125,135),1]

# indicators of the two groups of suspects to be excluded
sel_empty <- c(c(sel_empty),c(38,105,106,125,135))

# remove these suspects from the dataset
A <- A[-sel_empty,]

# create a vector with the actors’ names
actors <- A[,1]
actors <- droplevels(actors)
```

Let us now create the **binary adjacency matrix** indicating presence or absence of a co–attendance in at least one of the monitored summits; see Section 1.1 for additional details and motivations regarding this dichotomization process.

``` r
A <- as.matrix(A[,-1])
A <- A%*%t(A)
A <- (A>0)*1
diag(A) <- 0
rownames(A) <- colnames(A) <- c(1:dim(A)[1])
```

Recalling Section 1.1, information on presumed **locale membership** and **leadership role** can be retrieved, for each suspect of interest, from the **judicial acts of Operazione Infinito**. These node attributes are imputed manually by controlling the information reported in the judicial acts for each suspect in the vector `actors`.  Most of these information are contained in pages 21-32 of the acts, but for some suspects further knowledge is available in the subsequent pages. 

``` r
# ----------------------------------------------------
# Locale membership ("OUT": Suspects not belonging to La Lombardia. "MISS": Information not available)

Locale <- c("C","OUT","A","MISS","O","A","MISS","D","D","D","D","D","C","P","L","L","Q","MISS","B","OUT","B","B","I","MISS","OUT","D","A","O","N","N","H","OUT","D","E","G","G","L","A","OUT","Q","C","OUT","Q","L","C","MISS","C","C","F","C","OUT","D","A","B","B","E","M","MISS","C","C","C","B","H","C","C","E","E","E","E","C","MISS","L","A","A","E","E","C","E","E","E","C","MISS","OUT","C","C","E","G","A","A","B","I","I","A","B","B","OUT","I","A","G","N","E","D","F","OUT","OUT","C","D","C","MISS","MISS","C","MISS","E","E","C","MISS","OUT","B","L","A","D","D","O","MISS","B","D","O","D","D","A","A","I","C","MISS","MISS","MISS","A","A","F","E","C","Q","H","B","B","B")   

# ----------------------------------------------------
# Leadership role ("miss”: Information not available)

Role <- c("aff","aff","aff","miss","aff","boss","miss","boss","boss","aff","aff","aff","aff","aff","aff","aff","aff","miss","aff","boss","boss","boss","boss","miss","boss","boss","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","aff","aff","boss","aff","aff","aff","boss","aff","aff","miss","aff","aff","aff","aff","boss","aff","aff","aff","aff","aff","boss","miss","aff","aff","aff","aff","boss","boss","boss","aff","boss","aff","aff","boss","miss","aff","aff","boss","boss","aff","aff","aff","aff","aff","aff","miss","boss","aff","aff","aff","aff","aff","aff","boss","aff","boss","aff","aff","aff","aff","boss","boss","boss","boss","aff","aff","aff","aff","aff","aff","aff","boss","miss","miss","aff","miss","aff","aff","aff","miss","aff","aff","boss","aff","aff","aff","aff","miss","aff","aff","boss","aff","boss","aff","aff","aff","aff","miss","miss","miss","aff","aff","boss","aff","aff","aff","boss","aff","aff","boss")
```

As mentioned above, the overarching focus is on the criminal organization *La Lombardia*. Hence, **let us restrict our attention only to those suspects who are known to belong to one of the *locali* of such a criminal organization**.

``` r
# indicators of those suspects who are not known to be part of the organization La Lombardia 
sel_miss <- which(Locale=="MISS" | Locale=="OUT")

# remove such suspects from the dataset
Locale_temp <- Locale[-sel_miss]
Role_temp <- Role[-sel_miss]
actors_temp <- actors[-sel_miss]
A <- A[-sel_miss,-sel_miss]
rownames(A) <- colnames(A) <- c(1:dim(A)[1])

# interaction between Locale and Role
RoleLocale_temp <- paste(Role_temp,Locale_temp,sep="_")
```

To assess out–of–sample predictive performance, we **perform inference on the** `V = 84` **suspects affiliated to the `5` most populated locali**, and **hold out as a test set the** `34` **members of those smaller–sized locali** with ≤ 6 monitored affiliates; see Section 1.1 for additional motivations regarding the choice of such training and test sets.

``` r
table(Locale_temp)

# Create the dataset used for modeling and inference
sel <- which(Locale_temp=="A" | Locale_temp=="B" | Locale_temp=="C" | Locale_temp=="D" | Locale_temp=="E")
Y <- A[sel,sel]
RoleLocale <- RoleLocale_temp[sel]
Role <- Role_temp[sel]
Locale <- Locale_temp[sel]
actors <- actors_temp[sel]
rownames(Y) <- colnames(Y) <- c(1:dim(Y)[1])

# Create the test dataset used for assessing predictive performance
Y_test <- A[-sel,sel]
rownames(Y_test) <- c(1:dim(Y_test)[1])
colnames(Y_test) <- c(1:dim(Y_test)[2])
RoleLocale_test <- RoleLocale_temp[-sel]
Role_test <- Role_temp[-sel]
Locale_test <- Locale_temp[-sel]
actors_test <- actors_temp[-sel]
```

Finally, **save the data** in the file `crime_net.RData` and clean the workspace.

``` r
save(Y,Y_test,Locale,Locale_test,Role,Role_test,RoleLocale,RoleLocale_test,file="Application/crime_net.RData")

rm(list=ls())
```

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

load("Application/crime_net.RData")
V <- dim(Y)[1]

# note that Y must have diagonal equal to 0
diag(Y)
```

Setting the hyperparameters
================
For the hyperparameters of the `Beta(a,b)` **priors on the block probabilities** we follow common implementations of stochastic block models and consider the default values `a=1` and `b=1` to induce a **uniform** prior. Less straightforward is instead the choice of the hyperparameters for the **Gibbs-type priors on the random partition**. A possibility to address this issue is to specify such quantities in order to obtain a value for the expectation of the non-empty number of groups `H` that matches some prior knowledge. Below, we provide the **code to obtain such a quantity as a function of pre-specified hyperparameters for the four relevant examples of Gibbs-type priors** discussed in the article. 

``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------
sigma_dm <- 0   
H_dm <- 50 # Conservative upper bound 
beta_dm <- 12/H_dm 
round(expected_cl_py(V, sigma = sigma_dm, theta = beta_dm*H_dm, H = H_dm))

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
sigma_dp <- 0   
H_dp <- Inf 
alpha_dp <- 8
round(expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp))

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
sigma_py <- 0.725
H_py <- Inf 
alpha_py <- -0.350
round(expected_cl_py(V, sigma = sigma_py, theta = alpha_py, H = H_py))

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------
gamma <- 0.3
probs_gnedin <- HGnedin(V, 1:V, gamma = gamma)
round(sum(1:V*probs_gnedin))
```

Here we **set the hyperparameters so that the expectation of `H` is close to *20* under all the four priors**. This value is four times the number of *locali* in the network, which seems reasonably conservative.

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
Z_DM <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, beta_DM = 12/50, H_DM = 50)

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

my_prior <- "DP"
Z_DP <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

my_prior <- "PY"
Z_PY <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = -0.350, sigma_PY = 0.725)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

my_prior  <- "GN"
Z_GN <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.3)
```

Once the above steps have been done, **compute the logarithm of the marginal likelihoods** that will be used for comparing the performance of the different prior specifications, and **save the output** in the file `Posterior_No_Attributes.RData`.

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

# Save the output
save(Z_DP,l_y_DP,Z_PY,l_y_PY,Z_GN,l_y_GN,Z_DM,l_y_DM,file="Application/Posterior_No_Attributes.RData")
rm(Z_DP,l_y_DP,Z_PY,l_y_PY,Z_GN,l_y_GN,Z_DM,l_y_DM)
```

Implementation with node-specific attributes
------------------
As shown in Table 4 in the article, the **GN process** yields the best performance among the four relevant examples of unsupervised Gibbs-type priors discussed in the article. Hence, let us now **perform posterior computation for the supervised GN prior with node-specific attributes**. As clarified in Section 5, we consider as node attribute an external variable in which:

- the class of each **affiliate** corresponds to the associated locale
- the **bosses** share a common label indicating that such members have a leadership role within the organization
- a **subset of the affiliates of locale `D`** who are know from the judicail acts to cover a peripheral role within the locale are assigned a distinct class

``` r
N_iter  <- 50000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

# define the vector with node attributes
my_x <- c(as.factor(RoleLocale))
my_x[which(my_x>5)] <- 7

# actors with known peripheral roles in locale D
sel_periphery_D <- c(6,8,69,73)
my_x[sel_periphery_D]<- 6

# set hyperparameters for the Dirichlet-Multinomial cohesion function (see the article)
my_alpha_xi <- rep(1,7)

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

my_prior <- "GN"
Z_GN_x <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.3, x = my_x, alpha_xi = my_alpha_xi)
```

Also in this case we **compute the logarithm of the marginal likelihood** that will be used for assessing performance, and **save the output** in the file `Posterior_Attributes.RData`.

``` r
# compute the logarithm of the marginal likelihood under supervised GN prior

l_y_GN_x <- rep(0,N_iter)
for (t in 1:N_iter){
  l_y_GN_x[t] <- log_pY_z(Y,Z_GN_x[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

# save the output
save(Z_GN_x,l_y_GN_x,file="Application/Posterior_Attributes.RData")
rm(Z_GN_x,l_y_GN_x)
```

Posterior inference under ESBM [Table 4]
================
This section contains the **code to perform estimation, uncertainty quantification and model selection for ESBM** leveraging the samples from the collapsed Gibbs sampler. In particular, we **reproduce the analyses in Table 4 of the article**. To accomplish this goal let us first **import the MCMC samples** and define the `burn_in`. 

``` r
burn_in <- 10000
load("Application/Posterior_No_Attributes.RData")
load("Application/Posterior_Attributes.RData")
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

The above traceplots confirm that our Gibbs sampler has **satisfactory mixing and convergence**. Due to the stability of the chains for the quantity in Eq. [1], we can reliably compute the **logarithm of the marginal likelihoods** for the different priors and models via the harmonic mean in Eq. [14] (*see the first column in Table 4*).

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

As it can be noticed, the **Gnedin process performs slightly better** relative to the other priors. Moreover, the **external node attribute**, defined by a combination of *locale* membership and leadership *role*, **yields further improvements in the learning process**. For instance, we obtain a positive-strong evidence in favor of the supervised GN process relative to the unsupervised representation, when studying the **Bayes factor**.

``` r
2*(l_y_post_GN_x-l_y_post_GN)
```

As discussed in the article, accurate learning of the underlying number of groups is a fundamental goal. Hence, let us study the **quartiles of the posterior distribution for the number of non-empty groups** under the different priors and models (*see the third column in Table 4*). 

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

The above results seem to provide consistent evidence for the presence of either **14, 15 or 16 groups in the Infinito network**. 

To complete Table 4, let us obtain **point estimate** and **credible balls** for the group assignments of the different nodes. This is done by adapting the methods presented in Wade and Ghahramani (2018) and implemented in the `R` package `mcclust.ext`. To apply these strategies we also require an estimate of the **co-clustering matrix**, whose generic element `c[v,u]` encodes the relative frequency of MCMC samples in which nodes `v` and `u` are in the same cluster. Such an estimate can be obtained via the function `pr_cc()` in the source code `esbm.R` (*see the second column of Table 4 for the logarithm of the likelihood in Eq. [1] evaluated at the estimated partition. The `VI` distance between the estimated partition and the 95% credible bound is reported in the fourth column of Table 4*). 


``` r
# ------------------------------------
# DIRICHLET MULTINOMIAL
# ------------------------------------

c_Z_DM <- pr_cc(Z_DM[,(burn_in+1):N_iter])

# point estimate
memb_Z_DM_VI <- minVI(c_Z_DM,method="avg",max.k=20)
memb_Z_DM <- memb_Z_DM_VI$cl

# logarithm of likelihood in Eq. [1] evaluated at the estimated partition
log_pY_z(Y,memb_Z_DM,1,1)

# horizontal bound of the credible ball
credibleball(memb_Z_DM_VI$cl,t(Z_DM[,(burn_in+1):N_iter]))[[5]]

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

c_Z_DP <- pr_cc(Z_DP[,(burn_in+1):N_iter])

# point estimate
memb_Z_DP_VI <- minVI(c_Z_DP,method="avg",max.k=20)
memb_Z_DP <- memb_Z_DP_VI$cl

# logarithm of likelihood in Eq. [1] evaluated at the estimated partition
log_pY_z(Y,memb_Z_DP,1,1)

# horizontal bound of the credible ball
credibleball(memb_Z_DP_VI$cl,t(Z_DP[,(burn_in+1):N_iter]))[[5]]

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------

c_Z_PY <- pr_cc(Z_PY[,(burn_in+1):N_iter])

# point estimate
memb_Z_PY_VI <- minVI(c_Z_PY,method="avg",max.k=20)
memb_Z_PY <- memb_Z_PY_VI$cl

# logarithm of likelihood in Eq. [1] evaluated at the estimated partition
log_pY_z(Y,memb_Z_PY,1,1)

# horizontal bound of the credible ball
credibleball(memb_Z_PY_VI$cl,t(Z_PY[,(burn_in+1):N_iter]))[[5]]

# ------------------------------------
# GNEDIN PROCESS
# ------------------------------------

c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# logarithm of likelihood in Eq. [1] evaluated at the estimated partition
log_pY_z(Y,memb_Z_GN,1,1)

# horizontal bound of the credible ball
credibleball(memb_Z_GN_VI$cl,t(Z_GN[,(burn_in+1):N_iter]))[[5]]

# ------------------------------------

c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# logarithm of likelihood in Eq. [1] evaluated at the estimated partition
log_pY_z(Y,memb_Z_GN,1,1)

# horizontal bound of the credible ball
credibleball(memb_Z_GN_VI$cl,t(Z_GN_x[,(burn_in+1):N_iter]))[[5]]
```

Also these **results are in line with our previous discussion**. 

Comparison with state-of-the-art competitors
================
Let us now compare the **deviances** of **ESBM (with GN prior)** and **state–of–the–art competitors** in the `R` libraries `igraph` and `randnet`. Such alternative strategies include the **Louvain algorithm** [Blondel et al., 2008], **spectral clustering** [Von Luxburg, 2007] and **regularized spectral clustering** [Amini et al., 2013]. 

The deviances for the **unsupervised and supervised GN prior** can be easily computed as follows.

``` r
# ------------------------------------
# GNEDIN PROCESS UNSUPERVISED
# ------------------------------------
c_Z_GN <- pr_cc(Z_GN[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# deviance (D)
-log_pY_z(Y,memb_Z_GN,1,1)

# ------------------------------------
# GNEDIN PROCESS SUPERVISED
# ------------------------------------
c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])

# point estimate
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# deviance (D)
-log_pY_z(Y,memb_Z_GN,1,1)
```

Let us now focus on computing the deviance under the **Louvain algorithm**. To implement this strategy we rely on the function `cluster_louvain()` within the `R` library `igraph`.

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

# deviance (D)
-log_pY_z(Y,Louv,1,1)
```

In implementing **spectral clustering**, we first need to **specify the number of groups** `sel_H`. To do this, we consider a variety of model selection criteria available in the `R` library `randnet`, and set the number of groups equal to the median of the values of `H` estimated under the different strategies. 

``` r
set.seed(1)
H_select <- rep(0,8)

# Le and Levina (2015)
bhmc <- BHMC.estimate(Y,K.max=20)
H_select[1] <- bhmc$K

# Wang and Bickel (2017)
lrbic <- LRBIC(Y,Kmax=20)
H_select[2] <- lrbic$SBM.K

# Chen and Lei (2018)
ncv <- NCV.select(Y,max.K=20)
H_select[3] <- which.min(ncv$l2)
H_select[4] <- which.min(ncv$dev)

# Li et al. (2020)
ecv <- ECV.block(Y,max.K=20)
H_select[5] <- which.min(ecv$l2)
H_select[6] <- which.min(ecv$dev)

# Li et al. (2020)
ecv.R <- ECV.Rank(Y,20,weighted=FALSE,mode="undirected")
H_select[7] <- ecv.R$sse.rank
H_select[8] <- ecv.R$auc.rank

sel_H <- round(median(H_select))
```
Note that the outputs in `ncv` and `ecv` provide **empirical support in favor of SBM rather than degree-corrected SBM** in this specific application, thus further motivating the choice of **ESBM** for analyzing the *Infinito network*.

Once `sel_H` is available, we can obtain the deviances under **spectral clustering** and **regularized spectral clustering** as follows.

``` r
# ------------------------------------
# SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=0)$cluster

# estimated H
length(table(sc))

# deviance (D)
-log_pY_z(Y,sc,1,1)

# ------------------------------------
# REGULARIZED SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
r_sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=1)$cluster

# estimated H
length(table(r_sc))

# deviance (D)
-log_pY_z(Y,r_sc,1,1)
```

All the above **deviances are considerably higher relative to those provided by ESBM with GN prior**. Since `sel_H` is lower than the one obtained under the GN prior, let us also compute the deviance for spectral clustering with the same number of clusters `H = 15` inferred under the GN process.

``` r
# ------------------------------------
# SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
sc <- reg.SP(Y,K=15,lap=TRUE,tau=0)$cluster

# deviance (D)
-log_pY_z(Y,sc,1,1)

# ------------------------------------
# REGULARIZED SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
r_sc <- reg.SP(Y,K=15,lap=TRUE,tau=1)$cluster

# deviance (D)
-log_pY_z(Y,r_sc,1,1)
```

Results are still worse relative to those provided by **ESBM**, thereby confirming the **superior performance of the ESBM class also in this specific application**.

Predictive performance for held-out suspects
================
We study the performance of the **supervised GN process prior** (which yields the most accurate inference within the **ESBM** class for the **Infinito network**) in **predicting the group membership of the held-out suspects**. To accomplish this goal, let us first **create an augmented adjacency matrix** which adds to `Y` a last row and column that will contain the edges between a new incoming suspect in the held-out set and the `V=84` original suspects (Note: prediction is done for one node at-a-time). 

``` r
set.seed(1)

Y_augmented <- matrix(0,V+1,V+1)
Y_augmented[1:V,1:V] <- Y
```

Let us now compute a **plug–in estimate for the predictive probabilities of the cluster allocations** by applying Eq. [15] to each held-out suspect (one at-a-time). This can be done via the function `pred_esbm()` in the source code `esbm.R`. The final prediction of the group allocation for each suspect is that label with the highest predicted probability. 


``` r
# create the matrix which will be populated with the predictive probabilities 
Post_Prob <- matrix(0,dim(Y_test)[1],max(memb_Z_GN)+1)

# compute the predictive probabilities for the 34 held-out suspects (one at-a-time)
for (v in 1:dim(Y_test)[1]){
Y_augmented[V+1,1:V] <- Y_augmented[1:V,V+1] <- Y_test[v,]
Post_Prob[v,] <- pred_esbm(Y_augmented, prior="GN", z_hat=memb_Z_GN,a = 1, b = 1, gamma_GN = 0.3)
}

# predicted class membership for held-out suspects
Pred_memb <- apply(Post_Prob,1,which.max)
```

Execute the code below to **compute the predictive measures** discussed in the article.

``` r
# transform the adjacency matrix into an igraph object
net_Y <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# compute the betweenness and local transitivity of each group 
# by averaging those of the suspects within that group
z <- dummy(memb_Z_GN)
transivity_block <- c(c(t(z)%*%(transitivity(net_Y,type="local")))/(table(memb_Z_GN)))
betwenness_block <- c(c(t(z)%*%betweenness(net_Y))/(table(memb_Z_GN)))

# normalize these two measures
norm_transitivity <- (transivity_block-min(transivity_block))/(max(transivity_block)-min(transivity_block))
norm_betwenness <- (betwenness_block-min(betwenness_block))/(max(betwenness_block)-min(betwenness_block))

# assign to each held-out suspect the normalized betweenness 
# and local transitivity of the predicted group
norm_transitivity_pred <- norm_transitivity[Pred_memb]
norm_betwenness_pred <- norm_betwenness[Pred_memb]

# mean of the difference between such quantities for held-out affiliates
mean((norm_transitivity_pred-norm_betwenness_pred)[which(Role_test=="aff")])

# mean of the difference between such quantities for held-out bosses
mean((norm_transitivity_pred-norm_betwenness_pred)[which(Role_test=="boss")])
```


Graphical representations [Figures 1, 2, 4, 5 and 6]
================

Let us first **define the colors for the different categories of the attribute** `RoleLocale`. Such colors will be used for Figures 1, 2, 4, 5 and 6.

``` r
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3],
brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6])
```

The code to **reproduce Figure 1** is provided below.

``` r
# transform the adjacency matrix into an igraph object
net_Y <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# compute the node betweenness to be used for the size of the nodes
betw <- betweenness(net_Y)
# node sizes are proportional to their betweenness 
# Note: for graphical purposes, we consider a monotone transformation of such a measure
V(net_Y)$size <- sqrt(betw/1.5+mean(betw))*2

# node colors indicate the presumed locale membership and role
V(net_Y)$color <- adjustcolor(mycolors[c(as.factor(RoleLocale))], alpha.f = .7)

# node shapes represent leadership role
V(net_Y)$shape <- c("circle","square")[c(as.factor(Role))]

# additional graphical settings
V(net_Y)$frame.color <- "black"
V(net_Y)$label <- "" 
E(net_Y)$color <- brewer.pal(9,"Greys")[3]

# node positions are obtained via force–directed placement
set.seed(12)
l <- layout_with_fr(net_Y)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5)

# plot Figure 1
plot(net_Y, rescale=F, layout=l*1.5,edge.curved=.3,edge.width=0.5)
```

To **display Figure 2** execute the code below.

``` r
# ------------------------------------
# LOUVAIN ALGORITHM
# ------------------------------------

# transform the adjacency matrix into an igraph object
net <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# point estimate
Louv <- cluster_louvain(net)$membership

# to display the block structures, re-order the rows and columns of Y, and the elements 
# in RoleLocale according to the groupings estimated by the Louvain algorithm
sel <- order(Louv)
Louv <- Louv[sel]
Y_Louv <- Y[sel,sel]
RoleLocale_Louv <- RoleLocale[sel]

# plot the adjacency with the grouping structure defined by the Louvain algorithm
row_plotLouv <- as.data.frame(as.factor(matrix(RoleLocale_Louv,V,1)))
names(row_plotLouv) <- "RoleLocale_Louv"
rownames(Y_Louv) <- rownames(row_plotLouv)
names(mycolors) <- sort(unique(row_plotLouv$RoleLocale_Louv))

Adj_Louv <- pheatmap(Y_Louv,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plotLouv,annotation_names_row=F, show_rownames=F,show_colnames=F,legend=F,border_color=FALSE, annotation_legend=F, annotation_colors=list(RoleLocale_Louv = mycolors),gaps_row=c(which(diff(Louv)!=0)),gaps_col=c(which(diff(Louv)!=0)))


# ------------------------------------
# SPECTRAL CLUSTERING
# ------------------------------------

set.seed(1)

# point estimate
sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=0)$cluster

# for graphical purposed, set the order in which groups are displayed so that clusters 
# with nodes having similar attributes are shown close to each other
group_order <- c(1,4,9,3,2,6,7,10,5,8)

# to display the block structures, re-order the rows and columns of Y, and the elements 
# in RoleLocale according to the groupings estimated by spectral clustering 
sel <- which(sc==1)
for (k in 2:length(group_order)){
sel <- c(sel,which(sc==group_order[k]))	
}

sc <- sc[sel]
Y_sc <- Y[sel,sel]
RoleLocale_sc <- RoleLocale[sel]

# plot the adjacency with the grouping structure defined by spectral clustering
row_plotsc <- as.data.frame(as.factor(matrix(RoleLocale_sc,V,1)))
names(row_plotsc) <- "RoleLocale_sc"
rownames(Y_sc) <- rownames(row_plotsc)
names(mycolors) <- sort(unique(row_plotsc$RoleLocale_sc))

Adj_sp <- pheatmap(Y_sc,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plotsc, annotation_names_row=F, show_rownames=F, show_colnames=F, legend=F, border_color=FALSE,annotation_legend=F,annotation_colors=list(RoleLocale_sc = mycolors),gaps_row=c(which(diff(sc)!=0)),gaps_col=c(which(diff(sc)!=0)))

# ------------------------------------
# ESBM WITH SUPERVISED GN PRIOR
# ------------------------------------

set.seed(1)

# point estimate
c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# for graphical purposed, set the order in which groups are displayed so that clusters 
# with nodes having similar attributes are shown close to each other
group_order <- c(15,5,11,2,6,9,10,1,12,7,4,8,3,14,13)

# to display the block structures, re-order the rows and columns of Y, and the elements 
# in RoleLocale according to the groupings estimated under ESBM with supervised GN prior 
sel <- which(memb_Z_GN==15)
for (k in 2:length(group_order)){
sel <- c(sel,which(memb_Z_GN==group_order[k]))	
}

memb_Z_GN <- memb_Z_GN[sel]
Y_esbm <- Y[sel,sel]
RoleLocale_esbm <- RoleLocale[sel]

# plot the adjacency with the grouping structure defined by ESBM with supervised GN prior
row_plotesbm <- as.data.frame(as.factor(matrix(RoleLocale_esbm,V,1)))
names(row_plotesbm) <- "RoleLocale_esbm"
rownames(Y_esbm) <- rownames(row_plotesbm)
names(mycolors) <- sort(unique(row_plotesbm$RoleLocale_esbm))

Adj_esbm <- pheatmap(Y_esbm,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plotesbm, annotation_names_row=F,show_rownames=F, show_colnames=F, legend=F ,border_color=FALSE, annotation_legend=F,annotation_colors=list(RoleLocale_esbm = mycolors),gaps_row=c(which(diff(memb_Z_GN)!=0)),gaps_col=c(which(diff(memb_Z_GN)!=0)))

# ------------------------------------
# COMBINE THE DIFFERENT FIGURES
# ------------------------------------

g <- grid.arrange(Adj_Louv[[4]],Adj_sp[[4]],Adj_esbm[[4]],nrow=1,ncol=3,vp=viewport(width=1, height=1))
Fig_2 <- cowplot::ggdraw(g)+ theme(plot.background =element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(Fig_2)
```

To **reproduce Figure 4**, note that panel (a) has been already obtained in Figure 2 and is available in the object `Adj_esbm`. Hence, it is sufficient to obtain panel (b) via the following code.

``` r
# ------------------------------------
# EDGE PROBABILITY MATRIX [Panel (b)]
# ------------------------------------

set.seed(1)

# point estimate
c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# compute the matrix of estimated edge probabilities under ESBM with supervised GN prior 
# using the function edge_est() in the source code esbm.R
Y_edge_esbm <- edge_est(memb_Z_GN,Y,1,1)

# for graphical purposed, set the order in which groups are displayed so that clusters 
# with nodes having similar attributes are shown close to each other (same as for Adj_esbm)
group_order <- c(15,5,11,2,6,9,10,1,12,7,4,8,3,14,13)

# to display the block structures, re-order the rows and columns of Y, and the elements 
# in RoleLocale according to the groupings estimated under ESBM with supervised GN prior 
sel <- which(memb_Z_GN==15)
for (k in 2:length(group_order)){
sel <- c(sel,which(memb_Z_GN==group_order[k]))	
}

memb_Z_GN <- memb_Z_GN[sel]
Y_edge_esbm <- Y_edge_esbm[sel,sel]
RoleLocale_esbm <- RoleLocale[sel]

# plot the edge probability matrix with the grouping structure defined by ESBM under supervised GN prior
row_plotesbm <- as.data.frame(as.factor(matrix(RoleLocale_esbm,V,1)))
names(row_plotesbm) <- "RoleLocale_esbm"
rownames(Y_edge_esbm) <- rownames(row_plotesbm)
names(mycolors) <- sort(unique(row_plotesbm$RoleLocale_esbm))

Adj_edge_esbm <- pheatmap(Y_edge_esbm,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plotesbm, annotation_names_row=F, show_rownames=F,show_colnames=F, legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=list(RoleLocale_esbm = mycolors),gaps_row=c(which(diff(memb_Z_GN)!=0)),gaps_col=c(which(diff(memb_Z_GN)!=0)))

# ------------------------------------
# COMBINE THE DIFFERENT FIGURES
# ------------------------------------

g <- grid.arrange(Adj_esbm[[4]],Adj_edge_esbm[[4]],nrow=1,ncol=2,vp=viewport(width=1, height=1))
Fig_4 <- cowplot::ggdraw(g)+ theme(plot.background =element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(Fig_4)
```

The code to **reproduce Figure 5** is provided below. 

``` r
set.seed(1)

# point estimate
c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# set hyperparameters of the Beta prior for the block probabilities
a <- 1
b <- 1

# compute the matrix with estimated block probabilities
z <- dummy(memb_Z_GN)
H <- ncol(z)
Abs_Freq <- t(z)%*%Y%*%z
diag(Abs_Freq) <- diag(Abs_Freq)/2
Tot <- t(z)%*%matrix(1,V,V)%*%z
diag(Tot) <- (diag(Tot)-table(memb_Z_GN))/2
Block_freq <- (a+Abs_Freq)/(a+b+Tot)

# define the compositions with respect to locale affiliations and leadership role in each pie-chart
t_standard <- t(table(memb_Z_GN,RoleLocale))

# to provide more direct insights, the composition with respect to role and locale in the lower–sized 
# pie–charts is suitably re–weighted to account for the fact that bosses are less frequent in the 
# network relative to affiliates
t_standard[,which(apply(t_standard,2,sum)<5)] <- t_standard[,which(apply(t_standard,2,sum)<5)]/(apply(table(memb_Z_GN,RoleLocale),2,sum))
t_standard <- t(t_standard)

# define the colors in the pie-charts
values <- list()
for (h in 1:H){values[[h]] <- c(t_standard[h,])}
pie_colors <- list()
pie_colors[[1]] <- adjustcolor(c(mycolors), alpha.f = .7)

# transform the block probability matrix into an igraph object
net_Y <- graph.adjacency(Block_freq, mode=c("undirected"), weighted=TRUE, diag=FALSE)

# node sizes are proportional to cluster cardinality
V(net_Y)$size <- (c(table(memb_Z_GN)))*10

# edge sizes are proportional to the estimated block probabilities
# Note: for graphical purposes, the block probabilities below 0.1 are not displayed
E(net_Y)$width <- (E(net_Y)$weight*(E(net_Y)$weight>0.1))*3

# additional graphical settings
V(net_Y)$label <- NA
V(net_Y)$frame.color <- "black"
E(net_Y)$color <- "grey"

# node positions are obtained via force–directed placement
set.seed(2)
l <- layout_with_fr(net_Y)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5)

# plot Figure 5
plot(net_Y,rescale=F, layout=l*1.8,edge.curved=0.2,vertex.shape="pie", vertex.pie=values,vertex.pie.color=pie_colors,mark.groups=c(3,4,7,8,13,14), mark.col="#f0f0f0", mark.border=NA)
```

Finally, to **reproduce Figure 6** execute the following code.

``` r
# transform the adjacency matrix into an igraph object
net_Y <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# ------------------------------------
# LOUVAIN ALGORITHM
# ------------------------------------
set.seed(1)

# point estimate
Louv <- cluster_louvain(net_Y)$membership

# compute the betweenness and local transitivity of each group 
# by averaging those of the suspects within that group
z <- dummy(Louv)
transivity_block_Louv <- c(c(t(z)%*%(transitivity(net_Y,type="local")))/(table(Louv)))
betwenness_block_Louv <- c(c(t(z)%*%betweenness(net_Y))/(table(Louv)))

# ------------------------------------
# SPECTRAL CLUSTERING
# ------------------------------------
set.seed(1)

# point estimate
sc <- reg.SP(Y,K=sel_H,lap=TRUE,tau=0)$cluster

# compute the betweenness and local transitivity of each group 
# by averaging those of the suspects within that group
z <- dummy(sc)
transivity_block_sc <- c(c(t(z)%*%(transitivity(net_Y,type="local")))/(table(sc)))
betwenness_block_sc <- c(c(t(z)%*%betweenness(net_Y))/(table(sc)))

# ------------------------------------
# ESBM WITH SUPERVISED GN PRIOR
# ------------------------------------

set.seed(1)

# point estimate
c_Z_GN <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20)
memb_Z_GN <- memb_Z_GN_VI$cl

# compute the betweenness and local transitivity of each group 
# by averaging those of the suspects within that group
z <- dummy(memb_Z_GN)
transivity_block_esbm <- c(c(t(z)%*%(transitivity(net_Y,type="local")))/(table(memb_Z_GN)))
betwenness_block_esbm <- c(c(t(z)%*%betweenness(net_Y))/(table(memb_Z_GN)))


#------------------------------------------------------
# FIGURE
#------------------------------------------------------

# create the dataset containing betweenness and transitivity for the groups obtained by the
# three different methods (first 15 rows refer to the 15 clusters estimated under ESBM)

data_measure <- data.frame(cbind(c(betwenness_block_esbm,betwenness_block_Louv,betwenness_block_sc),c(transivity_block_esbm,transivity_block_Louv,transivity_block_sc)))
rownames(data_measure)<-c(1:dim(data_measure)[[1]])

# define a grouping variable coded as follow:
# 1 to 15: Ids of the groups estimated under ESBM
# 16: Unique Id for the groups estimated under Louvain and spectral clustering
data_measure$group <- as.factor(c(1:15,rep(16,14)))

# define a shape variable coded as follow:
# 1: ESBM groups mostly populated by affiliates according to the pie-charts in Figure 5 
# 2: ESBM groups mostly populated by bosses according to the pie-charts in Figure 5 
# 3: Louvain groups
# 4: Spectral clustering groups
shape_temp <- c(rep(0,15),rep(3,4),rep(4,10))
shape_temp[which(apply(t_standard,1,which.max)<=5)] <- 1
shape_temp[which(apply(t_standard,1,which.max)>5)] <- 2
data_measure$shape <- as.factor(shape_temp) 

# define a size variable which assigns the same dimension to points denoting 
# Louvain and spectral clustering groups, whereas for those corresponding to 
# ESBM groups the size is set proportional to the cardinality of such groups
data_measure$size <- c(table(memb_Z_GN),rep(1.2,14))

# define colors as follow: the color of each ESBM cluster is the one having 
# the highest frequency in the corresponding pie-chart in Figure 5, whereas
# Louvan and spectral clusters are displayed in gray

color_points <- c(mycolors[apply(t_standard,1,which.max)],brewer.pal(9,"Greys")[6]) 

# plot Figure 6
ggplot(data_measure, aes(x=X1, y=X2,color=group,fill=group)) + geom_point(aes(size=(size^1.5),shape=shape),stroke = 1)+theme_bw()+ theme(legend.position = "none")+
   scale_size_continuous(range = c(0,15))+ scale_color_manual(values=color_points)+ scale_fill_manual(values=alpha(color_points,0.3))+xlim(-5,260)+ylim(0.35,1.05)+xlab("Betweenness")+ylab("Local Transivity")+scale_shape_manual(values=c(21, 22,3,4))
```

Refer to Sections 1 and 5 in the article for detailed comments on the above Figures.
