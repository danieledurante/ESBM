# ESBM: Extended stochastic block models for networks

This repository is associated with the paper [**Extended stochastic block models**]() and aims at providing detailed materials and codes to perform posterior computation and inference under the general **ESBM** class. 

The documentation is organized in two main parts described below.  

- [`Data and Codes`]().  It contains [i] useful data to reproduce the [`Tutorial.md`](), [ii] commented source `R` functions in [`esbm.R`]() to perform posterior computation and inference under the **ESBM**, and [iii] the additional `cpp` file [`stirling.cpp`]() that is necessary to study the **Gibbs-type** priors discussed in the article.

- [`Tutorial.md`](). It contains a comprehensive tutorial to perform posterior computation and inference under the ESBM leveraging the methods and algorithms, presented in the article and implemented in the source code [`esbm.R`](). To accomplish this goal, we reproduce step-by-step the analysis of NETWORK 1 in the simulation study within the article.

All the above analyses are performed with a **MacBook Pro (OS X El Capitan, version 10.11.6)**, using a `R` version **3.6.3**. 

All the above functions rely on a **basic and reproducible `R` implementation**, mostly meant to provide a clear understanding of the computational routines and steps associated with the proposed model. **Optimized computational routines relying on C++ coding can be easily considered.** Generalizations to include additional priors in the Gibbs-type class and different types of edges and attributes require minor modifications on the functions in the file [`esbm.R`]().

The bill co-sponsorship network considered in the application is openly available at [`https://github.com/briatte/parlnet`](https://github.com/briatte/parlnet). More specifically, we studied a dichotomized version of the network `net_it_ca2013` with the data frame [`parlnet.rda`](https://github.com/briatte/parlnet/blob/master/parlnet.rda).
