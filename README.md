# ESBM: Extended stochastic block models for networks

This repository is associated with the article [**extended stochastic block models**](https://github.com/danieledurante/ESBM) and aims at providing detailed materials and codes to perform posterior computation and inference under the general **ESBM** class presented in the article.

The documentation is organized in two main parts described below.  

- [`Data and Codes`](https://github.com/danieledurante/ESBM/tree/master/Data%20and%20Codes).  It contains [i] useful data to reproduce the [`Tutorial.md`](https://github.com/danieledurante/ESBM/blob/master/Tutorial.md), [ii]  commented source `R` functions in [`esbm.R`](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/esbm.R) to perform posterior computation and inference under the **ESBM** class, and [iii]  the additional `cpp` file [`stirling.cpp`](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/stirling.cpp) that is necessary to study the **Gibbs-type** priors discussed in the article.

- [`Tutorial.md`](https://github.com/danieledurante/ESBM/blob/master/Tutorial.md). It contains a comprehensive tutorial to perform posterior computation and inference under the ESBM, leveraging the methods and algorithms presented in the article and implemented in the source code [`esbm.R`](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/esbm.R). To accomplish this goal, we reproduce step-by-step the analysis of `network 1` in the simulation study within the article.

All the above analyses are performed with a **MacBook Pro (OS X El Capitan, version 10.11.6)**, using a `R` version **3.6.3**. 

All the above functions rely on a **basic and reproducible `R` implementation**, mostly meant to provide a clear understanding of the computational routines and steps associated with the proposed model. **Optimized computational routines relying on C++ coding can be easily considered.** Generalizations to include additional priors in the Gibbs-type class, and different types of edges and attributes require minor modifications on the functions in the file [`esbm.R`](https://github.com/danieledurante/ESBM/blob/master/Data%20and%20Codes/esbm.R).

The **bill co-sponsorship network** considered in the application is openly available at [`https://github.com/briatte/parlnet`](https://github.com/briatte/parlnet). More specifically, we study a dichotomized version of the network `net_it_ca2013` within the data frame [`parlnet.rda`](https://github.com/briatte/parlnet/blob/master/parlnet.rda).
