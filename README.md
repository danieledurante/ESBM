# ESBM: Extended stochastic block models

This repository is associated with the article **extended stochastic block models with application to criminal networks** and aims at providing detailed materials and codes to implement the general **ESBM** class presented in the article, and **fully reproduce the results presented in Sections 1, 4 and 5**.

The documentation is organized in **three main folders** described below.  

- `Source`.  It contains all the **source** `R` **functions** [see `esbm.R`] which are required to perform posterior computation and inference under the **ESBM** class, and  the **additional** `cpp` **file** [see `stirling.cpp`] that is necessary to study the **Gibbs-type** priors discussed in the article.

- `Simulation`. It contains the three comprehensive tutorials [see `scenario_1.md`, `scenario_2.md` and `scenario_3.md`] to **fully reproduce step-by-step the results for the simulation scenarios 1, 2 and 3, respectively, presented in Section 4** of the article. The folder contains also the **simulated networks associated with these three scenarios** [see `network_1.RData`, `network_2.RData` and `network_3.RData`].  

- `Application`. It contains a comprehensive tutorial [see `application.md`] to **fully reproduce step-by-step the pre-processing and all the analysis of the *Infinito network*, presented in Sections 1 and 5** of the article. The folder contains also the **pre-processed network which is studied in the article** [see `crime_net.RData`]. Raw data are available at https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndranghetamafia2. 

The analyses are performed with a **iMac (macOS Sierra, version 10.12.6)**, using a `R` version **3.6.1**. 

All the above functions rely on a **basic and reproducible `R` implementation**, mostly meant to provide a clear understanding of the computational routines and steps associated with the proposed model. **Optimized computational routines relying on C++ coding can be easily considered.** Generalizations to include additional priors in the Gibbs-type class, and different types of edges and attributes require minor modifications on the functions in the file `esbm.R`.
