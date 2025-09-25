# Code for paper "Robust Estimation of Moderated Causal Excursion Odds Ratio in Micro-Randomized Trials"

Code to reproduce results in the paper [Robust Estimation of Moderated Causal Excursion Odds Ratio in Micro-Randomized Trials] by Jiaxin Yu and Tianchen Qian.

Jiaxin Yu
2025.09.24

## 1. Code structure

*  [Simulation1_functions.R](Simulation1_functions.R): the code contains functions including data generating process and nuisance model estimation that are necessary to run the simulations in Section 5.1 of the paper.
  
* [Simulation1.R](Simulation1.R): the code for running simulation study 1 in Section 5.1 of the paper.
  
* [Simulation2_functions.R](Simulation2_functions.R): contains functions including data generating process and nuisance model estimation that are necessary to run the simulations in Section 5.2 of the paper.
  
* [Simulation1.R](Simulation1.R): the code for running simulation study 2 in Section 5.2 of the paper.

* [DataAnalysis_functions.R](DataAnalysis_functions.R): the code contains functions for fitting the proposed estimators in data analysis.

* [DataAnalysis.R](DataAnalysis.R): the code for running the marginal and moderation analyses in the paper. The data used in the analysis can be found [here](https://osf.io/mtcfa).
  
## 2. Results in Paper and Corresponding Code

### 2.1 Simulations

* Simulation study in Section 5.1: the code is in [Simulation1.R](Simulation1.R).

* Simulation study in Section 5.2: the code is in [Simulation2.R](Simulation2.R).

### 2.2 Data Analysis of Drink Less Micro-Randomized Trial

* The Drink Less data is available [here](https://osf.io/mtcfa).

* Marginal and Moderation Analysis: the code is in [DataAnalysis.R](DataAnalysis.R)




