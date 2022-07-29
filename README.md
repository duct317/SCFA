# SCFA
Subtyping via Consensus Factor Analysis (SCFA) can efficiently remove noisy signals from consistent molecular patterns in order to reliably identify cancer subtypes and accurately predict risk scores of patients.
# How to install
- The package can be installed from this repository.
- Install devtools: `utils::install.packages('devtools')`
- Install the package using: `devtools::install_github('duct317/SCFA')`
- When the package is loaded, it will check for the necessary `libtorch`: `library(SCFA)`  
  `libtorch` can be installed using: `torch::install_torch()`
# Example 
Please follow the vignette to run the package on example data.
# Citation:
Tran, D., Nguyen, H., Le, U., Bebis, G., Luu, H. N., & Nguyen, T. (2020). A Novel Method for Cancer Subtyping and Risk Prediction Using Consensus Factor Analysis, <i>Frontiers in Oncology</i>, <i>10</i>, 1052. ([link](https://www.frontiersin.org/article/10.3389/fonc.2020.01052)) 
