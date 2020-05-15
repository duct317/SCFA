# SCFA
Subtyping via Consensus Factor Analysis (SCFA) can efficiently remove noisy signals from consistent molecular patterns in order to reliably identify cancer subtypes and accurately predict risk scores of patients.
# How to install
- The package can be installed from this repository.
- Install devtools: utils::install.packages('devtools')
- Install the package using: devtools::install_github('duct317/SCFA')
- Install keras package version 2.2.4.1 using: devtools::install_version('keras', version = '2.2.4.1', repos = 'https://cran.rstudio.com/')
- Install all dependencies: matrixStats, foreach, doParallel, igraph, Matrix, caret, clues, cluster, clusterCrit, psych, glmnet
- Install tensorflow and keras in python using: keras::install_keras(tensorflow = "1.13.1")
- For more information about installation of keras, please visit https://keras.rstudio.com/
# Example 
The docker contains the environment and scripts to run an example can be found in the Drive folder at https://tinyurl.com/scfaR.
