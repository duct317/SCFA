# SCFA
Subtyping via Consensus Factor Analysis (SCFA) can efficiently remove noisy signals from consistent molecular patterns in order to reliably identify cancer subtypes and accurately predict risk scores of patients.
# How to install
- The package can be installed from this repository.
- Install devtools: `utils::install.packages('devtools')`
- Install keras package version 2.2.4.1 using: `devtools::install_version('keras', version = '2.2.4.1', repos = 'https://cran.rstudio.com/')`
- Install tensorflow and keras in python using: `keras::install_keras(tensorflow = "1.10.0")`
- Install the package using: `devtools::install_github('duct317/SCFA')`
- For more information about installation of keras, please visit https://keras.rstudio.com/
# Example 
The docker contains the environment and scripts to run an example can be found at http://scfa.tinnguyen-lab.com.
# How to use the docker:
- Go to directory of the docker image
- Using gunzip to extract the tar file: `gunzip scfa_Docker.tar.gz`
- Load image to docker: `docker load -i scfa_Docker.tar`
- Start a container using this image: `docker run -e PASSWORD=<your-password> -d --name scfa_test -p <port>:8787 scfa`
	- <your-password> is custom password to log in the docker, for example 123456. <port> is an available port of the host, for example, 1234. 
	- The command for above example: `docker run -e PASSWORD=123456 -d --name scfa_test -p 1234:8787 scfa`
	- This command will create a container with rstudio server, this session can be accessed from: <ip-of-host>:1234, in case of personal computer: `localhost:1234`
- The user name for Rstudio is rstudio, password is the password set above.
- Inside this container, the Example.Rmd is used to run an example analysis on GBM dataset.
- If encountering error while running the docker on MacOS or Windows, please increase the memory limit for the docker to 8-16GB. The default setting of 2GB is too low for any analysis to perform well.

# How to use the package for new data 
The package includes these functions:
- SCFA: main function, generating subtyping result. The input is a list of data matrices. In each matrix, rows represent samples and columns represent genes/features.
- SCFA.class: predicting risk scores of test data. The inputs consist of list of train data matrices, train data label and list of test data matrices. 
- The result is reproducible by setting seed for these functions.
- More detail about parameters for each function could be found in the manual.

# Citation:
Tran, D., Nguyen, H., Le, U., Bebis, G., Luu, H. N., & Nguyen, T. (2020). A Novel Method for Cancer Subtyping and Risk Prediction Using Consensus Factor Analysis, <i>Frontiers in Oncology</i>, <i>10</i>, 1052. ([link](https://www.frontiersin.org/article/10.3389/fonc.2020.01052)) 
