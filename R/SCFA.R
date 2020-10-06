#' @import keras tensorflow
#' @importFrom BiocParallel SnowParam bplapply
#' @importFrom matrixStats colSums2 rowSds rowMeans2 rowMaxs rowMins colSds
#' @importFrom stats predict rnorm quantile
#' @importFrom psych fa
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom survival coxph
#' @title SCFA
#' @description The main function to perform subtyping. It takes a list of data matrices as the input and outputs the subtype for each patient.
#' @param dataList List of data matrices. In each matrix, rows represent samples and columns represent genes/features.
#' @param k Number of clusters, leave as default for auto detection.
#' @param max.k Maximum number of cluster
#' @param ncores Number of processor cores to use.
#' @param seed Seed for reproducibility, you still need to use set.seed function for full reproducibility.
#' @return A numeric vector containing cluster assignment for each sample.
#' @examples
#' #Load example data (GBM dataset)
#' data("GBM")
#' #List of one matrix (microRNA data)
#' dataList <- GBM$data
#' #Survival information
#' survival <- GBM$survival
#' library(survival)
#' #Generating subtyping result
#' set.seed(1)
#' subtype <- SCFA(dataList, seed = 1, ncores = 2L)
#' #Perform survival analysis on the result
#' coxFit <- coxph(Surv(time = Survival, event = Death) ~ as.factor(subtype), data = survival, ties="exact")
#' coxP <- round(summary(coxFit)$sctest[3],digits = 20)
#' print(coxP)
#' @export
SCFA <- function(dataList, k = NULL, max.k = 5, ncores = 10L, seed = NULL) {
    gen.fil <- TRUE
    all_data <- NULL
    all_clus <- NULL
    all_latent <- NULL
    counter <- 1
    for (data in dataList) {
        if (max(data) <= 1)
            data <- 10^data - 1

        if (ncol(data) > 50000) {
            col_mean_data = colMeans(data)
            idx <- order(col_mean_data, decreasing = TRUE)[seq(50000)]
            data <- data[, idx]
        }

        tmp <- SCFA.basic(data, k = k, max.k = max.k, ncores = ncores, gen.fil = gen.fil, seed = seed)


        all_data <- cbind(all_data, tmp$filter)
        all_clus <- c(all_clus, tmp$all.res)
        all_latent <- c(all_latent, tmp$all.latent)

        counter <- counter + 1
    }

    if (length(dataList) > 1) {
        tmp <- SCFA.basic(all_data, k = k, max.k = max.k, ncores = ncores, gen.fil = gen.fil, seed = seed)
        all_clus <- c(all_clus, tmp$all.res)
        all_latent <- c(all_latent, tmp$all.latent)
    }

    tmp1 <- list()
    tmp1$all <- all_clus
    cluster <- clustercom2(tmp1)
    return(as.numeric(cluster))
}


gene.filtering <- function(data, original_dim, batch_size, ncores.ind, seed)
{
    or <- list()
    param <- SnowParam(workers = 3, RNGseed = seed)
    FUN <- function(i, data, original_dim, batch_size, seed)
    {
        if (is.null(seed)) {
            config <- list()
            config$intra_op_parallelism_threads <- ncores.ind
            config$inter_op_parallelism_threads <- ncores.ind
            session_conf <- do.call(tf$ConfigProto, config)
            sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
            k_set_session(session = sess)
        } else {
            use_session_with_seed((seed + i))
        }

        if (nrow(data) > 2000) {
            ind <- sample.int(nrow(data), 2000, replace = FALSE)
            data.tmp <- data[ind, ]
            batch_size <- round(length(ind)/50)
        } else {
            data.tmp <- data
        }
        x <- layer_input(shape = c(original_dim))
        h <- layer_dense(x, 25, kernel_constraint = constraint_nonneg())
        x_decoded_mean <- layer_dense(h, original_dim)
        vae <- keras_model(x, x_decoded_mean)
        vae <- compile(vae, optimizer = tf$contrib$opt$AdamWOptimizer(1e-04, 0.001), loss = "mse")
        his <- fit(vae, data.tmp, data.tmp, shuffle = TRUE, epochs = 5, batch_size = batch_size, verbose = 0)
        W <- get_weights(get_layer(vae, index = 2))[[1]]
        Wsd <- rowSds(W)
        Wsd[is.na(Wsd)] <- 0
        Wsd <- (Wsd - min(Wsd))/(max(Wsd) - min(Wsd))
        Wsd
    }
    suppressMessages(suppressWarnings(
        or <- bplapply(seq(3), FUN, BPPARAM = param, data = data, original_dim = original_dim, batch_size = batch_size, seed = seed)
    ))
    or <- rowMeans2(as.matrix(data.frame(or)))
    or
}

generating.latent <- function(da, ncores, classification, seed)
{
    if(classification) re <- rep(10:15, 3) else re <- rep(5:10, 3)
    param <- SnowParam(workers = min(10, ncores), RNGseed = seed)
    FUN <- function(counter, da, re)
    {
        RhpcBLASctl::blas_set_num_threads(1)
        i <- re[counter]
        tmp <- da * matrix(rnorm(da, sd = 0.02, mean = 1), ncol = ncol(da))
        fit <- fa(t(tmp), nfactors = i, rotate = "varimax", scores = "tenBerge")
        fa <- fit$loadings[, seq(i)]
        fa
    }
    suppressMessages(suppressWarnings(
        latent <- bplapply(seq(length(re)), FUN, BPPARAM = param, da = da, re = re)
    ))
    latent
}

clustering.all <- function(latent, ncores, k, max.k, seed)
{
    result <- list()
    param <- SnowParam(workers = min(length(latent), ncores), RNGseed = seed)
    FUN <- function(x, k, max.k)
    {
        RhpcBLASctl::blas_set_num_threads(1)
        cluster <- clus(x, k = k, max.k = max.k)
        cluster
    }
    suppressMessages(suppressWarnings(
        result$all <- bplapply(latent, FUN, BPPARAM = param, k = k, max.k = max.k)
    ))
    result
}

SCFA.basic <- function(data = data, k = NULL, max.k = 5, ncores = 10L,
                       gen.fil = TRUE, classification = FALSE, seed = NULL) {
    non.zero.prop <- colSums2(data != 0)/nrow(data)
    data <- data[, non.zero.prop > 0]

    tmp.max <- rowMaxs(data)
    tmp.min <- rowMins(data)
    data <- (data - tmp.min)/(tmp.max - tmp.min)

    ncores.ind <- as.integer(max(1, floor(ncores/3)))
    original_dim <- ncol(data)
    batch_size <- max(round(nrow(data)/50), 2)
    gen.fil <- (gen.fil & (ncol(data) > 5000))
    n <- ifelse(gen.fil, min(5000, ncol(data)), ncol(data))

    # Feature selection
    if (gen.fil) {
        or <- gene.filtering(data = data, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, seed = seed)
        keep <- intersect(which(or > quantile(or, (1 - min(n, original_dim)/original_dim))), which(or > 0))
        da <- data[, keep]
        original_dim_reduce <- ncol(da)
        or.da <- da
    } else {
        da <- data
        original_dim_reduce <- ncol(da)
        or.da <- da
    }
    latent <- generating.latent(da = da, ncores = ncores, classification = classification, seed = seed)
    if(!classification)
    {
        result <- clustering.all(latent = latent, ncores = ncores, k = k, max.k = max.k, seed = seed)
        final <- clustercom2(result)
        g.en <- latent[[which.max(vapply(result$all, function(x) adjustedRandIndex(x, final), numeric(1)))]]
        list(cluster = final, latent = g.en, all.latent = latent,
             filter = or.da, keep = colnames(or.da), all.res = result$all)
    }
    else
    {
        list(all.latent = latent, filter = or.da)
    }
}


#' @import keras tensorflow
#' @importFrom matrixStats colSums2 rowSds rowMeans2 rowMaxs rowMins colSds
#' @importFrom stats predict rnorm quantile
#' @importFrom psych fa
#' @importFrom glmnet cv.glmnet
#' @title SCFA.class
#' @description Perform risk score prediction on input data. This function requires training data with survival information. The output is the risk scores of patients in testing set.
#' @param dataListTrain List of training data matrices. In each matrix, rows represent samples and columns represent genes/features.
#' @param trainLabel Survival information of patient in training set in form of Surv object.
#' @param dataListTest List of testing data matrices. In each matrix, rows represent samples and columns represent genes/features.
#' @param ncores Number of processor cores to use.
#' @param seed Seed for reproducibility, you still need to use set.seed function for full reproducibility.
#' @return A vector of risk score predictions for patient in test set.
#' @examples
#' #Load example data (GBM dataset)
#' data("GBM")
#' #List of one matrix (microRNA data)
#' dataList <- GBM$data
#' #Survival information
#' survival <- GBM$survival
#' library(survival)
#' #Split data to train and test
#' set.seed(1)
#' idx <- sample.int(nrow(dataList[[1]]), round(nrow(dataList[[1]])/2) )
#' survival$Survival <- survival$Survival - min(survival$Survival) + 1 # Survival time must be positive
#' trainList <- lapply(dataList, function(x) x[idx, ] )
#' trainSurvival <- Surv(time = survival[idx,]$Survival, event =  survival[idx,]$Death)
#' testList <- lapply(dataList, function(x) x[-idx, ] )
#' testSurvival <- Surv(time = survival[-idx,]$Survival, event =  survival[-idx,]$Death)
#' #Perform risk prediction
#' result <- SCFA.class(trainList, trainSurvival, testList, seed = 1, ncores = 2L)
#' #Validation using concordance index
#' c.index <- concordance(coxph(testSurvival ~ result))$concordance
#' print(c.index)
#' @export
SCFA.class <- function(dataListTrain, trainLabel, dataListTest, ncores = 10L, seed = NULL) {
    alpha <- 0.5
    nfold <- 5
    dataList <- list()
    for (i in seq(length(dataListTrain))) {
        dataList[[i]] <- rbind(dataListTrain[[i]], dataListTest[[i]])
    }
    all_data <- NULL
    all_latent <- NULL
    counter <- 1
    for (data in dataList) {
        if (max(data) <= 1)
            data <- 10^data - 1
        if (ncol(data) > 50000) {
            col_mean_data <- colMeans(data)
            idx <- order(col_mean_data, decreasing = TRUE)[seq(50000)]
            data <- data[, idx]
        }
        tmp <- SCFA.basic(data, ncores = ncores, classification = TRUE, seed = 1)
        all_data <- cbind(all_data, tmp$filter)
        all_latent <- c(all_latent, tmp$all.latent)
        counter <- counter + 1
    }
    if (length(dataList) > 1) {
        tmp <- SCFA.basic(all_data, ncores = ncores, classification = TRUE, seed = 1)
        all_latent <- c(all_latent, tmp$all.latent)
    }
    train.idx <- seq(nrow(dataListTrain[[1]]))
    all.predict <- matrix(ncol = length(all_latent), nrow = nrow(dataListTest[[1]]))
    for (i in seq(length(all_latent))) {
        data <- all_latent[[i]]
        train.x <- data[train.idx, ]
        test.x <- data[-train.idx, ]
        suppressWarnings(m <- cv.glmnet(train.x, trainLabel, family = "cox", nfolds = nfold, standardize = FALSE,
            alpha = alpha))
        pred <- predict(m, test.x, s = "lambda.min", type = "response")
        all.predict[, i] <- pred
    }
    if (length(which(colSds(all.predict) != 0)) >= 1)
        all.predict <- all.predict[, which(colSds(all.predict) != 0)]
    if (!is.null(dim(all.predict)))
        final.pred <- exp(rowMeans(log(all.predict), na.rm = TRUE)) else final.pred = all.predict
    return(final.pred)

}

#' @title GBM
#'
#' @description GBM dataset, including microRNA and survidal data.
#'
#' @format A list with two items:
#' \describe{
#'     \item{data}{List of microRNA data matrix.}
#'     \item{survival}{Survival information.}
#' }
"GBM"

