#' @import doParallel
#' @importFrom stats kmeans
#' @importFrom foreach %dopar% foreach %do%
clus <- function(data, k = NULL, max.k = 5)
{
  if (is.null(k)) k <- nclusterPar(data, max.k = max.k)
  if (nrow(data) < 1000 & (k <= 5))
  {
    k <- kmeans(data, k, nstart = 100, iter.max = 1e6)
    k$cluster
  } else
  {
    kknn <- try(specClust(data, k, nn = 7))
    while (class(kknn) == "try-error") {
      k <- k +1
      kknn <- try(specClust(data, k, nn = 7))
    }
    kknn$cluster
  }

}


nclusterPar <- function(data, max.k = 5)
{
  result <- foreach (j = 1:10) %do% {
    set.seed(j)
    idx <- sample(1:nrow(data),min(500, nrow(data)))
    to.test <- matrix(0, nrow = 10, ncol = 3)
    e <- 3
    for (i in 2:10) {

      kknn <- try(specClust(data[idx,], i, nn = 7))
      if (class(kknn) != "try-error")
      {
        to.test[i,1] <- kknn$betweenss/kknn$totss
        to.test[i,2] <- kknn$tot.withinss
      } else {
        e <- i + 2
      }


    }
    for (i in e:10) {
      to.test[i,3] <- (to.test[i,2]-to.test[i-1,2])/to.test[i-1,2]
    }
    c(which.max(to.test[,1]), which.max(to.test[,3]) )
  }
  result <- t(data.frame(result))
  min(max.k, floor(mean(result, na.rm = T)+0.5))

}


clustercom2 <- function(result)
{
  test <- matrix(0, ncol = length(result$all),nrow = length(result$all))
  for (i in 1:length(result$all)) {
    for (j in 1:length(result$all)) {
      if (i != j)  test[i,j] <- adjustedRandIndex(result$all[[i]],result$all[[j]])
    }
  }
  for (i in 1:length(result$all)) {
    test[i,i] <- mean(test[-i,i])
  }

  found <- FALSE
  if(sum(test < 0.7) >0)
  {
    i <- 2
  } else
  {
    i <- 1
  }
  while(!found )
  {
    k <- kmeans(test,i,nstart = 100, iter.max = 5000)$cluster
    max <- 0
    for (c in unique(k)) {
      score <- mean(test[which(k == c), which(k == c)])
      if (score > max  & length(which(k == c)) > 1)
      {
        max <- score
        idx <- which(k == c)
      }
    }
    if (max > 0.8) found <- TRUE
    if (i > 3) found <- TRUE
    i <- i + 1
  }



  res <- t(data.frame(result$all))
  res <- res[idx,]
  cl.max <- floor(mean(apply(res, 1, function(x) length(unique(x)))) )

  res <- data.frame(result$all)

  da = apply(res,1, paste, collapse="#")
  indUnique = which(!duplicated(da))
  indAll = match(da, da[indUnique])

  da <- res[indUnique,]

  test <- wMetaC(da, (cl.max+1), hmethod = "ward.D", minN.cluster = cl.max)
  test$finalC[indAll]


}

adjustedRandIndex <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}


#' @import cluster
#' @importFrom utils combn
#' @import clusterCrit
#' @importFrom stats as.dist cutree hclust median
#' @import Matrix
wMetaC <- function(nC, hmethod, enN.cluster, minN.cluster, maxN.cluster, sil.thre,
                   height.Ntimes) {
  # This is to obtain the weight matrix for each cluster solution for following
  # meta-clustering
  N = nrow(nC)  #number of points
  C = ncol(nC)  #number of clustering methods/times; or K

  AA = Reduce("+", apply(nC, 2, getA))  #sum over obtained matrices; execute along the column and then matrix sum
  AA = AA/C

  indA = Matrix::which(AA != 0, arr.ind = T)  #find non-zero indices of AA
  nd = vapply(AA[indA], function(x) x * (1 - x), numeric(1))

  newAA = sparseMatrix(i = indA[, 1], j = indA[, 2], x = nd, dims = c(N, N))

  w0 = 4/N * Matrix::rowSums(newAA)  #the weight for each point

  e = 0.01
  w1 = (w0 + e)/(1 + e)  #adjusted point weight

  x = as.vector(sapply(1:C, function(i) {
    paste(nC[, i], "_", i, sep = "")
  }))  #convert the matrix (N*C) to vector (concatenating them)

  newnC <- matrix(x, nrow = N, byrow = FALSE)  #reshape a vector to a matrix; by column

  R = unique(x)  #all unique labels
  allC = length(R)  #number of all unique labels

  cb = combn(allC, 2)  #all possible combinations (n*(n-1)/2)
  alls = apply(cb, 2, getss, R = R, x = x, w1 = w1)  #calculate the weight s for all combinations


  S0 = sparseMatrix(i = cb[1, ], j = cb[2, ], x = alls, dims = c(allC, allC))  #triangle part of the S
  S = S0 + t(S0) + diag(allC)


  if (missing(sil.thre)) {
    sil.thre = 0
  }
  hres = get_opt_hclust(S, hmethod, N.cluster = enN.cluster, minN.cluster, maxN.cluster,
                        sil.thre, height.Ntimes)  #solely using the silhouette index as the criteria

  tf = hres$f
  v = hres$v
  cat("The number of clusters before voting is: ", hres$optN.cluster, "\n")

  newnC[] <- vapply(newnC, function(q) tf[match(q, R)], numeric(1))  #apply to every element; reorganizing the clusters for different results

  finalC = apply(newnC, 1, function(d) names(sort(table(d), decreasing = TRUE)[1]))  #find the most repeated elements for each row

  N.cluster = length(unique(finalC))  #note that the number of clusters for meta-clustering is not determined by previous selection, but by the unique number in the final round.

  perc = 0.5
  if(N.cluster == 1){#better not to have only one cluster
    finalC = apply(newnC, 1, function(d){
      x = sort(table(d), decreasing = TRUE)[1:2]
      n0 = length(x[1])
      if(x[2] >= n0*perc){
        y = names(x[2])
      }else{
        y = names(x[1])
      }
      return(y)
    })

    N.cluster = length(unique(finalC))
  }
  cat("The optimal number of clusters for ensemble clustering is:", N.cluster,
      "\n")


  # For ease of visualization
  uC = unique(finalC)#unique clusters

  y0 = apply(newnC, 1, function(q){
    t = rep(0, N.cluster)
    for(i in c(1:N.cluster)){
      t[i] = length(which(q %in% uC[i]))
    }
    return(t)
  })#need to reorganize before counting
  #   print(dim(y0))
  y0 = t(y0)#transpose


  x0 = matrix(0, nrow = N, ncol = N.cluster)
  #   print(dim(x0))


  tw = 0.5
  #   print(uC)
  for(i in 1:N){
    xind = which(finalC[i]==uC)
    x0[i, xind] = 1#the correct clustering result
    allind = which(y0[i,]!=0)#all the counts
    diffind = setdiff(allind, xind)#some other counts which are not the correct cluster
    if(length(diffind) != 0){
      x0[i, diffind] = tw* y0[i, diffind]/y0[i, xind]#use a reduced weight
    }
  }


  out = list()  #declare
  out$finalC = finalC
  out$x0 = x0
  return(out)
}




getA <- function(rowColor) {
  # This is to obtain the weighted co-association matrix for clustering solution
  # rowColor
  N = length(rowColor)  #number of points

  L = levels(factor(rowColor))

  # find indices for each cluster, then all combinations of indices
  tmp = sapply(L, function(k) {
    r = which(rowColor %in% k)
    expand.grid(r, r)
  })

  # reshape to the indices
  allind = matrix(unlist(t(tmp)), ncol = 2, byrow = F)  #need transpose
  A = sparseMatrix(i = allind[, 1], j = allind[, 2], x = 1, dims = c(N, N))  #non-zero entries
  return(A)
}




getss <- function(pind, R, x, w1) {
  # This is to get the element of S
  pairk = lapply(pind, getnewk, R = R, x = x, N = length(w1))  #run for two indices

  intset = intersect(unlist(pairk[1]), unlist(pairk[2]))  #set intersection

  ss = 0
  if (length(intset) != 0) {
    uset = union(unlist(pairk[1]), unlist(pairk[2]))  #set union
    ss = sum(w1[intset])/sum(w1[uset])
  }
  return(ss)
}

getnewk <- function(k, R, x, N) {
  # This is to get the original index of the sample
  k1 = which(x %in% R[k])  #find samples with k-th cluster
  d1 = unlist(strsplit(R[k], "_"))  #the name contains only two parts; get the numbering part
  d = as.numeric(Matrix::tail(d1, n = 1))  #the last element of the split arrays
  newk1 = k1 - (d - 1) * N  #the index
  return(newk1)
}

get_opt_hclust <- function(mat, hmethod, N.cluster, minN.cluster, maxN.cluster, sil.thre,
                           height.Ntimes) {
  # if no agglomeration method for hierarchical clustering is provided
  if (missing(hmethod) || is.null(hmethod)) {
    hmethod = "ward.D"  #the default hierarchical clustering agglomeration method is 'ward.D'
  }

  # if no minimum number of clusters is provided
  if (missing(minN.cluster) || is.null(minN.cluster)) {
    minN.cluster = 2  #by default, we try the minimum number of clusters starting from 2
  }

  # if no maximum number of clusters is provided
  if (missing(maxN.cluster) || is.null(maxN.cluster)) {
    maxN.cluster = 40  #by default, we try the maximum number of clusters as large as 40 or the number of cells minus 1, whichever is smaller.
  }

  # if no threshold for the maximum Silhouette index is provided
  if (missing(sil.thre) || is.null(sil.thre)) {
    sil.thre = 0.35  #by default, we use 0.35 to determine whether we use Silhouette index as the criteria to determine the optimal number of clusters
  }

  # if no threshold for the height difference is provided
  if (missing(height.Ntimes) || is.null(height.Ntimes)) {
    height.Ntimes = 2  #by default, we select the first height which is (height.Ntimes) times larger than the immediate consecutive height
  }

  # just use simple criteria to determine whether they are feature vectors or
  # similarity matrix, and then we use different ways to measure the distance
  if (Matrix::isSymmetric(mat)) {
    # symmmetric matrix
    d = as.dist(1 - mat)
    flag1 = 1
  } else {
    d = as.dist(1 - cor(t(mat)))
    flag1 = 0
  }

  h = hclust(d, method = hmethod)  #ward to ward.D


  # if N.cluster is given, we simply use the given N.cluster for hierarchical
  # clustering
  #     if (!missing(N.cluster) && is.numeric(N.cluster)) {
  if (is.numeric(N.cluster)) {
    if (!is.numeric(N.cluster)) {
      stop("The given N.cluster is not a numeric!")
    }
    if (N.cluster%%1 != 0) {
      stop("The given N.cluster is not an integer!")
    }
    if (N.cluster < 2) {
      stop("The given N.cluster is less than 2, which is not suitable for clustering!")
    }

    v = cutree(h, k = N.cluster)  #for different numbers of clusters
    f = v  #the optimal clustering results
    sil = silhouette(v, d)
    msil = median(sil[, 3])
    ch0 = intCriteria(data.matrix(mat), as.integer(v), "Calinski_Harabasz")
    CHind = unlist(ch0, use.names = F)  #convert a list to a vector/value
    optN.cluster = N.cluster
  }


  hres = list()
  hres$f = f#optimal clustering results
  hres$v = v#different numbers of clustering results
  hres$maxsil = max(msil)
  hres$msil = msil
  hres$CHind = CHind
  hres$height = h$height
  hres$optN.cluster = optN.cluster

  return(hres)
}

#' @importFrom igraph arpack decompose graph
#' @importFrom stats cor dnorm qnorm
fast.table <- function (data)
{
  if(!is.data.frame(data))
    data = as.data.frame(data, stringsAsFactors = FALSE)
  da = do.call("paste", c(data, sep = "\r"))
  ind = !duplicated(da)
  levels = da[ind]
  cat <- factor(da,levels = levels)
  nl <- length(levels(cat))
  bin <- (as.integer(cat) - 1)
  pd <- nl
  bin <- bin[!is.na(bin)]
  if (length(bin)) bin <- bin + 1
  y <- tabulate(bin, pd)
  result=list(index = bin, weights = y, data = data[ind,])
  result
}


createFolds <- function(vec, k)
{
  tmp <- round(seq(1, max(vec), length.out = k+1))
  res <- list()
  for (i in 1:k) {
    res[[i]] <- tmp[i]:tmp[i+1]
  }
  res
}

mydist <- function(data, k = 20, distance = 2)
{
  m <- dim(data)[1]
  q <- dim(data)[2]
  D <- matrix(nrow = m, ncol = k)
  C <- matrix(nrow = m, ncol = k)
  folds <- createFolds(1:m, k = ceiling(m/1000))

  for (idx in folds) {
    tmp <- data[idx,]
    dis.tmp <- 1 - cor(t(tmp), t(data))

    for (i in 1:length(idx)) {
      tmp1 <- order(dis.tmp[i,])
      tmp1 <- tmp1[2:(k+1)]

      D[idx[i],] <- dis.tmp[i,tmp1]
      C[idx[i],] <- tmp1
    }



  }

  list(D,C)

}


getClosest = function(X, Y){
  m = nrow(Y)
  n = nrow(X)
  res = matrix(0, n, m)
  for(i in 1:m){
    tmp = (X-rep(Y[i,], each=n))**2
    res[,i] = rowSums(tmp)
  }
  apply(res, 2, which.min)
}


Laplacian <- function(DC, k, normalize="none"){
  normalize <- match.arg(normalize, c("none", "symmetric", "random-walk"))
  m <- dim(DC[[1]])[1]
  INDEX = matrix(c( rep(1:m,k), as.vector(DC[[2]])) , ncol=2)
  ind = which(INDEX[,2] < INDEX[,1])
  INDEX[ind, ] = INDEX[ind, c(2,1), drop=FALSE]
  INDEX2 = fast.table(INDEX)
  ind = which(!duplicated(INDEX2[[1]]))
  INDEX =INDEX2[[3]]
  i = c(INDEX[,1],INDEX[,2])
  j = c(INDEX[,2],INDEX[,1])
  X = as.vector(DC[[1]])[ind]
  x =  c(X, X)
  # graph.laplacian ??
  result <- sparseMatrix(i = i, j = j, x=x, dims = c(m,m))
  D = Matrix::rowSums(result)
  if(normalize=="none") return(Diagonal(x=D) - result)
  if(normalize=="symmetric"){
    TMP = Diagonal(x=1/sqrt(D))
    result = TMP %*% result %*% TMP
    return(Diagonal(m) - result)
  }
  if(normalize=="random-walk"){
    return(Diagonal(m) - Diagonal(x=1/D)%*%result)
  }
  result
}


AUC = function(y){
  l = length(y)
  x = 0:(l-1)
  y = y - y[1]
  res = numeric(0)

  for(i in 1:l){
    A = 0
    A = y[i]*(i-1)/2
    B = y[i] * (l-i)
    C = (y[l] - y[i]) *  (l-i) / 2
    res[i] = A+B+C
  }
  res
}




specClust <- function (data, centers=NULL, nn = 7, method = "symmetric", gmax=NULL, ...)
{
  call = match.call()
  if(is.data.frame(data)) data = as.matrix(data)
  # unique data points
  da = apply(data,1, paste, collapse="#")
  indUnique = which(!duplicated(da))
  indAll = match(da, da[indUnique])

  data2 = data
  data  = data[indUnique, ]
  n <- nrow(data)

  #data = scale(data, FALSE, TRUE)


  if(is.null(gmax)){
    if(!is.null(centers)) gmax = centers - 1L
    else gmax = 1L
  }
  test=TRUE
  DC.tmp = mydist(data, 30)
  while(test){

    if(nn > ncol(DC.tmp[[1]])) DC.tmp = mydist(data, nn*2)

    DC = list(DC.tmp[[1]][,1:nn], DC.tmp[[2]][,1:nn])
    sif <- rbind(1:n, as.vector(DC[[2]]))
    g <- graph(sif, directed=FALSE)
    g <- decompose(g, min.vertices=4)
    if (length(g) > 1) {
      #warning("graph not connected")
      if(length(g)>=gmax) nn = nn+2
      else test=FALSE
    }
    else test=FALSE
  }

  W <- DC[[1]]
  n <- nrow(data)
  wi <- W[,nn]
  SC <- matrix(1, nrow(W), nn)
  SC[] <-  wi[DC[[2]]] * wi
  W = W^2 / SC

  alpha=1/(2*(nn+1))
  qua=abs(qnorm(alpha))
  W = W*qua
  W = dnorm(W, sd = 1)

  DC[[1]] = W
  L = Laplacian(DC, nn, method)

  f <- function(x, extra) as.vector(extra %*% x)

  if(is.null(centers))kmax = 25
  else kmax = max(centers)

  U <- arpack(f, extra = L, options = list(n = n, which = "SM",
                                           nev = kmax, ncv = 2 * kmax, mode=1), sym = TRUE)
  ind <- order(U[[1]])
  U[[2]] = U[[2]][indAll, ind]
  U[[1]] = U[[1]][ind]
  if (is.null(centers)) {
    tmp = which.max(diff(U[[1]]))+1
    centers = which.min(AUC(U[[1]][1:tmp]))
  }
  if(method == "symmetric"){
    rs = sqrt(rowSums(U[[2]]^2))
    U[[2]] =  U[[2]]/rs
  }
  #result = kmeans(U[[2]], centers = centers, nstart = 50, iter.max = 100, ...)
  result = kmeans(scale(U[[2]], center = T), centers = centers, nstart = 500, iter.max = 1000, ...)
  archeType = getClosest(U[[2]][indAll, ], result$centers)
  result$eigenvalue = U[[1]]
  result$eigenvector = U[[2]]
  result$data = data2
  result$indAll = indAll
  result$indUnique = indUnique
  result$L = L
  result$archetype = archeType
  result$call = call
  class(result) = c("specClust", "kmeans")
  result
}
