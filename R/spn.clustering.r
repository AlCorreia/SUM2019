library(combinat)
library(clue)

spn.cluster.run <- function(data, nclasses, thr=0.001, verb=FALSE, ncores=1, max.it=10) {
  cat(paste(Sys.time(),'::init',"\n"))
  #Generate the initial class labels
  classcol=ncol(data)
  y_true <- data[, classcol]
  data <- cbind(data[, -classcol], clara2(data[, -classcol], nclasses, samples=20, pamLike=TRUE, metric='manhattan')$clustering)
  cat(paste(Sys.time(),'::start',"\n"))
  cnt <- 1
  while (cnt < max.it) {
      ##print(data[,classcol])
      cat(paste(Sys.time(),":::::SPN ", cnt, "\n"))
      spn <- spn.learn(data, verb=verb, thr=thr, classcol=classcol, class.unif=TRUE)
      # spn.print(spn)
      cat(paste(Sys.time(),":::Predict ", cnt, "\n"))
      res <- spn.predict(spn, data, classcol=classcol, run.rob=FALSE, ncores=ncores, use_memory=FALSE)
      ##print(res[,2])
      if (spn.cluster.convergence(data[,classcol], res[,2], nclasses=nclasses)) break
      data[,classcol] <- res[,2]
      cnt <- cnt + 1
      print(str(cnt))
  }
  res <- spn.predict(spn, data, classcol=classcol, run.rob=TRUE, ncores=ncores, use_memory=FALSE)
  if(cnt >= max.it) cat("Code has not converged within max iterations")
  res <- cbind(y_true, res)
  return(res)
}

spn.cluster.accuracy.perm <- function(real, pred, nclass) {
  permutations <- t(array(unlist(permn(1:nclass)), dim=c(nclass, gamma(nclass+1))))

  max <- -Inf
  for(p in 1:nrow(permutations)) {
    pred.perm <- c()
    for (i in 1:nclass)
      pred.perm[pred == i] <- permutations[p, i]
    sum <- sum(real == pred.perm)
    if (sum > max) max <- sum
  }

  return((max/length(real)))
}

spn.cluster.accuracy.hungarian <- function(real, pred, nclass) {
  contingency <- matrix(0, nrow=nclass, ncol=nclass)
  for (i in 1:length(real)) {
    contingency[real[i], pred[i]] <- contingency[real[i], pred[i]] + 1
  }

  #contingency <- base + table(real, pred)
  perm <- solve_LSAP(as.table(contingency), max=TRUE)

  pred.perm <- c()
  for (i in 1:nclass)
    pred.perm[pred == perm[i]]<- i

  sum <- sum(real == pred.perm)

  return((sum/length(real)))
}

spn.cluster.convergence <- function(old, new, nclasses, thr = 0.001) {
    if (nclasses < 8)
        v <- spn.cluster.accuracy.perm(old, new, nclasses)
    else
        v <- spn.cluster.accuracy.hungarian(old, new, nclasses)
    return(v > 1.0 - thr)
}
