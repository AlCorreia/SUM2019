clara2 <- function(data, k, ncat=NULL, metric = "manhattan", 
           stand = FALSE, samples = 5, sampsize = min(nrow(data), 40 + 2 * k),
           trace = 0, medoids.x = TRUE, keep.data = medoids.x, rngR = FALSE,
           pamLike = FALSE, correct.d = TRUE) {

  if(is.null(ncat)) {
    l <- spn.learncats(data)
    ncat <- l$ncat
    data <- l$data
  }

  #create total matrix and fill in -> sum over ncat
  if (metric != "jaccard") ncol <- sum(ncat) - sum(ncat == 2) + sum(ncat == 0)
  else ncol <- ncol(data)
  data.adapt <- data.frame(matrix(ncol=ncol,nrow=nrow(data)))

  cnt <- 1
  for (i in 1:length(ncat)) {
    if (ncat[i] > 2) {
      if (metric == "euclidean") {
        for (j in 1:ncat[i]) {
          data.adapt[,cnt] <- (0.5 * (data[,i] == j) * sqrt(2))
          cnt <- cnt + 1
        }
      } else if (metric == "manhattan") {
        for (j in 1:ncat[i]) {
          data.adapt[,cnt] <- ((data[,i] == j) * 0.5)
          cnt <- cnt + 1
        }
      } else if (metric == "jaccard") {
        ##To be implemented
        data.adapt[,cnt] <- (data[,i] == j)
        cnt <- cnt + 1
      }
    } else {
      data.adapt[,cnt] <- data[,i]
      cnt <- cnt + 1
    }
  }
  return(clara(data.adapt, k, metric, stand, samples, sampsize, trace, medoids.x, keep.data, rngR, pamLike, correct.d))
}
