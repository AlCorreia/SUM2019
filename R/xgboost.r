library(rlang)
library(tictoc)
library(RWeka)
library(xgboost)
library(stringr)

source('spn.learn.r')


getDummies <- function(data) {
    for(var in colnames(data)){
        data[[var]] <- as.integer(as.factor(data[[var]]))
    }
    return(data)
}


xg.cv <- function(data, classcol=ncol(data), nfolds=nrow(data), nruns=1, verb=FALSE, strat=TRUE, ncores=1, tfname=tempfile('outXG.',tmpdir='.',fileext='.tmp')) {
    n <- nrow(data)
    nfolds <- min(n,nfolds)
    res <- c()
    summ <- spn.learncats(data, classcol=classcol)
    if(summ$ncat[classcol] < 2) stop('class must be discrete')
    for(j in 1:nruns) {
        folds <- 1:n
        if(nfolds < n) {
            if(strat) {
                for(i in unique(data[,classcol])) {
                    nn <- sum(data[,classcol]==i)
                    ind <- rep(1:nfolds, ceiling(nn/nfolds))[1:nn]
                    folds[data[,classcol]==i] <- sample(ind, nn)
                }
            } else {
                ind <- rep(1:nfolds, ceiling(n/nfolds))[1:n]
                folds <- sample(ind, n)
            }
        }
        for(i in 1:nfolds) {
            cur <- which(folds == i)

            if(ncores > 1) {
                cat(paste(Sys.time(),'::learn',"\n"),file=tfname,append=TRUE)
            }
            xgb <- xgboost(data = data.matrix(data[-cur,-ncol(data)]),
             verbose=0,
             label = data[-cur, classcol]-1,
             eta = 0.1,
             max_depth = 15,
             nround=25,
             subsample = 0.5,
             colsample_bytree = 0.5,
             seed = 1,
             eval_metric = "mlogloss",
             objective = "multi:softmax",
             num_class = summ$ncat[classcol],
             nthread = 3
            )
            validation <- cbind(data[cur, classcol], predict(xgb, data.matrix(data[cur, -ncol(data)]))+1)
            res <- rbind(res, validation)
        }
    }
    write.csv(res, paste('Xgboost/', name, '_xg.csv', sep=""))
    return(res)
}



files <- list.files(path="../data/arff", full.names = TRUE)
for(file in files[-c(1:14)]){
    name <- tools::file_path_sans_ext(basename(file))
    print(name)
    data <- read.arff(file)
    data <- getDummies(data)
    xgres <- xg.cv(data)
}
