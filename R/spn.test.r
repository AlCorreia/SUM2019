## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.
require("parallel")
require("foreach")
require("doParallel")
source('spn.predict1.r')
source('spn.predict2.r')

# Wrapper for spn.predict2
spn.predict <- function(spn, data, classcol=ncol(data), verb=FALSE, run.rob=FALSE,ncores=detectCores()-1,tfname=tempfile('outSPN.',tmpdir='.',fileext='.tmp'), use_memory, len.missing=0) {
    nclass <- spn$ncat[classcol]
    if(nclass < 2) stop('class must be discrete')
    nr <- nrow(data)
    if(ncores > 1) {
        cat(paste(Sys.time(),'::setup',ncores,"cores\n"),file=tfname,append=TRUE)
        cl <- makeCluster(ncores)
        registerDoParallel(cl, cores = ncores)
        res <- foreach(i = 1:nr, .packages = c("cluster"), .combine=rbind) %dopar% {
            if(i %% 20 == 1)
                cat(paste(Sys.time(),":::Predict",i-1,"\n"),file=tfname,append=TRUE)
            spn.predict2(spn, as.numeric(data[i,]), classcol=classcol, verb=verb, run.rob=run.rob, use_memory=use_memory, len.missing=len.missing)
        }
        stopCluster(cl)
    } else {
        res <- c()
        for(i in 1:nr) {
            res <- rbind(res, spn.predict2(spn, as.numeric(data[i,]), classcol=classcol, verb=verb, run.rob=run.rob, use_memory=use_memory, len.missing=len.missing))
        }
    }
    logprobNames <- c()
    for(i in 1:nclass){
      logprobNames <- c(logprobNames, paste("logprob", i))
    }
    nloops <- ncol(res) - nclass - 6
    loopNames <- c()
    for(i in 1:nloops){
      loopNames <- c(loopNames, paste("loop", i))
    }
    if(run.rob)
      colnames(res) <- c("y_true", "y_pred", "prob1-prob2", "robustness", logprobNames, "count", "time", loopNames)
    else {
      colnames(res) <- c("y_true", "y_pred", "prob1-prob2", "robustness", logprobNames, "count", "time")
    }
    return(res)
}

spn.cv <- function(data,
                   classcol=ncol(data),
                   nfolds=nrow(data),
                   thr=0.01,
                   nruns=1,
                   verb=FALSE,
                   strat=TRUE,
                   height=1000000,
                   run.rob=TRUE,
                   ncores=detectCores()-1,
                   tfname=tempfile('outSPN.',tmpdir='.',fileext='.tmp'),
                   use_memory=FALSE,
                   len.missing=0,
                   seed=NULL,
                   path=NULL) {
    if(!is.null(seed)) set.seed(seed)
    n <- nrow(data)
    nfolds <- min(n,nfolds)
    res <- c()
    res_mem <- c()
    if(!is.null(classcol)){
      summ <- spn.learncats(data, classcol=classcol)
    }
    else {
      summ <- spn.learncats(data, classcol=ncol(data))
    }
    if(summ$ncat[classcol] < 2) stop('class must be discrete')
    for(j in 1:nruns) {
        folds <- 1:n
        if(nfolds < n) {
            if(strat) {
                for(i in unique(data[, classcol])) {
                    nn <- sum(data[, classcol]==i)
                    ind <- rep(1:nfolds, ceiling(nn/nfolds))[1:nn]
                    folds[data[, classcol]==i] <- sample(ind, nn)
                }
            } else {
                ind <- rep(1:nfolds, ceiling(n/nfolds))[1:n]
                folds <- sample(ind, n)
            }
        }
        for(i in 1:nfolds) {
            cat(paste(Sys.time(),':: Run',j,'Fold',i,'training'),'\n')
            cur <- which(folds == i)
            if(ncores > 1) {
                cat(paste(Sys.time(),'::learn',"\n"),file=tfname,append=TRUE)
            }
            spn <- spn.learn(summ$data[-cur,],ncat=summ$ncat,maxv=summ$maxv,minv=summ$minv,verb=verb,classcol=classcol,thr=thr,height=height)
            cat(paste(Sys.time(),':: Run',j,'Fold',i,'testing'),'\n')
            res <- rbind(res, spn.predict(spn, data[cur, , drop=FALSE], classcol, run.rob=run.rob, ncores=ncores,tfname=tfname, use_memory=use_memory, len.missing=len.missing))
        }
    }
    if(!is.null(path)) write.csv(res, path)
    return(res)
}


spn.cv_acc <- function(data,
                    name,
                    classcol=ncol(data),
                    nfolds=nrow(data),
                    thr=0.01,
                    nruns=1,
                    verb=FALSE,
                    strat=TRUE,
                    height=1000000,
                    run.rob=TRUE,
                    ncores=detectCores()-1,
                    tfname=tempfile('outSPN.',tmpdir='.',fileext='.tmp'),
                    use_memory=FALSE,
                    path=NULL) {
    n <- nrow(data)
    nfolds <- min(n,nfolds)
    res_regular <- c()
    res_select <- c()
    res_xg <- c()
    summ <- spn.learncats(data, classcol=classcol)

    if(summ$ncat[classcol] < 2) stop('class must be discrete')
    for(j in 1:nruns) {
        folds <- 1:n
        if(nfolds < n) {
            if(strat) {
                for(i in unique(data[, classcol])) {
                    nn <- sum(data[, classcol]==i)
                    ind <- rep(1:nfolds, ceiling(nn/nfolds))[1:nn]
                    folds[data[, classcol]==i] <- sample(ind, nn)
                }
            } else {
                ind <- rep(1:nfolds, ceiling(n/nfolds))[1:n]
                folds <- sample(ind, n)
            }
        }
        for(i in 1:nfolds) {
            cat(paste(Sys.time(),':: Run',j,'Fold',i,'training'),'\n')
            cur <- which(folds == i)
            if(ncores > 1) {
                cat(paste(Sys.time(),'::learn',"\n"),file=tfname,append=TRUE)
            }

            ptm <- proc.time()
            spn_select <- spn.learn(summ$data[-cur,],ncat=summ$ncat,maxv=summ$maxv,minv=summ$minv,verb=verb,classcol=classcol,thr=thr,height=height)
            time <- proc.time() - ptm
            time <- as.numeric(time['sys.self'] + time['user.self'])

            cur_res_select <- spn.predict(spn_select, data[cur, , drop=FALSE], classcol, run.rob=run.rob, ncores=ncores,tfname=tfname, use_memory=FALSE)
            res_select <- rbind(res_select, cur_res_select)

            if(exists("spn_sel_info")) {
              spn_sel_info <- rbind(spn_sel_info, spn.get.info(spn_select, time, sum(cur_res_select[, 'count']), sum(cur_res_select[, 'time'])))
            } else {
              spn_sel_info <- spn.get.info(spn_select, time, sum(cur_res_select[, 'count']), sum(cur_res_select[, 'time']))
            }

            ptm <- proc.time()
            spn_regular <- spn.learn(summ$data[-cur,],ncat=summ$ncat,maxv=summ$maxv,minv=summ$minv,verb=verb,classcol=NULL,thr=thr,height=height)
            time <- proc.time() - ptm
            time <- as.numeric(time['sys.self'] + time['user.self'])

            cur_res_regular <- spn.predict(spn_regular, data[cur, , drop=FALSE], classcol, run.rob=run.rob, ncores=ncores,tfname=tfname, use_memory=FALSE)
            res_regular <- rbind(res_regular, cur_res_regular)

            if(exists("spn_reg_info")) {
              spn_reg_info <- rbind(spn_reg_info, spn.get.info(spn_regular, time, sum(cur_res_regular[, 'count']), sum(cur_res_regular[, 'time'])))
            } else {
              spn_reg_info <- spn.get.info(spn_regular, time, sum(cur_res_regular[, 'count']), sum(cur_res_regular[, 'time']))
            }

            # Train xgboost model
            xgb <- xgboost(data = data.matrix(data[-cur,-ncol(data)]),
             verbose=0,
             label = data[-cur, classcol]-1,
             eta = 0.1,
             max_depth = 15,
             nround=25,
             subsample = 0.5,
             colsample_bytree = 0.5,
             eval_metric = "mlogloss",
             objective = "multi:softmax",
             num_class = summ$ncat[classcol],
             nthread = 3
            )

            ##spn.print(spn)
            cat(paste(Sys.time(), ':: Run', j, 'Fold', i, 'testing'),'\n')
            validation <- cbind(data[cur, classcol], predict(xgb, data.matrix(data[cur, -ncol(data)]))+1)
            res_xg <- rbind(res_xg, validation)
        }
    }
    write.csv(spn_sel_info, paste(path, '_slInfo.csv', sep=""))
    write.csv(spn_reg_info, paste(path, '_rgInfo.csv', sep=""))
    write.csv(res_select, paste(path, '_sel.csv', sep=""))
    write.csv(res_regular, paste(path, '_reg.csv', sep=""))
    write.csv(res_xg, paste(path, '_xg.csv', sep=""))
}


spn.computeall <- function(spn, eps=0.2, use_memory=FALSE) {
    s <- 0
    pncat <- prod(spn$ncat)
    if(pncat > 1) {
        for(i in 1:pncat) {
            cfg=list()
            cfg$scope=1:length(spn$ncat)
            cfg$value=number2binary(i, length(spn$ncat))+1
            z=spn.value.int(spn, cfg, eps)
            vv <- spn.value(spn, cfg)
            v <- exp(vv[1])
            s <- s + v
            #print(paste(paste(cfg$value,collapse=' '),v))
            print(paste(paste(cfg$value,collapse=' '),exp(z[[1]][1]), v, exp(z[[2]][1])))
            vmi <- spn.value.min(spn, cfg, eps, use_memory=use_memory)$res
            vma <- spn.value.max(spn, cfg, eps, use_memory=use_memory)$res
            if(abs(exp(z[[1]][1])-exp(vmi)) > 1e-8) stop('error min')
            if(abs(exp(z[[2]][1])-exp(vma)) > 1e-8) stop('error max')
        }
        print(s)
    } else {
        print('only works for fully discrete data')
    }
}
