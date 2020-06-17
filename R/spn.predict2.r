source('spn.value.r')
source('spn.utils.r')
source('spn.value.max.r')
source('spn.value.min.r')
library(rlang)
library(tictoc)


# This is the same as predict1 but optmised for class-selective SPNs.
spn.predict2 <- function(spn, data, classcol=length(data), verb=FALSE, run.rob=FALSE, use_memory=FALSE, len.missing=0) {
    nclass <- spn$ncat[classcol]
    if(nclass < 2) stop('class must be discrete')
    res <- c()
    cfg <- list()
    len <- length(data)
    cfg$scope <- sort(c(sample(setdiff(1:len,classcol),len-1-len.missing),classcol))
    cfg$value <- data[cfg$scope]
    maxclass <- 1
    maxlogpr <- -Inf
    logprs <- c()
    count <- 0
    ptm <- proc.time()
    for(j in 1:nclass) {
        cfg$value[which(classcol==cfg$scope)] <- j
        res_count <- spn.value.max(spn, cfg, use_memory=use_memory)
        logpr <- res_count$res
        count <- count + res_count$count
        logprs <- c(logprs,logpr)
        if(logpr > maxlogpr) {
            maxclass <- j
            maxlogpr <- logpr
        }
    }
    if(max(logprs)==-Inf) {
        logprs <- log(1/nclass + 0*(1:nclass))
        maxlogpr <- logprs[1]
        maxclass <- sample(1:nclass,1)
    }
    loop_count <- c(count)

    rob <- 0
    if(run.rob) {
        rob.l <- 0
        rob.u <- 1
        rob <- 0.5
        while(rob.l < rob.u - 0.005) {
            cfg$value[which(classcol==cfg$scope)] <- maxclass
            res_count <- spn.value.min(spn, cfg, eps=rob, eps.gauss=rob, use_memory=use_memory)
            v.min <- res_count$res
            count <- res_count$count
            for(j in order(logprs,decreasing=TRUE)) {
                if(j != maxclass) {
                    cfg$value[which(classcol==cfg$scope)] <- j
                    res_count <- spn.value.max(spn, cfg, eps=rob, eps.gauss=rob, use_memory=use_memory)
                    v.max <- res_count$res
                    count <- count + res_count$count
                    if(v.max >= v.min) {
                        rob.u <- rob
                        break
                    }
                }
            }
            if(rob < rob.u) rob.l <- rob
            rob <- (rob.l + rob.u)/2
            loop_count <- c(loop_count, count)
        }
    }
    slogprs <- logsumexp(logprs)
    time <- proc.time() - ptm
    timing <- as.numeric(time['sys.self'] + time['user.self'])
    if(verb && ncores<2) print(paste(paste(data,collapse=' '), data[classcol], maxclass, exp(maxlogpr-slogprs), rob))

    if(run.rob) {
      res <- c(as.numeric(data[classcol]),maxclass,exp(maxlogpr-slogprs),rob, logprs, sum(loop_count), timing, loop_count)
    } else {
      res <- c(as.numeric(data[classcol]),maxclass,exp(maxlogpr-slogprs),rob, logprs, loop_count, timing)
    }
    return(res)
}
