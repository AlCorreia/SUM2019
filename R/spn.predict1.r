source('spn.value.r')
source('spn.utils.r')
source('spn.value.int.r')

spn.predict1 <- function(spn, data, classcol=ncol(data), verb=FALSE, run.rob=FALSE) {
    nclass <- spn$ncat[classcol]
    if(nclass < 2) stop('class must be discrete')
    res <- c()
    cfg <- list()
    len <- length(data)
    cfg$scope <- 1:len
    cfg$value <- data
    maxclass <- 1
    maxlogpr <- -Inf
    logprs <- c()
    for(j in 1:nclass) {
        cfg$value[classcol] <- j
        logpr <- signed.logunbuild(spn.value(spn, cfg))
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

    cfg$scope <- setdiff(1:len,classcol)
    cfg$value <- data[-classcol]
    rob <- 0
    if(run.rob) {
        f <- rep_len(0, nclass)
        f[maxclass] <- 1
        rob.l <- 0
        rob.u <- 1
        while(rob.l < rob.u - 0.005) {
            upper <- c()
            lower <- c()
            rob <- (rob.l + rob.u)/2
            ok <- TRUE
            for(j in setdiff(1:nclass,maxclass)) {
                f[j] <- -1
                z=spn.value.int(spn, cfg, eps=rob, query=classcol, f=f, eps.gauss=rob)
                f[j] <- 0
                if(signed.nonpositive(z$min)) {
                    ok <- FALSE
                    break
                }
            }
            if(ok) rob.l <- rob
            else rob.u <- rob
        }
    }
    slogprs <- logsumexp(logprs)
    if(verb && ncores<2) print(paste(paste(data,collapse=' '), data[classcol], maxclass, exp(maxlogpr-slogprs), rob))

    return(c(data[classcol], maxclass, exp(maxlogpr-slogprs), rob, logprs))
}
