## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.

spn.evidence <- function(spn, evi, query) {
    s <- spn.evidence.aux(spn$root, evi, query)
    print(signed.unbuild(s[[1]]))
    spn$root <- s[[3]]
    spn$ncat <- spn$ncat[spn$root$scope]
    return(spn)
}
spn.evidence.aux <- function(node, evi, query) {
    if(length(node)==0) return(list(signed.build(1),FALSE,node))
    if(node$type == 'leaf-indicator') {
        return(list(spn.value.aux(node, evi),node$scope==query,node))
    }
    if(node$type == 'leaf-gaussian') {
        return(list(spn.value.aux(node, evi),FALSE,node)) ## log(1)
    }
    hasquery <- FALSE
    if(node$type == 'prod') {
        print('prod')
        ## for product nodes, return the sum of the result of the children (sum since they are logs)
        l <- length(node$children)
        vals <- cbind(rep_len(0,l), rep_len(0,l))
        for(nod in 1:l) {
            res <- spn.evidence.aux(node$children[[nod]], evi, query)
            vals[nod,] <- res[[1]] ## log product is sum of logs
            hasquery <- hasquery || res[[2]]
            if(!res[[2]]) node$children[[nod]] <- list()
            else node$children[[nod]] <- res[[3]]
        }
        pr <- signed.prod(vals)
        if(!signed.nonzero(pr)) {
            print('prod superend')
            return(list(pr,FALSE,node))
        }
        i <- 0
        scope <- c()
        children <- list()
        for(nod in node$children) {
            if(length(nod)>0) {
                i <- i+1
                children[[i]] <- nod
                scope <- union(scope,nod$scope)
            }
        }
        print('prod end')
        if(i == 1) return(list(pr,hasquery,children[[1]]))
        node$children <- children
        node$scope <- scope
        return(list(pr,hasquery,node))
    }
    if(node$type == 'sum') {
        print('sum')
        ## for sum nodes, combine the results from the children with the appropriate weights
        l <- length(node$children)
        vals <- signed.build(node$weight/sum(node$weight))
        yet <- 0
        scope <- c()
        for(nod in 1:l) {
            res <- spn.evidence.aux(node$children[[nod]], evi, query)
            hasquery <- hasquery || res[[2]]
            if(!signed.nonzero(res[[1]])) {
                node$weight[nod] <- 0
                node$children[[nod]] <- list()
            } else {
                node$children[[nod]] <- res[[3]]
                node$weight[nod] <- node$weight[nod] * signed.unbuild(res[[1]])
                scope <- union(scope,res[[3]]$scope)
                yet <- yet+1
                nyet <- res[[3]]
                vyet <- nod
            }
            vals[nod,] <- signed.prod(vals[nod,,drop=FALSE], res[[1]])
	}
        if(yet == 1) {
            print('sum superend')
            return(list(vals[vyet,,drop=FALSE],hasquery,nyet))
        }
        ## since vals are logs, we need to combine them with log-sum-exp
        print('sum end')
        node$scope <- scope
	return(list(signed.sum(vals),hasquery,node))
    }
}
