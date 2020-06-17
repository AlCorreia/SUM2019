## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.

spn.evidence.int <- function(spn, evi, query, eps) {
    s <- spn.evidence.int.aux(spn$root, evi, query, eps)
    spn$root <- s[[3]]
    spn$ncat <- spn$ncat[spn$root$scope]
    return(spn)
}
spn.evidence.int.aux <- function(node, evi, query, eps) {
    if(length(node)==0) return(list(signed.build(1),FALSE,node))
    if(node$type == 'leaf-indicator') {
        return(list(spn.value.int.aux(node, evi, eps=eps),node$scope==query,node))
    }
    if(node$type == 'leaf-gaussian') {
        return(list(spn.value.int.aux(node, evi, eps=eps),FALSE,node)) ## log(1)
    }
    if(node$type == 'prod') {
        print(paste('prod',paste(node$scope,collapse=' ')))
        ## for product nodes, return the sum of the result of the children (sum since they are logs)
        l <- length(node$children)

        vals.max <- cbind(rep_len(0,l), rep_len(0,l))
        vals.min <- vals.max
        dec <- -1
        for(nod in 1:l) {
            res <- spn.evidence.int.aux(node$children[[nod]], evi, query, eps) ## log product is sum of logs
            vals.min[nod,] <- res[[1]]$min
            vals.max[nod,] <- res[[1]]$max
            node$children[[nod]] <- res[[3]]
            if(res[[2]]) dec <- nod
        }
        pr <- signed.pair(signed.prod(vals.min),signed.prod(vals.max))
        if(!signed.nonzero(pr[[2]]) || dec == -1) {
            print('prod superend')
            node$children <- list()
            node$scope <- c()
            return(list(pr,FALSE,node))
        } else {
            maxv <- signed.prod(vals.max[-dec,,drop=FALSE])
            minv <- signed.prod(vals.min[-dec,,drop=FALSE])
            pr <- signed.pair(minv,maxv)
            print(paste('prod end',TRUE,paste(node$children[[dec]]$scope,collapse=' ')))
            return(list(pr,TRUE,node$children[[dec]]))
        }
    }
    if(node$type == 'sum') {
        hasquery <- FALSE
        print('sum')
        ## for sum nodes, combine the results from the children with the appropriate weights
        l <- length(node$children)

        if(l != 2) stop('only implemented for sum nodes with exactly 2 children')

        vals.max <- cbind(rep_len(0,l), rep_len(0,l))
        vals.min <- vals.max
        vals.up <- signed.build((node$weight/sum(node$weight)) * (1-eps) + eps)
        vals.low <- signed.build((node$weight/sum(node$weight)) * (1-eps))
        node$up <- rep_len(0,l)
        node$low <- rep_len(0,l)
        yet <- 0
        scope <- c()
        for(nod in 1:l) {
            res <- spn.evidence.int.aux(node$children[[nod]], evi, query, eps)
            vals.max[nod,] <- res[[1]]$max
            vals.min[nod,] <- res[[1]]$min

            hasquery <- hasquery || res[[2]]
            if(!signed.nonzero(res[[1]]$max) || !res[[2]]) {
                node$weight[nod] <- 0
                node$children[[nod]] <- list()
            } else {
                node$children[[nod]] <- res[[3]]
                node$up[nod] <- signed.unbuild(res[[1]]$max)
                node$low[nod] <- signed.unbuild(res[[1]]$min)
                scope <- union(scope,res[[3]]$scope)
                yet <- yet+1
                nyet <- res[[3]]
                vyet <- res[[1]]
            }
        }
        if(yet == 1) {
            print('sum superend')
            return(list(vyet,hasquery,nyet))
        }

        val1min <- signed.sum(signed.prod(vals.min[1,,drop=FALSE], vals.up[1,,drop=FALSE]),signed.prod(vals.min[2,,drop=FALSE], vals.low[2,,drop=FALSE]))
        val2min <- signed.sum(signed.prod(vals.min[1,,drop=FALSE], vals.low[1,,drop=FALSE]),signed.prod(vals.min[2,,drop=FALSE], vals.up[2,,drop=FALSE]))
        val1max <- signed.sum(signed.prod(vals.max[1,,drop=FALSE], vals.up[1,,drop=FALSE]),signed.prod(vals.max[2,,drop=FALSE], vals.low[2,,drop=FALSE]))
        val2max <- signed.sum(signed.prod(vals.max[1,,drop=FALSE], vals.low[1,,drop=FALSE]),signed.prod(vals.max[2,,drop=FALSE], vals.up[2,,drop=FALSE]))
        ## since vals are logs, we need to combine them with log-sum-exp
        node$scope <- scope
        print(paste('sum end',hasquery,paste(node$scope,collapse=' ')))
	return(list(signed.pair(signed.min(val1min,val2min), signed.max(val1max,val2max)),hasquery,node))
    }


}
