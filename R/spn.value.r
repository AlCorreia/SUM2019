## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.

## INPUT: spn (an SPN), evi (a config of variables for which we want to compute the marginal probability for that instantiation),
##        query (a variable to apply the function f), f (the function to be applied)
## OUTPUT: log probability that vars in scope are value
spn.value <- function(spn, evi, query=NULL, f=NULL) {
    ## quick check about the configuration
    if(sum(spn$ncat[evi$scope]*spn$ncat[evi$scope] < spn$ncat[evi$scope]*evi$value) > 0 || sum((spn$ncat[evi$scope] > 1) & (evi$value < 1)) > 0 || sum(spn$ncat[evi$scope]*(round(evi$value) != evi$value)) > 0) {
        stop('Invalid evidence configuration')
    }
    if(!is.null(query) || !is.null(f)) {
        if(spn$ncat[query] < 2 || length(f) != spn$ncat[query])
            stop('Invalid query variable or function')
    }
    for(i in evi$scope) {
        if(spn$ncat[i] == 0) {
            evi$value[i] <- (evi$value[i] - spn$minv[i]) / (spn$maxv[i] - spn$minv[i])
        }
    }
    ## get the answer recursively
    return(spn.value.aux(spn$root, evi, query, f))
}
## auxiliary function that does the job
spn.value.aux <- function(node, evi, query=NULL, f=NULL) {
    if(node$type == 1) { #'leaf-indicator'
        if(!is.null(query) && node$scope == query) {
            ## if it is a queried variable, apply the function
            return(signed.build(f[node$value]))
        }
        ## for leaf nodes, return log(1) unless the var of this leaf appears in the evi config and is not compatiable with it
        pos <- which(node$scope == evi$scope)
        if(length(pos) > 0 && evi$value[pos] != node$value)
            return(signed.build(0)) ## log(0)
        return(signed.build(1)) ## log(1)
    }
    if(node$type == 2) { #'leaf-gaussian'
        ## for leaf nodes, return log(1) unless the var of this leaf appears in the config
        pos <- which(node$scope == evi$scope)
        if(length(pos) > 0) {
            ## print(list(val=evi$value[pos],mean=node$value[1],sd=node$value[2]))
            eviv <- as.numeric(evi$value[pos])
            return(signed.build(min(1,dnorm(eviv,mean=node$value[1],sd=node$value[2])))) ## log(Normal(val,mean=node$value[1],sd=node$value[2]))
        }
        return(signed.build(1)) ## log(1)
    }
    if(node$type == 3) { #'prod'
        ## for product nodes, return the sum of the result of the children (sum since they are logs)
        l <- length(node$children)
        vals <- signed.build(node$weight)
        for(nod in 1:l) {
            res <- spn.value.aux(node$children[[nod]], evi, query, f) ## log product is sum of logs
            if(!signed.nonzero(res)) return(res)
            vals[nod,] <- signed.exp(res, node$weight[nod])
        }
        return(signed.prod(vals))
    }
    if(node$type == 4) { #'sum'
        ## for sum nodes, combine the results from the children with the appropriate weights
        l <- length(node$children)
        vals <- signed.build(node$weight/sum(node$weight))
        for(nod in 1:l) {
            if(node$weight[nod] != 0) {
                res <- spn.value.aux(node$children[[nod]], evi, query, f)
                vals[nod,] <- signed.prod(vals[nod,,drop=FALSE], res)
            }
	}
        ## weird case: only one child, then return its value
        if(length(node$children) == 1) return(vals)
        ## since vals are logs, we need to combine them with log-sum-exp
	return(signed.sum(vals))
    }
}
