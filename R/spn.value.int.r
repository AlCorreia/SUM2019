## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.

## INPUT: spn (an SPN), evi (a config of variables for which we want to compute the marginal probability for that instantiation),
##        query (a variable to apply the function f), f (the function to be applied)
## OUTPUT: log probability that vars in scope are value
spn.value.int <- function(spn, evi, eps, query=NULL, f=NULL, eps.gauss=0) {
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
    return(spn.value.int.aux(spn$root, evi, query, f, eps, eps.gauss))
}
## auxiliary function that does the job
spn.value.int.aux <- function(node, evi, query=NULL, f=NULL, eps=0, eps.gauss=0) {
    if(node$type == 1) { #'leaf-indicator'
        if(!is.null(query) && node$scope == query) {
            ## if it is a queried variable, apply the function
            return(signed.pair(signed.build(f[node$value])))
        }
        ## for leaf nodes, return log(1) unless the var of this leaf appears in the evi config and is not compatiable with it
        pos <- which(node$scope == evi$scope)
        if(length(pos) > 0 && evi$value[pos] != node$value)
            return(signed.pair(signed.build(0))) ## log(0)
        return(signed.pair(signed.build(1))) ## log(1)
    }
    if(node$type == 2) { #'leaf-gaussian'
        ## for leaf nodes, return log(1) unless the var of this leaf appears in the config
        pos <- which(node$scope == evi$scope)
        if(length(pos) > 0) {
##            v1 <- min(1,dnorm(evi$value[pos],mean=(node$value[1]),sd=node$value[2]))
##            v2 <- v1
            eviv <- as.numeric(evi$value[pos])
            v1 <- min(1,dnorm(eviv,mean=(node$value[1]-eps.gauss/2),sd=node$value[2]))
            v2 <- min(1,dnorm(eviv,mean=(node$value[1]+eps.gauss/2),sd=node$value[2]))
            return(signed.pair(signed.build(min(v1,v2)),
                               signed.build(max(v1,v2)))) ## log(Normal(val,mean=node$value[1],sd=node$value[2]))
        }
        return(signed.pair(signed.build(1))) ## log(1)
    }
    if(node$type == 3) { #'prod'
        ## for product nodes, return the sum of the result of the children (sum since they are logs)
        l <- length(node$children)
        vals.max <- cbind(rep_len(0,l), rep_len(0,l))
        vals.min <- vals.max
        dec <- 0
        for(nod in 1:l) {
            res <- spn.value.int.aux(node$children[[nod]], evi, query, f, eps) ## log product is sum of logs
            vals.min[nod,] <- res$min
            vals.max[nod,] <- res$max
            if(!signed.nonzero(res$max)) return(signed.pair(res$max,res$max))
            if(!is.null(query) && query %in% node$children[[nod]]$scope) dec <- nod
        }
        ##        print('prod')
        ##        print(vals)
        ##        print(signed.prod(vals))
        if(dec == 0)
            return(signed.pair(signed.prod(vals.min),signed.prod(vals.max)))
        if(vals.min[dec,2] < 0)
            minv <- signed.prod(signed.prod(vals.max[-dec,,drop=FALSE]),vals.min[dec,,drop=FALSE])
        else minv <- signed.prod(vals.min)

        if(vals.max[dec,2] < 0)
            maxv <- signed.prod(signed.prod(vals.min[-dec,,drop=FALSE]),vals.max[dec,,drop=FALSE])
        else maxv <- signed.prod(vals.max)

        return(signed.pair(minv,maxv))
    }
    if(node$type == 4) { #'sum'
        ## for sum nodes, combine the results from the children with the appropriate weights
        l <- length(node$children)

        vals.max <- cbind(rep_len(0,l), rep_len(0,l))
        vals.min <- vals.max
        vals.up <- signed.build(node$weight/sum(node$weight) * (1-eps) + eps)
        vals.low <- signed.build(node$weight/sum(node$weight) * (1-eps))

        nonz <- c()
        for(nod in 1:l) {
            res <- spn.value.int.aux(node$children[[nod]], evi, query, f, eps)
            vals.max[nod,] <- res$max
            vals.min[nod,] <- res$min
            if(signed.nonzero(res$max) || signed.nonzero(res$min)) nonz <- c(nonz, nod)
        }
        if(length(nonz) == 0) return(signed.pair(signed.build(0),signed.build(0)))
        if(length(nonz) > 2) stop('only implemented for sum nodes with at most 2 nonzero children')
        if(length(nonz) == 1)
            return(signed.pair(signed.prod(vals.min[nonz,,drop=FALSE], vals.low[nonz,,drop=FALSE]), signed.prod(vals.max[nonz,,drop=FALSE], vals.up[nonz,,drop=FALSE])))

        val1min <- signed.sum(signed.prod(vals.min[nonz[1],,drop=FALSE], vals.up[nonz[1],,drop=FALSE]),signed.prod(vals.min[nonz[2],,drop=FALSE], vals.low[nonz[2],,drop=FALSE]))
        val2min <- signed.sum(signed.prod(vals.min[nonz[1],,drop=FALSE], vals.low[nonz[1],,drop=FALSE]),signed.prod(vals.min[nonz[2],,drop=FALSE], vals.up[nonz[2],,drop=FALSE]))
        val1max <- signed.sum(signed.prod(vals.max[nonz[1],,drop=FALSE], vals.up[nonz[1],,drop=FALSE]),signed.prod(vals.max[nonz[2],,drop=FALSE], vals.low[nonz[2],,drop=FALSE]))
        val2max <- signed.sum(signed.prod(vals.max[nonz[1],,drop=FALSE], vals.low[nonz[1],,drop=FALSE]),signed.prod(vals.max[nonz[2],,drop=FALSE], vals.up[nonz[2],,drop=FALSE]))
        return(signed.pair(signed.min(val1min,val2min), signed.max(val1max,val2max)))
    }
}
