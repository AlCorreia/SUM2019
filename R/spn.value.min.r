## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.

## INPUT: spn (an SPN), evi (a config of variables for which we want to compute the marginal probability for that instantiation),
## OUTPUT: log probability that vars in scope are value
spn.value.min <- function(spn, evi, eps=0, eps.gauss=0, use_memory=FALSE) {
    ## quick check about the configuration
    if(sum(spn$ncat[evi$scope]*spn$ncat[evi$scope] < spn$ncat[evi$scope]*evi$value) > 0 || sum((spn$ncat[evi$scope] > 1) & (evi$value < 1)) > 0 || sum(spn$ncat[evi$scope]*(round(evi$value) != evi$value)) > 0) {
        stop('Invalid evidence configuration')
    }
    for(i in evi$scope) {
        if(spn$ncat[i] == 0) {
            evi$value[i] <- (evi$value[i] - spn$minv[i]) / (spn$maxv[i] - spn$minv[i])
        }
    }
    ## get the answer recursively
    return(spn.value.min.aux(spn$root, evi, eps, eps.gauss, use_memory))
}
## auxiliary function that does the job
spn.value.min.aux <- function(node, evi, eps, eps.gauss, use_memory) {
    evi <- reduce.cfg(node$scope, evi)
    if(length(evi$scope)==0) return(list("res"=log(1), "count"=1))
    if(node$type == 1) { #'leaf-indicator'
        ## for leaf nodes, return log(1) unless the var of this leaf appears in the evi config and is not compatiable with it
        pos <- which(node$scope == evi$scope)
        if(length(pos) > 0 && evi$value[pos] != node$value) {
            return(list("res"=log(0), "count"=1))
        }
        return(list("res"=log(1), "count"=1))
    }
    if(node$type == 2) { #'leaf-gaussian'
        ## for leaf nodes, return log(1) unless the var of this leaf appears in the config
        pos <- which(node$scope == evi$scope)
        if(length(pos) > 0) {
            eviv <- as.numeric(evi$value[pos])
            return(list("res"=dnorm(eviv,mean=(node$value[1]-eps.gauss/2),sd=node$value[2]), "count"=1))
        }
        return(list("res"=log(1), "count"=1))
    }
    # Check if evidence has been evaluated before
    if(use_memory) {
      evi.id <- paste(paste(evi$scope, evi$value, collapse=":"), paste("min", eps, eps.gauss), collapse="-")
    }
    if(use_memory && exists(evi.id, envir = node$memory)) {
      return(list("res"=env_get(node$memory, evi.id), "count"=0))
    } else {
      if(node$type == 3) { #'prod'
          ## for product nodes, return the sum of the result of the children (sum since they are logs)
          l <- length(node$children)
          res <- 0
          count <- 1
          for(nod in 1:l) {
              res_count <- spn.value.min.aux(node$children[[nod]], evi, eps, eps.gauss, use_memory)
              count <- count + res_count$count
              res <- res + res_count$res ## log product is sum of logs
              if(res == -Inf) break
          }
          if(use_memory) {
            env_poke(node$memory, evi.id, res) # Add result to hashtable
          }
          return(list("res"=res, "count"=count))
      }
      if(node$type == 4) { #'sum'
          ## for sum nodes, combine the results from the children with the appropriate weights
          l <- length(node$children)
          count <- 1
          vals.min <- rep_len(0,l)
          for(nod in 1:l) {
              res_count <- spn.value.min.aux(node$children[[nod]], evi, eps, eps.gauss, use_memory)
              count <- count + res_count$count
              vals.min[nod] <- res_count$res
          }
          vals.low <- node$weight * (1-eps)
          ##print(vals.low)
          total <- 1 - sum(vals.low)
          ##print(total)
          vals <- vals.low

          i.set <- order(vals.min)
          for(i in i.set) {
              if(total > eps) {
                  vals[i] <- vals[i] + eps
                  total <- total - eps
              } else {
                  vals[i] <- vals[i] + total
                  break
              }
          }
          res <- logsumexp(log(vals) + vals.min)
          if(use_memory) {
            env_poke(node$memory, evi.id, res) # Add result to hashtable
          }
          return(list("res"=res, "count"=count))
      }
   }
}
