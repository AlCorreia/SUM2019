## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.
setRefClass("node",
            fields=list(children="list", scope="vector", weight="vector", n="numeric", type="numeric", value="vector", len="numeric", size="numeric",id="numeric", memory="environment")
            )
setRefClass("spn",
            fields=list(root="node",ncat="vector",maxv="vector",minv="vector")
            )

## INPUT: data (the dataset, columns are variables, rows are iid samples), ncat is an array with the number of categories per variable (each variable var must have values in 1:ncat[var]), thr is a threshold for chisq.test, classcol is the column of the class variable, verbose is a Boolean
## OUTPUT: an SPN inspired by the fitting idea of Gens&Domingos 2013
spn.learn <- function(data, ncat=NULL, thr=0.001, classcol=NULL, verb=FALSE, minv=NULL, maxv=NULL, height=1000000, class.unif=FALSE) {
    ## clean up of data and computation of number of categories per variable (unless already given)
    if(is.null(ncat) || is.null(maxv) || is.null(minv)) {
        l <- spn.learncats(data, verb, classcol=classcol)
        ncat <- l$ncat
        data <- l$data
        maxv <- l$maxv
        minv <- l$minv
    }
    ## print(data)
    ## in this implementation, the root node is a sum node related to the class variable (kind of a more discriminative approach), if that is given
    scope <- 1:ncol(data)
    root <- spn.learn.aux(data, ncat, scope=scope, thr=thr, nclusters=2, verb=verb, height=height, classcol=classcol, classvalues=data[,classcol], last.prod=TRUE, class.unif=FALSE)
    return(new("spn",root=root,ncat=ncat,maxv=maxv,minv=minv))
}
spn.learncats <- function(data.ori,verb=FALSE,classcol=ncol(data)) {
    data <- matrix(0, nrow=nrow(data.ori), ncol=ncol(data.ori))
    ncat <- rep.int(0, ncol(data))
    maxv <- ncat
    minv <- ncat
    if(is.null(classcol)) classcol <- ncol(data)
    for(i in 1:ncol(data)) {
        if(i != classcol && (min(data.ori[,i]) < 1 || sum(data.ori[,i] != round(data.ori[,i])) > 0 || length(unique(data.ori[,i])) > min(30,length(data.ori[,i]/3)))) {
            warning(paste('column',i,'considered Gaussian\n'),immediate.=TRUE)
            ## make them be in [0,1], good for comparison metrics during clustering
            ##            m <- min(data.ori[,i])
            ##            data[,i] <- (data.ori[,i] - m) / (max(data.ori[,i]) - m)
            ## make them standardized as a normal(0,1), can still be good for clustering
            ## data[,i] <- (data.ori[,i] - mean(data.ori[,i])) / sd(data.ori[,i])
            maxv[i] <- max(data.ori[,i])
            minv[i] <- min(data.ori[,i])
            if(maxv[i]==minv[i]) stop('cannot have constant continuous variable in the data')
            data[,i] <- (data.ori[,i] - minv[i])/(maxv[i] - minv[i])
        } else {
            data[,i] <- as.integer(data.ori[,i])
            ncat[i] <- max(data[,i])
        }
    }

    return(list(data=data,ncat=ncat,maxv=maxv,minv=minv))
}
## learning is done recursively, splitting the data horizontally by clustering and vertically by "independence" tests
spn.learn.aux <- function(data, ncat, scope, thr, nclusters, verb, height, classcol, last.prod=FALSE, classvalues, class.unif=FALSE) {
    n <- length(scope)
    m <- nrow(data)
    if(verb) cat(paste(Sys.time(),":: New node with",m,"points and",n,"vars, scope",paste(scope,collapse=' '),"\n"))
    if (n > 1) {
        if((height <= 0) || (!last.prod)) { ## just to speed up, since never a product node will have a product node as child (easy to see that)
            ## let us build an undirected graph where two nodes (variables) are connected if they are dependent
            deplist <- list()
            depfunc <- function(i) { if(deplist[[i]] == i) return(i) else return(depfunc(deplist[[i]])); }
            for(i in 1:n) deplist[[i]] <- i
            if(height > 0) {
                if(verb) cat(paste(Sys.time(),":: Finding components\n"))
                for(i in 1:(n-1)) {
                    for(j in (i+1):n) {
                        fatheri <- depfunc(i)
                        deplist[[i]] <- fatheri
                        fatherj <- depfunc(j)
                        deplist[[j]] <- fatherj
                        if(fatheri != fatherj) {
                            v <- 1
                            unii <- length(unique(data[,scope[i]]))
                            unij <- length(unique(data[,scope[j]]))
                            if(unii > 1 && unij > 1) {
                                if(m > 4 && ncat[scope[i]] == 0 && ncat[scope[j]] == 0) {
                                    ## both continuous
                                    ##                                print(cor.test(data[,scope[i]],data[,scope[j]],method='kendall'))
                                    v <- suppressWarnings(cor.test(data[,scope[i]],data[,scope[j]],method='kendall'))$p.value
                                }
                                if(m > 4*unij && ncat[scope[i]] == 0 && ncat[scope[j]] > 1) {
                                    ## i continuous, j discrete
                                    ##                                print(kruskal.test(data[,scope[i]],data[,scope[j]]))
                                    v <- suppressWarnings(kruskal.test(data[,scope[i]],data[,scope[j]]))$p.value
                                }
                                if(m > 4*unii && ncat[scope[i]] > 1 && ncat[scope[j]] == 0) {
                                    ## i discrete, j continuous
                                    ##                                print(kruskal.test(data[,scope[j]],data[,scope[i]]))
                                    v <- suppressWarnings(kruskal.test(data[,scope[j]],data[,scope[i]]))$p.value
                                }
                                if(m > unii*unij*2 && ncat[scope[i]] > 1 && ncat[scope[j]] > 1) {
                                    ## both discrete
                                    ##                                print(chisq.test(data[,scope[i]],data[,scope[j]]))
                                    v <- suppressWarnings(chisq.test(data[,scope[i]],data[,scope[j]]))$p.value
                                }
                                if(v < thr) {
                                    deplist[[fatherj]] <- fatheri
                                }
                            }
                        }
                    }
                }
            }
            clu <- list(no=0,membership=1:n)
            for(i in 1:n) {
                clu$membership[i] <- depfunc(i)
            }
            clu$unique <- unique(clu$membership)
            clu$no <- length(clu$unique)
            for(i in 1:n) {
                clu$membership[i] <- which(clu$membership[i] == clu$unique)
            }
            ## g <- graph_from_adjacency_matrix(dep,mode='undirected')
            ## and from such graph, take the components to form the children of the product node, as long as there are more than one component (single component means that we cannot split vertically at this moment, and in that case we move on to split horizontally
            ## clu <- clusters(g)
            if(clu$no > 1) {
                if(verb) cat(paste(Sys.time(),":: Found independency (",clu$no,"groups). Separating independent sets\n"))
                for(i in 1:clu$no) {
                ##     if(TRUE) {
                ##         ## GINI VERSION PER COLUMN
                ##         g <- 0
                ##         for(ii in scope[which(clu$membership == i)]) {
                ##             g <- c(g,min(gini(data[,ii]),1,na.rm=T))
                ##         }
                ##     } else {
                ##         ## PVALUE PER COLUMN WRT CLASS VARIABLE - THIS IS NOT USEFUL RIGHT NOW, AS THE CLASS SPLITS THE DATA IN THE BEGINNING!!
                ##         g <- 0
                ##         for(ii in scope[which(clu$membership == i)]) {
                ##             if(length(unique(classvalues))==1) g <- c(g, 0)
                ##             else {
                ##                 if(ncat[scope[ii]] == 0) g <- c(g, kruskal.test(data[,ii],classvalues)$p.value)
                ##                 else g <- c(g, chisq.test(data[,ii],classvalues)$p.value)
                ##             }
                ##         }
                ##     }
                ##     clu$weight <- c(clu$weight, 1/(1+mean(g)))
                    clu$weight <- c(clu$weight, 1)
                }
                # prodnode <- new("node",children=list(), scope=scope, n=m, type=3, len=clu$no, weight=rep_len(1,clu$no), size=1,id=round(10000000*runif(1))) #'prod'
                prodnode <- new("node",children=list(), scope=scope, n=m, type=3, len=clu$no, weight=as.vector(clu$weight), size=1,id=round(10000000*runif(1))) #'prod'
                sizes <- rep_len(0,clu$no)
                children <- list()
                for(i in 1:clu$no) {
                    children[[i]] <- spn.learn.aux(data, ncat, scope=scope[which(clu$membership == i)], thr=thr, nclusters=nclusters, verb=verb, last.prod=TRUE, height=height-1, classcol=classcol, classvalues=classvalues)
                    prodnode$size <- prodnode$size + children[[i]]$size
                    sizes[i] <- children[[i]]$size
                }
                o <- order(sizes)
                for(i in 1:clu$no) {
                    prodnode$children[[i]] <- children[[o[i]]]
                    prodnode$weight[i] <- clu$weight[o[i]]
                }
                return(prodnode)
            }
            if(verb) cat(paste(Sys.time()," :: No independencies found\n"))
        }
    } else {
        ## single variable in the scope
        ## has it a single value? If so, then place an indicator function (as leaf node)
        if(ncat[scope] > 1) {
            if(verb) cat(paste(Sys.time(),":: Creating new indicator leaf\n"))
            ## if (length(unique(data[,scope]))==1) {
            ##     return(new("node",scope=scope, value=data[1,scope], type=1,children=list(), len=0, size=1,id=round(10000000*runif(1)))) #'leaf-indicator'
            ## }
            ##### if multiple values for the variable are still present in the data, then use a sum node with indicator functions as children
            ncategory <- ncat[scope]
            sumnode <- new("node",children=list(), scope=scope, weight=rep_len(0,ncategory), n=m, type=4, len=ncategory, size=ncategory+1,id=round(10000000*runif(1))) #'sum'

            for(i in 1:ncategory) {
                members <- which(data[,scope]==i)
                sumnode$children[[i]] <- new("node",scope=scope, value=i,  type=1,children=list(), len=0, size=1,id=round(10000000*runif(1))) #'leaf-indicator'
                sumnode$weight[i] <- length(members)+0.01  ##+0.5 is regularisation
            }
            ss <- sum(sumnode$weight)
            for(i in 1:ncategory) {
                sumnode$weight[i] <- sumnode$weight[i]/ss
            }
            return(sumnode)
        } else {
            if(verb) cat(paste(Sys.time(),":: Creating new Gaussian leaf\n"))
            return(new("node",scope=scope, value=c(mean(data[,scope]),max(1e-4,sd(data[,scope]),na.rm=TRUE)), type=2,children=list(), len=0, size=1,id=round(10000000*runif(1)))) #'leaf-gaussian'
        }
    }

    ## we were not able to cut vertically, so we run clustering of the data (each row is a point), and then we use the result of clustering to create the children of a sum node
    sumnode <- new("node",children=list(), scope=scope, weight=vector(), n=m, len=0, type=4, size=1,id=round(10000000*runif(1))) #'sum'
    do.unif <- FALSE
    ##do.empty <- FALSE
    nclusters.t <- nclusters
    if(!is.null(classcol) && !(classcol %in% scope)) classcol <- NULL
    if(!is.null(classcol) && ncat[classcol] > 1) {
        if(verb) cat(paste(Sys.time(),":: Creating sum node with class as cluster\n"))
        clu.ind <- data[,classcol]
        do.unif <- class.unif
        classcol <- NULL
        nclusters.t <- length(unique(clu.ind))
        ##do.empty <- TRUE
    } else {
        if(verb) cat(paste(Sys.time(),":: Clara clustering with",nrow(data),"rows and",nclusters,"clusters\n"))
        if(nclusters >= nrow(data)) clu.ind <- 1:nrow(data)
        else {
            ## clara is very efficient as clustering, but not necessarily very accurate - can you live with that? :-)
            ## clu.ind <- clara(data[,scope,drop=FALSE],nclusters,samples=20,pamLike=TRUE,metric='manhattan')$clustering
            clu.ind <- clara2(data[,scope,drop=FALSE],nclusters,ncat=ncat[scope],samples=20,pamLike=TRUE,metric='manhattan')$clustering
        }
    }
    j <- 0
    ## after clusters are found, build the children using that partition of the data
    sizes <- c()
    weights <- c()
    children <- list()
    for(i in 1:nclusters.t) {
        members <- which(clu.ind == i)
        if(length(members) > 0) {
            j <- j + 1
            children[[j]] <- spn.learn.aux(data[members,,drop=FALSE], ncat, scope=scope, thr=thr, nclusters=nclusters, verb=verb, last.prod=FALSE, height=height-1, classcol=classcol, classvalues=classvalues[members])
            if(do.unif) weights <- c(weights,1)
            else
                weights <- c(weights, length(members))
            sizes <- c(sizes, children[[j]]$size)
            sumnode$size <- sumnode$size + children[[j]]$size
        }
        ## else {
        ##     if(do.empty) {
        ##         j <- j + 1
        ##         children[[j]] <- new("node",scope=classcol, value=i, type=1,children=list(), len=0, size=1,id=round(10000000*runif(1)))
        ##         if(do.unif) weights <- c(weights,1)
        ##         else
        ##             weights <- c(weights, 0.5)
        ##         sizes <- c(sizes, 1)
        ##         sumnode$size <- sumnode$size + 1
        ##     }
        ## }
    }
    o <- order(sizes)
    sumnode$weight <- as.vector(weights)
    for(i in 1:j) {
        sumnode$children[[i]] <- children[[o[i]]]
        sumnode$weight[i] <- weights[o[i]]/sum(weights)
    }
    return(sumnode)
}
