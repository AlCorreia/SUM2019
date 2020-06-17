## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.

## INPUT: spn (an SPN)
spn.get.info <- function(spn, ltime, count, itime) {
  info <- list("SPN nodes"=spn.nnodes(spn),
               "sum"=spn.nsumnodes(spn),
               "prod"=spn.nprodnodes(spn),
               "leaf_indicator"=spn.nleafindicatornodes(spn),
               "leaf_gaussian"=spn.nleafgaussiannodes(spn),
               "params"=spn.nparams(spn),
               "height"=spn.height(spn),
               "learning time"=ltime,
               "inference time"=itime,
               "count"=count)
  return(info)
}

spn.print <- function(spn) {
    ## show some general stats about the spn
    cat(paste("SPN nodes=",spn.nnodes(spn),"sum=",spn.nsumnodes(spn),"prod=",spn.nprodnodes(spn),"leaf_indicator=",spn.nleafindicatornodes(spn),"leaf_gaussian=",spn.nleafgaussiannodes(spn),"params=",spn.nparams(spn),"height=",spn.height(spn),"\n"))
    cat("states=",paste(spn$ncat),",maxv=",paste(spn$maxv),",minv=",paste(spn$minv),'\n')
    ## and print the spn in a nice way, recursively
    return(spn.print.aux(spn$root, ""))
}

## recursively function to print an spn; spaces will show how deep the nodes are when printing
spn.print.aux <- function(node, spaces='') {
    if(length(node) == 0) return()
    if(node$type == 2) { #'leaf-gaussian'
        cat(paste(spaces,"leaf-gaussian var",node$scope,"=N(",node$value[1],",",node$value[2],")\n"))
    }
    if(node$type == 1) { #'leaf-indicator'
        cat(paste(spaces,"leaf-indicator var",node$scope,"=",node$value,"\n"))
    }
    if(node$type == 3) { #'prod'
        cat(paste(spaces,"prod scope",paste(node$scope,collapse=" "),"with n=",node$n," w=",paste(node$weight,collapse=" "),"\n"))
        for(nod in node$children) {
            spn.print.aux(nod, paste(spaces," "))
        }
    }
    if(node$type == 4) { #'sum'
        cat(paste(spaces,"sum scope",paste(node$scope,collapse=" "),"with n=",node$n," w=",paste(node$weight,collapse=" "),"\n"))
        if(is.element("up", names(node)))
            cat(paste(spaces,"      up=",paste(node$up,collapse=" ")," low=",paste(node$low,collapse=" "),"\n"))
        for(nod in node$children) {
            spn.print.aux(nod, paste(spaces," "))
	}
    }
}

spn.height <- function(spn) {
    return(spn.height.aux(spn$root))
}
spn.height.aux <- function(node) {
    if(length(node) == 0) return(0)
    if(node$type == 1 || node$type == 2) return(0) #'leaf-'
    h <- 0
    for(nod in node$children) {
        h <- max(h, spn.height.aux(nod))
    }
    return(h+1)
}
spn.nnodes <- function(spn) {
    return(spn.nnodes.aux(spn$root,type=0)) #'all'
}
spn.nsumnodes <- function(spn) {
    return(spn.nnodes.aux(spn$root,type=4)) #'sum'
}
spn.nprodnodes <- function(spn) {
    return(spn.nnodes.aux(spn$root,type=3)) #'prod'
}
spn.nleafnodes <- function(spn) {
    return(spn.nleafgaussiannodes(spn) + spn.nleafindicatornodes(spn))
}
spn.nleafgaussiannodes <- function(spn) {
    return(spn.nnodes.aux(spn$root,type=2)) #'leaf-gaussian'
}
spn.nleafindicatornodes <- function(spn) {
    return(spn.nnodes.aux(spn$root,type=1)) #'leaf-indicator'
}
spn.nnodes.aux <- function(node,type) {
    if(length(node) == 0) return(0)
    h <- 1*(type == 0 || type == node$type) #'all'
    for(nod in node$children) {
        h <- h + spn.nnodes.aux(nod,type)
    }
    return(h)
}

spn.nparams <- function(spn) {
    return(spn.nparams.aux(spn$root))
}
spn.nparams.aux <- function(node) {
    if(length(node) == 0) return(0)
    if(node$type == 2) return(2) #'leaf-gaussian'
    if(node$type == 1) return(0) #'leaf-indicator'
    if(node$type == 3) #'sum'
        h <- length(node$weight)
    else
        h <- 0
    for(nod in node$children) {
        h <- h + spn.nparams.aux(nod)
    }
    return(h)
}
