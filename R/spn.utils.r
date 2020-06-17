## Copyright (c) 2019 C. P. de Campos (cassiopc@acm.org). All rights reserved.

logsumexp <- function(x) {
    p <- which.max(x)
    if(x[p] == -Inf) return(-Inf)
    return(log1p(sum(exp(x[-p] - x[p]))) + x[p])
}

signed.sum <- function(x,y=NULL) {
    if(!is.null(y))
        x <- rbind(x,y)
    p <- signed.which.max(x)
    if(x[p,1] == -Inf) return(x[p,])
    if(is.vector(x)) return(x)
    r <- sum(exp(x[-p,1] - x[p,1])*(x[-p,2]*x[p,2]))
    if (r < -1.0) return(cbind(x[p,1] + log(-1.0-r), -x[p,2]))
    return(cbind(x[p,1] + log1p(r), x[p,2]))
}
signed.which.max <- function(x) {
    if(is.vector(x)) return(1)
    return(which.max(x[,1]))
}
signed.prod <- function(x,y=NULL) {
    if(is.null(y)) {
        if(is.vector(x)) return(x)
        return(cbind(sum(x[,1]),prod(x[,2])))
    }
    x <- as.vector(x)
    y <- as.vector(y)
    return(cbind(x[1]+y[1],x[2]*y[2]))
}
signed.exp <- function(x,y) { ## y is not signed.number, is a normal number
    x <- as.vector(x)
    return(cbind(x[1]*y,x[2]))
}
signed.build <- function(x) {
    return(cbind(log(abs(x)), sign(x)))
}
signed.max <- function(x,y) {
    x <- as.vector(x)
    y <- as.vector(y)
    if(x[2] < y[2]) return(y)
    if(y[2] < x[2]) return(x)
    if(x[1] < y[1]) return(y)
    return(x)
}
signed.min <- function(x,y) {
    x <- as.vector(x)
    y <- as.vector(y)
    if(x[2] < y[2]) return(x)
    if(y[2] < x[2]) return(y)
    if(x[1] < y[1]) return(x)
    return(y)
}
signed.nonpositive <- function(x) {
    x <- as.vector(x)
    return(x[2] <= 0)
}
signed.nonnegative <- function(x) {
    x <- as.vector(x)
    return(x[2] >= 0)
}
signed.nonzero <- function(x) {
    x <- as.vector(x)
    return(x[1] > -Inf && x[2] != 0)
}
signed.unbuild <- function(x) {
    x <- as.vector(x)
    return(x[2]*exp(x[1]))
}
signed.logunbuild <- function(x) {
    x <- as.vector(x)
    if(x[2] < 0) stop('log of negative cannot unbuild')
    return(x[1])
}
signed.pair <- function(x,y=NULL) {
    if(is.null(y))
        return(list(min=x,max=x))
    else return(list(min=x,max=y))
}

number2binary <- function(number, noBits) {
    binary_vector = rev(as.numeric(intToBits(number)))
    if(missing(noBits)) {
        return(binary_vector)
    } else {
        binary_vector[-(1:(length(binary_vector) - noBits))]
    }
}

## some tests
if(FALSE) {
    x=signed.build(c(-1,-2,0,1,2,3))
    signed.unbuild(x)
    y=c(-1,-2,1,2,3)
    x=signed.build(y)
    prod(y)
    exp(signed.prod(x))
}


getDummies <- function(data) {
    for(var in colnames(data)){
        data[[var]] <- as.integer(as.factor(data[[var]]))
        uni <- sort(unique(data[[var]]))
        for(i in 1:length(data[[var]])) {
            data[[var]][i] <- which(data[[var]][i]==uni)
        }
    }
    return(data)
}
