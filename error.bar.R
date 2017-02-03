error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

## Usage
y <- rnorm(500, mean=1)
y <- matrix(y,100,5)
y.means <- apply(y,2,mean)
y.sd <- apply(y,2,sd)
barx <- barplot(y.means, names.arg=1:5,ylim=c(0,1.5), col="blue", axis.lty=1, xlab="Replicates", ylab="Value (arbitrary units)")
error.bar(barx,y.means, 1.96*y.sd/10)
