"plot.spotseg" <-
function(x, ...) {

col <- NULL
title <- NULL
file <- NULL

plotSpotImage <-
function(z, col = NULL, title = NULL, one = FALSE){

        if (is.null(col)) col <- c("lightyellow", "gray", "black", "yellow",
                                   "blue", "red")

        if (!one) par(mfrow = c(1,1), pty = "m")
       
	nrowz <- nrow(z)
        ncolz <- ncol(z)

        if (!one) {
          PIN <- par()$pin
          if (ncolz < nrowz) 
            par(pin=c(PIN[1]*(ncolz/nrowz), PIN[1]))
          else  
	    par(pin=c(PIN[1], PIN[1]*(nrowz/ncolz)))
        }
	
## reverse the y-axis of z, and transpose so that the image is the same as the raw image
#z[z < 0] <- 0

z[z < 0] <- 3 + abs(z[z < 0])
u <- sort(unique(as.vector(z)))
if (length(col) < max(u)) stop("not enough colors")
L <- length(u)
breaks <- rep(0, L+1)
breaks[1] <- u[1] - .5
breaks[L+1] <- u[L] + .5
if (L > 1) {
  for (i in 2:L) breaks[i] <- (u[i-1] + u[i])/2
}
col <- col[u]

##col <- col[z]
##col <- array(col,dim(z))
##col <- t(col[nrowz:1,])
if (FALSE) {
# show a pixel
z[1,1] <- min(breaks) - .5
col <- c("blue", col)
breaks <- c(min(breaks) - 1, breaks)
}
image(t(z[nrowz:1,]), main=paste(title), xlab="", ylab="", col=col, 
      breaks = breaks, axes=FALSE)

## set some components to NULL to prevent warnings
##      spar[c("cin", "cra", "csi", "cxy", "din", "gamma")] <- NULL
        if (!one) par(pin = PIN)

      invisible()
}

   if (!is.null(file))
     postscript(file, horizontal=FALSE)

   plotSpotImage(x,col,title)    

   invisible()
}

