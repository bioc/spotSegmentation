plotBlockImage <-
function(z, title = NULL, one = FALSE){

	nrowz <- nrow(z)
        ncolz <- ncol(z)

        if (!one) {
          par(mfrow = c(1,1), pty = "m")

          PIN <- par()$pin
          if (ncolz < nrowz) 
            par(pin=c(PIN[1]*(ncolz/nrowz), PIN[1]))
          else  
	    par(pin=c(PIN[1], PIN[1]*(nrowz/ncolz)))
        }

	
## reverse the y-axis of z, and transpose so that the image is the same as the raw image
if (any(z < 0)) stop("z has negative components")


image(t(z[nrowz:1,]), main=paste(title), xlab="", ylab="", col=gray((100:0)/100), axes=FALSE)

## set some components to NULL to prevent warnings
##      spar[c("cin", "cra", "csi", "cxy", "din", "gamma")] <- NULL

      if (!one) par(pin = PIN)
      invisible()
}
