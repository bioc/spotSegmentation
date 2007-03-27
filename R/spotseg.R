"spotseg" <-
function (chan1, chan2, rowcut, colcut, R = NULL, C = NULL, threshold = 100, 
    hc = FALSE, show = FALSE) 
{
    spotseg1 <- function(spot, i, j, threshold = 100, ccl = TRUE, 
        hc = FALSE, show = FALSE) {
        dis <- function(spot, k) {
            L <- spot == k
            n <- sum(as.numeric(L))
            IJ <- matrix(NA, n, 2)
            IJ[, 1] <- ((row(spot) - (nrow(spot) + 1)/2))[L]
            IJ[, 2] <- ((col(spot) - (ncol(spot) + 1)/2))[L]
            sum(apply(IJ, 1, vecnorm))/n
        }
        cir <- function(spot, k) {
            L <- spot == k
            n <- sum(as.numeric(L))
            IJ <- matrix(NA, n, 2)
            IJ[, 1] <- row(spot)[L]
            IJ[, 2] <- col(spot)[L]
            IJ <- sweep(IJ, MARGIN = 2, FUN = "-", STATS = colMeans(IJ))
            d <- apply(IJ, 1, vecnorm)
            sum(as.numeric(d <= sqrt(n/pi)))/n
        }
        orderBYmean <- function(summary, q) {
            cl <- summary$classification
            u <- sort(unique(cl))
            G <- length(u)
            mu <- rep(NA, G)
            names(mu) <- as.character(u)
            for (k in 1:G) {
                n <- sum(as.numeric(cl == u[k]))
                mu[k] <- sum(spot[cl == u[k]]/n)
            }
            ord <- order(mu)
            cl <- -cl
            for (k in 1:G) cl[cl == -u[k]] <- ord[k]
            if (G == 2) {
                cl[cl == 2] <- 3
                nm <- names(mu)
                nm[nm == "2"] <- "3"
                names(mu) <- nm
            }
            mu <- mu[ord]
            list(cl = cl, mu = mu)
        }
        concomp <- function(binaryImage, nNeighbors = 4) {
            if (nNeighbors != 4 && nNeighbors != 8) 
                stop("4 or 8 neighbors")
            nrowImage <- nrow(binaryImage)
            ncolImage <- ncol(binaryImage)
            nPixels <- nrowImage * ncolImage
            M <- sum(as.numeric(binaryImage))
            K <- (1:nPixels)[as.logical(binaryImage)]
            binaryImage[K] <- 1:M
            Krev <- rev(K)
            change <- TRUE
            while (change) {
                change <- FALSE
                for (k in K) {
                  j <- ceiling(k/nrowImage)
                  i <- k - nrowImage * (j - 1)
                  l <- binaryImage[i, j]
                  if (nNeighbors == 4) {
                    iPairs <- c(i, i - 1, i, i + 1, i)
                    jPairs <- c(j - 1, j, j, j, j + 1)
                  }
                  else {
                    iPairs <- c(i - 1, i, i + 1, i - 1, i, i + 
                      1, i - 1, i, i + 1)
                    jPairs <- c(j - 1, j - 1, j - 1, j, j, j, 
                      j + 1, j + 1, j + 1)
                  }
                  I <- (iPairs > 0 & iPairs <= nrowImage)
                  I <- (jPairs > 0 & jPairs <= ncolImage) & I
                  N <- ((jPairs - 1) * nrowImage + iPairs)[I]
                  m <- min((binaryImage[N])[as.logical(binaryImage[N])])
                  binaryImage[i, j] <- m
                  change <- change || m != l
                }
                for (k in Krev) {
                  j <- ceiling(k/nrowImage)
                  i <- k - nrowImage * (j - 1)
                  l <- binaryImage[i, j]
                  if (nNeighbors == 4) {
                    iPairs <- c(i, i - 1, i, i + 1, i)
                    jPairs <- c(j - 1, j, j, j, j + 1)
                  }
                  else {
                    iPairs <- c(i - 1, i, i + 1, i - 1, i, i + 
                      1, i - 1, i, i + 1)
                    jPairs <- c(j - 1, j - 1, j - 1, j, j, j, 
                      j + 1, j + 1, j + 1)
                  }
                  I <- (iPairs > 0 & iPairs <= nrowImage)
                  I <- (jPairs > 0 & jPairs <= ncolImage) & I
                  N <- ((jPairs - 1) * nrowImage + iPairs)[I]
                  m <- min((binaryImage[N])[as.logical(binaryImage[N])])
                  binaryImage[i, j] <- m
                  change <- change || m != l
                }
            }
            tbl <- sort(unique(binaryImage[K]))
            I <- 2:length(tbl)
            for (i in I) binaryImage[binaryImage == tbl[i]] <- i
            binaryImage
        }
        plotSpotImage <- function(z, col = NULL, title = NULL, 
            one = FALSE) {
            if (is.null(col)) 
                col <- c("lightyellow", "gray", "black", "yellow", 
                  "blue", "red")
            if (!one) 
                par(mfrow = c(1, 1), pty = "m")
            nrowz <- nrow(z)
            ncolz <- ncol(z)
            if (!one) {
                PIN <- par()$pin
                if (ncolz < nrowz) 
                  par(pin = c(PIN[1] * (ncolz/nrowz), PIN[1]))
                else par(pin = c(PIN[1], PIN[1] * (nrowz/ncolz)))
            }
            z[z < 0] <- 3 + abs(z[z < 0])
            u <- sort(unique(as.vector(z)))
            if (length(col) < max(u)) 
                stop("not enough colors")
            L <- length(u)
            breaks <- rep(0, L + 1)
            breaks[1] <- u[1] - 0.5
            breaks[L + 1] <- u[L] + 0.5
            if (L > 1) {
                for (i in 2:L) breaks[i] <- (u[i - 1] + u[i])/2
            }
            col <- col[u]
            if (FALSE) {
                z[1, 1] <- min(breaks) - 0.5
                col <- c("blue", col)
                breaks <- c(min(breaks) - 1, breaks)
            }
            image(t(z[nrowz:1, ]), main = paste(title), xlab = "", 
                ylab = "", col = col, breaks = breaks, axes = FALSE)
            if (!one) 
                par(pin = PIN)
            invisible()
        }
        if (show) 
            plotBlockImage(sqrt(spot), title = "spot image", 
                one = TRUE)
        spot <- as.matrix(spot)
        if (hc) {
            BIC <- mclustBIC(as.vector(spot), G = 1:3, 
            initialization = list(hcPairs = hcE(data = as.vector(spot))))
        }
        else {
            BIC <- mclustBIC(as.vector(spot), G = 1:3)
        }
        spot[] <- orderBYmean(summary(BIC, 
                              data = as.vector(spot)), q)$cl
        G <- max(as.vector(spot))
        K <- unique(as.vector(spot))
        G <- length(K)
        if (G == 1) {
            spot[] <- 1
            if (show) 
                plotSpotImage(spot, title = "after clustering and final labeling", 
                  one = TRUE)
            if (show) 
                frame()
            if (show) 
                frame()
            if (show) 
                frame()
            return(spot)
        }
        if (show) 
            plotSpotImage(spot, title = "after clustering", one = TRUE)
        if (ccl) {
            for (k in K) {
                CC <- concomp(spot == k)
                tabl <- table(CC)[-1]
                namt <- names(tabl)
                M <- as.numeric(namt[tabl >= threshold])
                spot[spot == k & !as.logical(match(CC, M, nomatch = 0))] <- -k
            }
        }
        if (show) 
            plotSpotImage(spot, title = "after CC thresholding", 
                one = TRUE)
        if (any(spot > 0)) {
            tab <- table(spot[spot > 0])
            lab <- as.numeric(names(tab))
            if (length(tab) == 1) {
                if (lab == 2) 
                  spot[spot == 2] <- 1
            }
            else {
                fground <- lab[length(lab)]
                bground <- lab[1]
                if (fground <= bground) 
                  stop("thresholding bug")
                spot[spot > 0 & spot != fground & spot != bground] <- -9
                spot[spot == fground] <- 3
                spot[spot == bground] <- 1
            }
        }
        spot[spot < -1] <- 2
        spot[spot == -1] <- 1
        dis1 <- dis3 <- 0
        if (any(spot == 3)) {
            dis1 <- dis(spot, 1)
            dis3 <- dis(spot, 3)
            if (dis3 > dis1) 
                spot[spot == 3] <- 2
        }
        if (show) 
            plotSpotImage(spot, title = "final labeling", one = TRUE)
        spot
    }
    if (show) {
        par(mfrow = c(2, 2), pty = "m")
        if (length(R) == 1 || length(C) == 1) {
            par(ask = FALSE)
        }
        else {
            par(ask = TRUE)
        }
    }
    else {
        par(mfrow = c(1, 1), pty = "m")
    }
    m <- length(rowcut)
    n <- length(colcut)
    s <- chan1 + chan2
    L <- list(channel1 = list(foreground = list(mean = matrix(NA, 
        m, n), median = matrix(NA, m, n)), background = list(mean = matrix(NA, 
        m, n), median = matrix(NA, m, n))), channel2 = list(foreground = list(mean = matrix(NA, 
        m, n), median = matrix(NA, m, n)), background = list(mean = matrix(NA, 
        m, n), median = matrix(NA, m, n))))
    R <- if (is.null(R)) 
        2:m
    else R + 1
    C <- if (is.null(C)) 
        2:n
    else C + 1
    for (i in R) {
        I <- rowcut[i - 1]:(rowcut[i] - 1)
        for (j in C) {
            J <- colcut[j - 1]:(colcut[j] - 1)
            s[I, J] <- spotseg1(s[I, J], i = i - 1, j = j - 1, 
                threshold = threshold, hc = hc, show = show)
            if (show) {
                print(c(i, j) - 1)
                if (FALSE) {
                  pick <- menu("continue;0:exit", title = NULL)
                  if (!pick) 
                    return(NULL)
                }
            }
            fore <- s[I, J] == 3
            back <- s[I, J] == 1
            if (any(fore)) {
                L$channel1$foreground$mean[i, j] <- mean(chan1[I, 
                  J][fore], na.rm = TRUE)
                L$channel1$foreground$median[i, j] <- median(chan1[I, 
                  J][fore], na.rm = TRUE)
                L$channel2$foreground$mean[i, j] <- mean(chan2[I, 
                  J][fore], na.rm = TRUE)
                L$channel2$foreground$median[i, j] <- median(chan2[I, 
                  J][fore], na.rm = TRUE)
            }
            if (any(back)) {
                L$channel1$background$mean[i, j] <- mean(chan1[I, 
                  J][back], na.rm = TRUE)
                L$channel1$background$median[i, j] <- median(chan1[I, 
                  J][back], na.rm = TRUE)
                L$channel2$background$mean[i, j] <- mean(chan2[I, 
                  J][back], na.rm = TRUE)
                L$channel2$background$median[i, j] <- median(chan2[I, 
                  J][back], na.rm = TRUE)
            }
        }
    }
    L$channel1$foreground$mean <- L$channel1$foreground$mean[2:m, 
        2:n]
    L$channel1$foreground$median <- L$channel1$foreground$median[2:m, 
        2:n]
    L$channel1$background$mean <- L$channel1$background$mean[2:m, 
        2:n]
    L$channel1$background$median <- L$channel1$background$median[2:m, 
        2:n]
    L$channel2$foreground$mean <- L$channel2$foreground$mean[2:m, 
        2:n]
    L$channel2$foreground$median <- L$channel2$foreground$median[2:m, 
        2:n]
    L$channel2$background$mean <- L$channel2$background$mean[2:m, 
        2:n]
    L$channel2$background$median <- L$channel2$background$median[2:m, 
        2:n]
    s <- s[rowcut[1]:(rowcut[m] - 1), colcut[1]:(colcut[n] - 
        1)]
    structure(s, summaryStatistics = L, class = "spotseg")
}

