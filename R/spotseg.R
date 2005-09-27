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
        EMclust1 <- function(data, G, eps, tol, itmax, warnSingular = FALSE, 
            ...) {
            equalPro <- FALSE
            dimData <- dim(data)
            oneD <- is.null(dimData) || length(dimData[dimData > 
                1]) == 1
            if (!oneD) 
                stop("data must be a vector")
            if (missing(eps)) 
                eps <- .Mclust$eps
            if (missing(tol)) 
                tol <- .Mclust$tol
            if (missing(itmax)) 
                itmax <- .Mclust$itmax
            itmax[is.infinite(itmax)] <- .Machine$integer.max
            data <- as.vector(data)
            n <- length(data)
            p <- 1
            emModelNames <- c("E", "V")
            m <- length(emModelNames)
            if (missing(G)) {
                G <- 1:9
            }
            else {
                G <- sort(G)
            }
            if (any(G) <= 0) 
                stop("G must be positive")
            l <- length(G)
            Glabels <- as.character(G)
            BIC <- matrix(0, nrow = l, ncol = m, dimnames = list(Glabels, 
                emModelNames))
            maxG <- max(G)
            paramList <- rep(list(list(E = list(mu = rep(0, maxG), 
                sigmasq = 0, pro = rep(0, maxG)), V = list(mu = rep(0, 
                maxG), sigmasq = rep(0, maxG), pro = rep(0, maxG)))), 
                nrow(BIC))
            names(paramList) <- Glabels
            if (G[1] == 1) {
                for (modelName in emModelNames) {
                  mle <- mvn(modelName = modelName, data = data)
                  BIC[1, modelName] <- bic(modelName = modelName, 
                    loglik = mle$loglik, n = n, d = p, G = 1, 
                    equalPro = equalPro)
                  paramList[[1]][[modelName]] <- c(mle[c("mu", 
                    "sigmasq")], list(pro = 1))
                }
                G <- G[-1]
                Glabels <- Glabels[-1]
            }
            if (l > 1) {
                for (i in seq(along = G)) {
                  k <- G[i]
                  cl <- rep(0, length(data))
                  Q <- quantile(data, probs = (0:k)/k)
                  Q[1] <- Q[1] - 1
                  for (j in 1:k) cl[data > Q[j] & data <= Q[j + 
                    1]] <- j
                  z <- unmap(cl)
                  for (modelName in emModelNames) {
                    mle <- me(modelName = modelName, data = data, 
                      z = z, eps = eps, tol = tol, itmax = itmax, 
                      equalPro = equalPro, warnSingular = warnSingular)
                    BIC[Glabels[i], modelName] <- bic(modelName = modelName, 
                      loglik = mle$loglik, n = n, d = p, G = G[i], 
                      equalPro = equalPro)
                    paramList[[Glabels[i]]][[modelName]] <- mle[c("mu", 
                      "sigmasq", "pro")]
                  }
                }
            }
            structure(BIC, eps = eps, tol = tol, itmax = itmax, 
                warnSingular = warnSingular, paramList = paramList, 
                args = as.list(match.call())[-1], class = "EMclust1")
        }
        summary.EMclust1 <- function(object, data, G, modelNames, 
            ...) {
            x <- object
            paramList <- attr(object, "paramList")
            n <- if (is.null(dimData <- dim(data))) 
                length(data)
            else dimData[1]
            eps <- attr(x, "eps")
            tol <- attr(x, "tol")
            itmax <- attr(x, "itmax")
            equalPro <- attr(x, "equalPro")
            warnSingular <- attr(x, "warnSingular")
            class(x) <- attr(x, "args") <- attr(x, "warnSingular") <- NULL
            if (missing(G)) {
                G <- as.numeric(dimnames(x)[[1]])
            }
            else {
                G <- sort(G)
            }
            Glabels <- as.character(G)
            if (missing(modelNames)) 
                modelNames <- dimnames(x)[[2]]
            x <- x[Glabels, modelNames, drop = FALSE]
            X <- is.na(x)
            if (all(X)) 
                stop("none of the selected models could be fitted")
            x[X] <- -.Machine$double.xmax
            l <- nrow(x)
            m <- ncol(x)
            best <- max(x)
            rowsBest <- (matrix(rep(1:l, m), l, m)[x == best])[1]
            colsBest <- (matrix(rep(1:m, rep(l, m)), l, m)[x == 
                best])[1]
            namesBest <- dimnames(x[rowsBest, colsBest, drop = FALSE])
            bestG <- namesBest[[1]]
            maxG <- max(G)
            minG <- min(G)
            if (minG != maxG) {
                if (bestG == maxG) {
                }
                else if (minG != 1 && bestG == min(G)) {
                }
            }
            bestModel <- namesBest[[2]]
            if (min(l, m) > 1) {
                M <- modelNames[modelNames != bestModel]
                y <- x[, M]
                other <- max(y)
                otherG <- (matrix(rep(Glabels, m - 1), l, m - 
                  1)[y == other])[1]
                otherModel <- (matrix(rep(M, rep(l, m - 1)), 
                  l, m - 1)[y == other])[1]
                y <- x[, bestModel]
                w <- y[y != best]
                if (length(w) == l - 1) {
                  same <- max(w)
                  sameG <- (Glabels[y == same])[1]
                }
                else {
                  same <- best
                  sameG <- (Glabels[y == same])[2]
                }
                nam1 <- paste(bestModel, bestG, sep = ",")
                nam2 <- paste(bestModel, sameG, sep = ",")
                nam3 <- paste(otherModel, otherG, sep = ",")
                bestBICs <- c(nam1 = best, nam2 = same, nam3 = other)
                names(bestBICs) <- c(nam1, nam2, nam3)
            }
            else if (l != 1) {
                w <- x[x != best]
                if (length(w) == l - 1) {
                  same <- max(w)
                  sameG <- (Glabels[x == same])[1]
                }
                else {
                  same <- best
                  sameG <- (Glabels[x == same])[2]
                }
                nam1 <- paste(bestModel, bestG, sep = ",")
                nam2 <- paste(bestModel, sameG, sep = ",")
                bestBICs <- c(nam1 = best, nam2 = same)
                names(bestBICs) <- c(nam1, nam2)
            }
            else if (m != 1) {
                M <- (1:m)[modelNames == bestModel]
                y <- x[, -M]
                other <- max(y)
                otherG <- (matrix(rep(Glabels, m - 1), l, m - 
                  1)[y == other])[1]
                otherModel <- (matrix(rep(modelNames[-M], rep(l, 
                  m - 1)), l, m - 1)[y == other])[1]
                nam1 <- paste(bestModel, bestG, sep = ",")
                nam3 <- paste(otherModel, otherG, sep = ",")
                bestBICs <- c(best, other)
                names(bestBICs) <- c(nam1, nam3)
            }
            else {
                nam1 <- paste(bestModel, bestG, sep = ",")
                bestBICs <- best
                names(bestBICs) <- nam1
            }
            bestBICs[bestBICs == -.Machine$double.xmax] <- NA
            out <- paramList[[bestG]][[bestModel]]
            if (as.numeric(bestG) == 1) {
                structure(c(list(bic = bestBICs, G = bestG, classification = rep(1, 
                  n), uncertainty = rep(0, n)), out), class = "summary.EMclust1")
            }
            else {
                z <- do.call("estep", c(list(data = data, modelName = bestModel), 
                  out))$z
                structure(c(list(bic = bestBICs, G = bestG, classification = map(z), 
                  uncertainty = 1 - apply(z, 1, max), modelName = bestModel), 
                  out), class = "summary.EMclust1")
            }
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
            BIC <- EMclust(as.vector(spot), G = 1:3, hcPairs = hc(modelName = "E", 
                data = as.vector(spot)))
        }
        else {
            BIC <- EMclust1(as.vector(spot), G = 1:3)
        }
        spot[] <- orderBYmean(summary(BIC, data = as.vector(spot)), 
            q)$cl
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

