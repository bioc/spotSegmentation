"spotgrid" <-
function (chan1, chan2, rows = NULL, cols = NULL, span = NULL, 
    show = FALSE) 
{
    signal <- chan1+chan2
    s <- min(signal)
    r <- min(signal[signal > 0])
    if (s <= 0) {
       signal <- signal - s + 1
    }
    else if (r < 1) {
        signal <- signal + 1
    }
    signal <- log(signal)
    spotgridPeaks <- function(series, span = 3, seed = 0) {
        if (!(span%%2)) 
            span <- span + 1
        d <- sort(diff(sort(series)))
        if (!d[1]) {
            d <- d[as.logical(d)][1]
            set.seed(seed)
            noise <- runif(length(series), min = -d/2, max = d/2)
            series <- series + noise
        }
        zmaxcol <- apply(embed(series, span)[, span:1], 1, function(x) {
            if (length(y <- seq(along = x)[x == max(x)]) == 1) 
                y
            else 0
        })
        halfspan <- (span - 1)/2
        c(rep(FALSE, halfspan), zmaxcol == 1 + halfspan)
    }
    isum <- colSums(signal)
    gridcomp <- function(isum, nspots, span) {
        N <- length(isum)
        if (!(span%%2)) 
            span <- span + 1
        Peaks <- spotgridPeaks(isum, span)
        Vales <- spotgridPeaks(-isum, span)
        iPeaks <- (1:N)[Peaks]
        iVales <- (1:N)[Vales]
        nPeaks <- length(iPeaks)
        nVales <- length(iVales)
        if (FALSE) {
            plot(1:N, isum, type = "l", xlab = "", ylab = "")
            points(iPeaks, isum[iPeaks], pch = "P")
            points(iVales, isum[iVales], pch = "V")
        }
        V <- rep(0, nPeaks - 1)
        for (i in 2:nPeaks) {
            K <- iPeaks[i - 1]:iPeaks[i]
            I <- isum[K]
            J <- K[I == min(I)]
            if (length(J) > 1) {
                M <- match(J, iVales, nomatch = 0)
                J <- if (any(M)) 
                  (iVales[M])[1]
                else J[1]
            }
            V[i - 1] <- J
        }
        P <- rep(0, nVales - 1)
        for (i in 2:nVales) {
            K <- iVales[i - 1]:iVales[i]
            I <- isum[K]
            J <- K[I == max(I)]
            if (length(J) > 1) {
                M <- match(J, iPeaks, nomatch = 0)
                J <- if (any(M)) 
                  (iPeaks[M])[1]
                else J[1]
            }
            P[i - 1] <- J
        }
        iPeaks <- unique(sort(c(iPeaks, P)))
        iVales <- unique(sort(c(iVales, V)))
        nPeaks <- length(iPeaks)
        nVales <- length(iVales)
        if (nPeaks < nspots) {
##            warning("fewer peaks than spots")
            nspots <- nPeaks
        }
        peakvals <- rep(0, nPeaks)
        for (i in 1:nPeaks) {
            m <- iPeaks[i]
            J <- iVales[iVales < m]
            K <- iVales[iVales > m]
            lSum <- rSum <- 0
            d <- 0
            if (length(J)) {
                lSum <- isum[J[length(J)]]
                d <- d + 1
            }
            if (length(K)) {
                rSum <- isum[K[1]]
                d <- d + 1
            }
            peakvals[i] <- isum[m] - (lSum + rSum)/d
        }
        i <- 0
        smax <- sum(sort(peakvals))
        k <- nPeaks - nspots
        if (k) {
            smax <- s <- smax - sum(peakvals[-(1:nspots)])
            for (j in 1:k) {
                d <- (peakvals[j + nspots] - peakvals[j])
                s <- s + d
                if (s > smax) {
                  smax <- s
                  i <- j
                }
            }
        }
        span <- max(span, max(diff(iPeaks[i + 1:nspots])))
        j <- 1
        while (TRUE) {
            if (iVales[j] > iPeaks[i + 1] && iVales[j] < iPeaks[i + 
                2]) 
                break
            j <- j + 1
        }
        index <- rep(0, nspots + 1)
        index[2:nspots] <- iVales[j:(j + nspots - 2)]
        halfspan <- floor(span/2)
        index[1] <- iPeaks[i + 1] - halfspan
        index[nspots + 1] <- iPeaks[nspots + i] + halfspan
        if (FALSE) {
            par(ask = T)
            plot(1:N, isum, type = "l", xlab = "", ylab = "")
            points(iPeaks, isum[iPeaks], pch = "P")
            points(iVales, isum[iVales], pch = "V")
            abline(v = index[1], col = "red")
            abline(v = index[length(index)], col = "red")
            for (i in 1:length(index)) {
                abline(v = index[i], col = "red")
            }
        }
        index
    }
    if (show) {
        plotBlockImage(signal)
    }
    rowcut <- colcut <- NA
    if (!is.null(rows) && rows > 0) {
        if (is.null(span)) {
            span <- floor(nrow(signal)/rows)
            if (!(span%%2)) 
                span <- span - 1
        }
        rowcut <- gridcomp(rowSums(signal), rows, if (length(span) == 
            2) 
            span[1]
        else span)
        if (length(rowcut) < rows+1) {
          span <- min(diff(rowcut))
          rowcut <- gridcomp( rowSums(signal), rows, span) 
        }
        if (length(rowcut) < rows+1) warning("fewer peaks than spots")
        if (show && (is.null(cols) || cols <= 0)) {
            chan1[] <- 0
            chan1[rowcut, ] <- 1
            contour(z = t(chan1[nrow(signal):1, ]), nlevels = 1, 
                levels = 1, drawlabels = FALSE, col = "red", 
                add = TRUE)
        }
    }
    if (!is.null(cols) && cols > 0) {
        if (is.null(span)) {
            span <- floor(ncol(signal)/cols)
            if (!(span%%2)) 
                span <- span - 1
        }
        colcut <- gridcomp(colSums(signal), cols, if (length(span) == 
            2) 
            span[2]
        else span)
        if (length(colcut) < cols+1) {
          span <- min(diff(colcut))
          colcut <- gridcomp( colSums(signal), cols, span) 
        }
        if (length(colcut) < cols+1) warning("fewer peaks than spots")
        if (show && (is.null(rows) || rows <= 0)) {
            chan1[] <- 0
            chan1[, colcut] <- 1
            contour(z = t(chan1[nrow(signal):1, ]), nlevels = 1, 
                levels = 1, drawlabels = FALSE, col = "red", 
                add = TRUE)
        }
    }
    if (show && !is.na(rowcut) && !is.na(colcut)) {
        chan1[] <- 0
        chan1[rowcut, colcut[1]:colcut[length(colcut)]] <- 1
        chan1[rowcut[1]:rowcut[length(rowcut)], colcut] <- 1
        contour(z = t(chan1[nrow(chan1):1, ]), nlevels = 1, levels = 1, 
            drawlabels = FALSE, col = "red", add = TRUE)
    }
    list(rowcut = rowcut, colcut = colcut)
}

