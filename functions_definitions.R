###############################################################################
# MCM estimators functions definitions when the cure status is partially known
###############################################################################

# The generalized product-limit estimator of S(t|x) 
survfitcurePK <- function(x, t, d, xinu, dataset = NULL, 
                          x0, h, local = TRUE) {
  dfr <-
    if (missing(dataset)){
      na.omit(data.frame(x, t, d, xinu))
    }else{
      na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)), 
                          deparse(substitute(d)), deparse(substitute(xinu)))])}
  names(dfr) <- c("x", "t", "d", "xinu")
  n_row <- dim(dfr)[1]
  ord.dfr <- dfr[order(dfr$t, - dfr$d, dfr$xinu),]
  dfr$x <- as.numeric(ord.dfr$x)
  dfr$t <- as.numeric(ord.dfr$t)
  dfr$d <- as.integer(ord.dfr$d)
  dfr$xinu <- as.integer(ord.dfr$xinu)
  
  ordx0 <- order(x0)
  x0 <- as.numeric(x0[ordx0])
  lx0 <- length(x0)
  
  if (!local && missing(h)) warning("Option 'local = FALSE' overridden: with missing 'h'")
  
  if (missing(h)) {
    h <- survfitcurePKboot(dfr = dfr, x0 = x0)
  } else{
    if (local) {
      stopifnot("When 'local = TRUE', 'x0' and 'h' must have the same length" = lx0 == length(h))
      h <- as.numeric(h[ordx0])
    }
    else {
      warning("Option 'local = FALSE' global bandwidth(s) is used. Authors recommend 'local = TRUE'.")
      h <- as.numeric(h)
    }
  }
  lh <- length(h)
  if(local) {
    surv <- matrix(0, n_row, lx0)
    for (ix0 in 1:lx0) {
      res <- 1
      xh <- (x0[ix0] - dfr$x)/h[ix0]
      K <- (abs(xh) < 1) * 0.75 * (1 - xh^2)
      BCx <- sum(K * (dfr$xinu == 1))
      sumBx <- rev(cumsum(rev(K * (dfr$xinu == 0))))
      for (i in 1:n_row) {
        denom <- sumBx[i] + BCx
        if (!is.na(denom) & denom > 0 & (dfr$d[i] == 1))
          res <- res * (1 - K[i]/denom)
        surv[i, ix0] <- res
      }
    }
  } else{
    if(lh == 1L){
      surv <- matrix(0, n_row, lx0)
      for (ix0 in 1:lx0) {
        res <- 1
        xh <- (x0[ix0] - dfr$x)/h
        K <- (abs(xh) < 1) * 0.75 * (1 - xh^2)
        BCx <- sum(K * (dfr$xinu == 1))
        sumBx <- rev(cumsum(rev(K * (dfr$xinu == 0))))
        for (i in 1:n_row) {
          denom <- sumBx[i] + BCx
          if (!is.na(denom) & denom > 0 & (dfr$d[i] == 1))
            res <- res * (1 - K[i]/denom)
          surv[i, ix0] <- res
        }
      }
    }
    
    if(lh > 1L){
      surv <- array(0, dim = c(n_row, lh, lx0))
      res <- matrix(1, lh, lx0)
      for (ix0 in 1:lx0) {
        for (ih in 1:lh) {
          xh <- (x0[ix0] - dfr$x)/h[ih]
          K <- (abs(xh) < 1) * 0.75 * (1 - xh^2)
          BCx <- sum(K * (dfr$xinu == 1))
          sumBx <- rev(cumsum(rev(K * (dfr$xinu == 0))))
          for (i in 1:n_row) {
            denom <- sumBx[i] + BCx
            if (!is.na(denom) & denom > 0 & (dfr$d[i] == 1))
              res[ih, ix0] <- res[ih, ix0] * (1 - K[i]/denom)
            surv[i, ih, ix0] <- res[ih, ix0]
          }
        }
      }
    }
    
  }
  
  return(list(S = surv, p = 1-res))
  
}


# Uncondional estimator for the latency function

survfitcurePK_un <- function (dataset = dataset) {
  dfr <- dataset
  names(dfr) <- c("t", "d", "xinu")
  n_row <- dim(dfr)[1]
  ord.dfr <- dfr[order(dfr$t, - dfr$d, dfr$xinu),]
  t <- as.numeric(ord.dfr$t)
  d <- as.integer(ord.dfr$d)
  xinu <- as.integer(ord.dfr$xinu)
  
  n <- nrow(dfr)

  cum.nu <- cumsum(xinu)  # Number of known cures up to time Ti
  S <- rep(1, n)
  for (i in 2:n) {
    
    if(d[i]==0) {S[i] <- S[i-1]}
    
    if(d[i]==1) {S[i] <- S[i-1] * (1 - 1/(n - i + 1 + cum.nu[i-1]))}
  }
  
  p <- 1 - min(S)
  return(list(S, p, t))
  
}
###############################################################################
# Bootstrap bandwidth calculations
###############################################################################
# Bootstrap parameters
B <-  50; hbound <-  c(0.1, 3); hl <-  10L; nnfrac = 0.25

## The functions to compute Nadaraya-Watson weights with Epanechnikov kernel.
kernel <- function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
weights.function <- function(x, x0, K, h) {K((x - x0)/h)}

survfitcurePKinhbootmise <- function(dfr, x0, g){
  ord.dfr <- dfr[order(dfr$t, - dfr$d, dfr$xinu),]
  n_row <- dim(dfr)[1]
  n_col <- dim(dfr)[2]
  hbound <- IQR(dfr$x) / 1.349 * hbound
  lhgrid <- hl
  steph <- (hbound[2] / hbound[1])^(1 / lhgrid)
  hgrid <- as.numeric(hbound[1] * steph^seq(0, lhgrid, length.out = lhgrid))
  out <- survfitcurePK(x, t, d, xinu, ord.dfr, x0, g)
  Sg <- out$S
  pg <-  out$p
  dfrboot <-  data.frame(matrix(0, nrow = n_row, ncol = n_col,
                                dimnames = list(NULL, c("x", "t", "d", "xinu"))))
  if(pg == 0) {
    S0g <- pg*Sg
    tmax <- quantile(ord.dfr[, 2], probs = 1 - 0.9)[[1]]
  } else {
    S0g <- ifelse(((Sg - (1 - pg))/pg) > 0 & ord.dfr$t > tail(ord.dfr$t[ord.dfr$d == 1], 1), 0,
                   ifelse(((Sg - (1 - pg)) / pg) < 0, 0, ((Sg - (1 - pg)) / pg)))
    tmax <- data.table::first(ord.dfr$t[S0g == data.table::first(DescTools::Closest(S0g, 1 - 0.9, na.rm = T))])
    
  }
  
  tgrid <- seq(0, tmax, length.out = 100L)
  weights.matrix <- outer(dfr$x, dfr$x, weights.function, K = kernel, h = g)
  
  dfrboot[ , 1] <-  dfr$x
  for (i in 1:n_row) {
    dfrboot[i, 2:4] <- dfr[sample(nrow(dfr), size = 1, prob = weights.matrix[i, ]), 2:4]
  }
  ord.dfrboot <- dfrboot[order(dfrboot$t, - dfrboot$d, dfrboot$xinu), ]
  bootmse <- bootmise <- numeric(lhgrid)
  
  for (ih in 1:lhgrid) {
    out1 <- survfitcurePK(x, t, d, xinu, dfrboot, x0, hgrid[ih])
    
    # Bootstrap MISE
    Sh_boot <- out1$S
    Sh.boot.grid <- approxfun(x = c(0, ord.dfrboot$t), y = c(1, Sh_boot), method = "constant", f = 0, rule = 2, ties = "ordered")
    Sg.grid <- approxfun(x = c(0, ord.dfr$t), y = c(1, Sg), method = "constant", f = 0, rule = 2, ties = "ordered")
    diff <- (Sh.boot.grid(tgrid) - Sg.grid(tgrid))^2
    diff.grid <- approxfun(x = tgrid, y = diff, method = "constant", f = 0, rule = 2, ties = "ordered")
    se_boot <- function(t) diff.grid(t)
    
    bootmise[ih] <- bootmise[ih] + as.numeric(integrate(se_boot, 0 + .Machine$double.eps^0.3, tmax, rel.tol = 0.1, subdivisions = 2000)[1])
    
    # Bootstrap MSE
    ph <- out1$p
    bootmse[ih] <- bootmse[ih]+(ph - pg)^2
  }
  return(list(bootmise = bootmise, bootmse = bootmse))
}


survfitcurePKboot <- function(dfr, x0, ncpus = 2L) {
  names(dfr) <- c("x", "t", "d", "xinu")
  dfr$x <- as.numeric(dfr$x)
  dfr$t <- as.numeric(dfr$t)
  dfr$d <- as.integer(dfr$d)
  dfr$xinu <- as.integer(dfr$xinu)
  x0 <- as.numeric(x0)
  lx0 <- length(x0)
  n_row <- dim(dfr)[1]
  n_col <- dim(dfr)[2]
  hbound <- IQR(dfr$x) / 1.349 * hbound
  lhgrid <- hl
  steph <- (hbound[2] / hbound[1])^(1 / lhgrid)
  hgrid <- as.numeric(hbound[1] * steph^seq(0, lhgrid, length.out = lhgrid))
  k <- floor(n_row*nnfrac)
  
  
  #mise <- matrix(0, lhgrid, lx0)
  min_mise <- min_mse <- hboot_S <- hboot_p <- numeric(lx0)
  
  pilot <- npcure::hpilot(dfr$x, x0, nnfrac)
  
  ## Compute  MISE with parallelization --> %dopar%
  
  cl <- parallel::makeCluster(ncpus)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(b) setTxtProgressBar(pb, b)
  opts <- list(progress = progress)
  out <- foreach::foreach(b = 1 : B,
                          .packages = c("npcure", "DescTools","data.table", "doParallel"),
                          .export = ls(globalenv()),
                          .options.snow = opts,
                          .verbose = TRUE)%dopar%{
                            survfitcurePKinhbootmise(dfr, x0, g = pilot)
                          }
  miseres <- Reduce("+", lapply(out, "[[", 1))
  mise <- miseres / B
  min_mise <- min(mise)
  hboot_S <- hgrid[mean(which(min_mise == mise))]
  
  mseres <- Reduce("+", lapply(out, "[[", 2))
  mse <- mseres / B
  min_mse <- min(mse)
  hboot_p <- hgrid[mean(which(min_mse == mse))]
  
  close(pb)
  parallel::stopCluster(cl)
  
  
  return(list(hboot_S = hboot_S, hboot_p = hboot_p))
}

# Latency bootstrap computation
latencyinbootfun <- function(dfr = dfr,   x0 = x0, hS = hS, hp = hp){
  
  ord.dfr <- dfr[order(dfr$t, - dfr$d, dfr$xinu),]
  n_row <- dim(dfr)[1]
  n_col <- dim(dfr)[2]
  
  # Bandwidth interval
  hbound <- IQR(dfr$x) / 1.349 * hbound
  lhgrid <- hl
  steph <- (hbound[2] / hbound[1])^(1 / lhgrid)
  hgrid <- as.numeric(hbound[1] * steph^seq(0, lhgrid, length.out = lhgrid))
  
  dfrboot <-  data.frame(matrix(0, nrow = n_row, ncol = n_col,
                                dimnames = list(NULL, c("x", "t", "d", "xinu"))))
  
  ## latency estimator with pilot bandwidth computed from Shct(t|x) and phc(x)
  Sg <- survfitcurePK(x, t, d, xinu, dfr, x0, hS)$S
  pg <-  survfitcurePK(x, t, d, xinu, dfr, x0, hp)$p
  
  if(pg == 0) {
    S0g <- pg*Sg
    tmax <- quantile(ord.dfr[, 2], probs = 1 - 0.9)[[1]]
  } else {
    S0g <- ifelse(((Sg - (1 - pg))/pg) > 0 & ord.dfr$t > tail(ord.dfr$t[ord.dfr$d == 1], 1), 0,
                  ifelse(((Sg - (1 - pg)) / pg) < 0, 0, ((Sg - (1 - pg)) / pg)))
    tmax <- data.table::first(ord.dfr$t[S0g == data.table::first(DescTools::Closest(S0g, 1 - 0.9, na.rm = T))])
    
  }  
  
  tgrid <- seq(0, tmax, length.out = 100L)
  weights.matrix <- outer(dfr$x, dfr$x, 
                          weights.function, K = kernel, h = hS)
  
  dfrboot[ , 1] <-  dfr$x
  for (i in 1:n_row) {
    dfrboot[i, 2:4] <- dfr[sample(nrow(dfr), size = 1, prob = weights.matrix[i, ]), 2:4]
  }
  ord.dfrboot <- dfrboot[order(dfrboot$t, - dfrboot$d, dfrboot$xinu), ]
  bootmise <- matrix(0, lhgrid, lhgrid) 
  
  
  for (ih in 1:lhgrid) {
    for (ij in 1:lhgrid) {
      
      Sh <- survfitcurePK(x, t, d, xinu, dfr, x0, hgrid[ij])$S
      ph <- survfitcurePK(x, t, d, xinu, dfr, x0, hgrid[ih])$p
      
      if(ph == 0) {
        S0h <- ph*Sh
      } else {
        S0h <- ifelse(((Sh-(1-ph))/ph) > 0 & ord.dfrboot[,2] > tail(ord.dfrboot[, 2][ord.dfrboot[, 3] == 1], 1), 0,
                      ifelse(((Sh - (1 - ph))/ph) < 0,0,((Sh - (1 - ph)) / ph)))
      }
      
      
      S0h.grid <- approxfun(x = c(0, ord.dfrboot[, 2]), y = c(1, S0h),
                            method = "constant", f = 0, rule = 2, ties = "ordered")
      S0g.grid <- approxfun(x = c(0, ord.dfr$t), y = c(1, S0g), 
                            method = "constant", f = 0, rule = 2, ties = "ordered")
      
      diff <- (S0h.grid(tgrid) - S0g.grid(tgrid))^2
      diff.grid <- approxfun(x = tgrid, y = diff, method = "constant", f = 0, rule = 2, ties = "ordered")
      
      se <- function(t) {diff.grid(t)}
      bootmise[ij, ih] <-  as.numeric(integrate(se, 0 + .Machine$double.eps^0.3, 
                                                tmax, rel.tol = 0.1, subdivisions = 2000)[1])
      
    } # for each ij
  } # for each ih
  return(bootmise = bootmise)
}

latencybootfun <- function(dfr, x0, ncpus = 2L){
  names(dfr) <- c("x", "t", "d", "xinu")
  dfr$x <- as.numeric(dfr$x)
  dfr$t <- as.numeric(dfr$t)
  dfr$d <- as.integer(dfr$d)
  dfr$xinu <- as.integer(dfr$xinu)
  x0 <- as.numeric(x0)
  lx0 <- length(x0)
  n_row <- dim(dfr)[1]
  n_col <- dim(dfr)[2]
  
  # Bandwidth interval
  hbound <- IQR(dfr$x) / 1.349 * hbound
  lhgrid <- hl
  steph <- (hbound[2] / hbound[1])^(1 / lhgrid)
  hgrid <- as.numeric(hbound[1] * steph^seq(0, lhgrid, length.out = lhgrid))
  
  ## Compute pilot bandwidths to be used below
  hSp <- survfitcurePKboot(dfr, x0)
  
  ## Compute  MISE with parallelization --> %dopar%
  cl <- parallel::makeCluster(ncpus)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(b) setTxtProgressBar(pb, b)
  opts <- list(progress = progress)
  out <- foreach::foreach(b = 1 : B,
                          .packages = c("npcure", "DescTools","data.table", "doParallel"),
                          .export = ls(globalenv()),
                          .options.snow = opts,
                          .verbose = TRUE)%dopar%{
                            latencyinbootfun(dfr, x0, hS = hSp$hboot_S, hp = hSp$hboot_p)
                          }
  miseres <- Reduce("+", out)
  mise <- miseres / B
  min_mise <- min(mise)
  ind <- colMeans(which(min_mise == mise, arr.ind = T))
  hboot <- hgrid[ind]
  
  close(pb)
  parallel::stopCluster(cl)
  
  return(hboot = hboot)
}









