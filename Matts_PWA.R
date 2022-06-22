pwa <- function(pw, filt = FALSE, plot = FALSE) {
  
  # Load functions
  fsg721 <- function(x) {
    # 1st derivative with SG filter
    #2nd order polynomial
    C = c(0.107143, 0.071429, 0.035714)
    B = integer(7)
    for (i in 1:3) {
      B[i] = C[i]
    }
    B[4] = 0.0
    for (i in 5:7) {
      B[i] = -C[8 - i]
    }
    A = c(1, 0)
    s = length(x)
    dx = signal::filter(B, A, x)
    dx = c(dx[7], dx[7], dx[7], dx[7:s], dx[s], dx[s], dx[s])
  }
  
  Tintersect <- function(wf, plot = FALSE) {
    
    k1 <- wf - min(wf[1:which.max(wf)])
    xvar <- (1:length(k1) - 1)
    spl <- smooth.spline(k1 ~ xvar)
    newx <- which.max(diff(k1))
    pred0 <- predict(spl, x = newx, deriv = 0)
    pred1 <- predict(spl, x = newx, deriv = 1)
    yint <- pred0$y - (pred1$y * newx)
    xint <- (-yint / pred1$y)
    
    if(isTRUE(plot)) {
      plot(xvar, k1, ylim=c(min(k1)-20, max(k1)))
      abline(h=min(k1), col="red", lty=3)
      lines(spl, col="red") 
      lines(xvar, yint + pred1$y*xvar, col="green", lwd=2)
      points(pred0,col="red", pch=8, lwd=2) 
      points(xint, 0, col="red", pch=8, lwd=2) 
      abline(v=xint, lty=3)
    }
    
    return(xint)
    
  }
  
  dicrotic <- function(pw, plot = FALSE) {
    
    # Get derivatives
    dp1 <- fsg721(pw)
    dp2 <- fsg721(fsg721(pw))
    dp3 <- fsg721(fsg721(fsg721(pw)))
    
    
    # FIND DICROTIC DEPRESSION ------------------------------------------------
    
    # Isolate notch area with 2nd and 3rd derivatives
    nni <- which.min(dp1)
    end <- length(pw)
    
    # End index without potential perturbation at end diastole  
    end2 <- end * .9
    
    # Dicrotic notch from local dp2 max
    dic <- which.max(dp2[nni:end2]) + nni - 1

    # plot(pw, type="l", lwd=2)
    # par(new=T)
    # plot(dp2, type='o',col="grey")
    # abline(v = dic, h = 0)
    
    
    # FIND DICROTIC PEAK ------------------------------------------------------
    
    end3 <- ((end - dic) * .6) + dic # 60% of diastolic duration
    #abline(v=end3, lty=2, col=2)
    
    if(sum(dp2[dic:end3] < 0) < 1) {
      dia <- 9999
    } else {
      dia <- which.min(dp2[dic:end3]) + dic - 1
    }
    
    # plot(pw, type="l", lwd=2)
    # par(new=T)
    # plot(dp2, type='o',col="grey")
    # abline(v = c(dic, dia), h = 0)
    
    
    # PLOTS -------------------------------------------------------------------
    
    if(isTRUE(plot)) {
      plot(pw, type = "l", lwd=2, ylab="BP (mmHg)")
      abline(v=c(dic, dia), col="grey", lty=3, lwd=2)
      mtext(c("Ed", "P3"), side = 3, at = c(dic,dia))
    }
    
    return(data.frame(dicrotic_notch = dic, 
                      dicrotic_peak = dia))
    
  }
  
  RootSpline1 <- function (x, y, y0 = 0, verbose = TRUE) {
    if (is.unsorted(x)) {
      ind <- order(x)
      x <- x[ind]; y <- y[ind]
    }
    z <- y - y0
    ## which piecewise linear segment crosses zero?
    k <- which(z[-1] * z[-length(z)] <= 0)
    ## analytical root finding
    xr <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
    ## make a plot?
    if (verbose) {
      plot(x, y, "l"); abline(h = y0, lty = 2)
      points(xr, rep.int(y0, length(xr)))
    }
    ## return roots
    xr
  }
  
  low.pass <- function(y, fq, do.plot = FALSE) {
    
    # Second order low pass filter
    # Removes high frequency components below fq
    # y = a numeric vector, typically a tree-ring series.
    # fq = a numeric vector giving frequency or period of the filter.
    # Rp = a numeric value giving the dB for the passband ripple.
    
    if (any(is.na(y))) stop("y contains NA")
    
    ## n = a numeric value giving the order of the filter. 
    ## Larger numbers create steeper fall off.
    n = 4
    
    if (any(fq>1)) {
      f <- 1/fq
      p <- fq
    } else {
      p <- 1/fq
      f <- fq
    }
    
    # sort f in case it's passed in backwards
    f <- sort(f)
    
    filt <- signal::butter(
      n = n,
      W = f * 2,
      type = "low",
      plane = "z"
    )
    
    # remove mean
    yAvg <- mean(y)
    y <- y - yAvg
    
    # pad the data to twice the max period
    pad <- max(p) * 2
    ny <- length(y)
    
    # pad the data
    yPad <- c(y[pad:1], y, y[ny:(ny - pad)])
    
    # run the filter
    yFilt <- signal::filtfilt(filt, yPad)
    
    # unpad the filtered data
    yFilt <- yFilt[(pad + 1):(ny + pad)]
    
    # return with mean added back in
    filt.sig <- yFilt + yAvg
    
    if(isTRUE(do.plot)){
      
      # plot results
      plot(filt.sig,
           type = "l",
           lwd = 2)
    }
    
    # return filtered signal
    return(filt.sig)
    
  }
  
  # Low pass waveform
  if (isTRUE(filt)) {
    pw <- low.pass(pw, 10, do.plot = F)
  }
  
  # Create derivatives
  d1 <- fsg721(pw)
  d2 <- fsg721(fsg721(pw))
  d3 <- fsg721(fsg721(fsg721(pw)))
  d4 <- fsg721(fsg721(fsg721(fsg721(pw))))
  
  # Some additional calcs
  end <- length(pw)
  foot <- Tintersect(pw) - 1
  if(foot < 1) {foot <- 1}
  maxpi <- which.max(pw)
  notchdat <- dicrotic(pw)
  notch <- notchdat$dicrotic_notch
  notchpeak <- notchdat$dicrotic_peak
  
  # plot(pw,type="o")
  # abline(v=c(foot, notch, notchpeak))
  
  # Create time
  time <- (0:(length(pw)-1)) / 200
  X <- 1:length(pw)
  
  # get zero crossing of 4th derivative
  zero_cross <- RootSpline1(X[round(foot):end], d4[round(foot):end], 
                            verbose = F)
  
  # index of p1
  p1i <- zero_cross[2]
  
  # plot
  # plot(pw[foot:end], type='l', lwd=2, col='grey',yaxt='n',ylab="")
  # par(new=T)
  # plot(d4,type="l",xaxt='n',lwd=2,ylab="4th derivative")
  # abline(h=0,v=p1i,col=2,lwd=1.5)
  
  # Find p2 from 3rd derivative
  p2i <- which.min(d3[maxpi:(notch - 5)]) + maxpi
  
  # Depending type of pressure waveform p1 or p2 will aprox equal max p
  # Find which is closest to max P
  distp1 <- abs(maxpi - p1i)
  distp2 <- abs(maxpi - p2i)
  
  # Which is closer to max P, p1 | p2
  if(distp1 > distp2) {
    p2i <- which.max(pw)
  } else if(distp2 > distp1) {
    p1i <- which.max(pw)
  }
  
  # Calculate augmentation index
  ap <- (pw[p2i] - pw[p1i])
  pp <- (pw[maxpi] - pw[foot])
  aix <- (ap / pp) * 100
  
  # Determine Murgo waveform type
  type <- NA
  if(aix > 12) {
    type <- "A"
  } else  if (aix <= 12 & aix >= 0) {
    type <- "B"
  } else if (aix < 0) {
    type <- "C"
  }
  
  # plot(pw,type="o")
  # abline(v=c(p1i,p2i), col=2:3)
  
  # Find dp/dt max
  dpdt.max <- which.max(d1)
  
  # Find inflection after p1 (local max of 2nd derivative per SphygmoCor)
  p1i2 <- 9999
  if(p1i < maxpi-3) {
    p1i2 <- which.max(d2[p1i:maxpi]) + (p1i - 1)
  }
  
  # Find p1 from 1st derivative (local min as per Kelly et al. 10.1161/01.CIR.80.6.1652)
  # This method was simplified to the 4th derivative method but works fine here
  
  # p1i1 <- 9999
  # if(p1i < maxpi-3) {
  #   p1i1 <- which.min(d1[dpdt.max:p1i2]) + (dpdt.max - 1)
  # }
  
  # plot(pw,type="o"); abline(v=c(p1i, p1i1, p1i2), col=2)
  # par(new=T)
  # plot(d1,type="o"); abline(h=0)
  
  # Plot results
  if (isTRUE(plot)) {
    
    plot(time, pw,
         type = 'l',
         lwd = 3,
         ylab = "Pressure (mmHg)",
         xlab = "Time (s)")
    grid(NULL,NULL, lty = 3, col = "lightgrey") 
    legend("topright", type, bty = 'n')
    
    points(x = c(time[foot], 
                 time[p1i], 
                 time[p1i2], 
                 time[p2i], 
                 time[notch], 
                 time[notchpeak]),
           y = c(pw[foot], 
                 pw[p1i], 
                 pw[p1i2], 
                 pw[p2i], 
                 pw[notch], 
                 pw[notchpeak]),
           pch = "|",
           col = 2,
           lwd = 3,
           cex = 1.7)
    
  }
  
  
  df <- data.frame(
    # Index of values
    MaxP_index = maxpi,
    Foot_index = foot,
    P1_index = p1i,
    P2_index = p2i,
    Ed_index = notch,
    P3_index = notchpeak,
    #P1x_index = p1i1,
    DpDt_index = dpdt.max,
    # Values in unit seconds
    MaxP_sec = time[maxpi],
    Foot_sec = time[foot],
    P1_sec = time[p1i],
    P2_sec = time[p2i],
    Ed_sec = time[notch],
    P3_sec = time[notchpeak],
    #P1x_sec = time[p1i1],
    DpDt_sec = time[dpdt.max],
    # Values in unit mmHg
    MaxP_mmhg = pw[maxpi],
    Foot_mmhg = pw[foot],
    P1_mmhg = pw[p1i],
    P2_mmhg = pw[p2i],
    Ed_mmhg = pw[notch],
    P3_mmhg = pw[notchpeak],
    #P1x_mmhg = pw[p1i1],
    DpDt_mmhg = pw[dpdt.max],
    AP_mmHg = ap,
    AIX = aix,
    Type = type
  )
  
  # round values in df
  df[,1:26] <- round(df[,1:26], 3) # round values in df
  
  # print values to console
  for(i in 1:length(df)){
    print(paste0(names(df[i]),": ", df[1,i]), quote = F)
  }
  
  return(df)
  
}
