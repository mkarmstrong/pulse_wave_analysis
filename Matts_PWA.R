pwa <- function(pw, filt = F, plot = FALSE) {
  
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
    end2 <- length(pw) - round(length(nni:end)/1.7)
    Max_dp2_dbp <- which.max(dp2[nni:end2]) + nni - 2
    Min_dp3_dbp <- which.min(dp3[nni:end2]) + nni
    #resid <- lm(pw[Max_dp2_dbp:Min_dp3_dbp] ~ time(pw[Max_dp2_dbp:Min_dp3_dbp]))$residuals
    
    # De-trend notch area
    narea <- pw[Max_dp2_dbp:Min_dp3_dbp]
    resid <- NA
    for(i in 1:length(narea)) {
      resid[i] <- narea[i] - narea[1]
    }
    
    #plot(resid)
    
    # Find notch
    dic <- unname(which.min(resid)) + Max_dp2_dbp - 1
    
    if(dic >= Min_dp3_dbp-1) {
      #dic <- (length(narea)/2) + Max_dp2_dbp
      dic <- Max_dp2_dbp + 3
    }
    
    # # Testing plots
    # plot(pw);abline(v=c(Max_dp2_dbp, dic, Min_dp3_dbp, end2),
    #                 lty=c(3,1,3,3),
    #                 col=c("grey","red","grey","grey"),
    #                 lwd=2)
    
    
    # FIND DICROTIC PEAK ------------------------------------------------------
    
    above <- dp1[(dic+2):end2] > 0
    dia <- NA
    if(TRUE %in% above) {
      dia <- which(diff(above) < 0) + (dic+2)
      dia <- dia[which.max(pw[dia])]
    } else {
      dia <- which.min(dp2[dic:end2]) + dic
    }
    
    # plot(pw); abline(v=dia)
    # plot(pw[dic:end2])
    # par(new=T)
    # plot(dp1[dic:end2]); abline(h=0)
    # abline(v=dia-(dic))
    
    
    # PLOTS -------------------------------------------------------------------
    
    if(isTRUE(plot)) {
      plot(pw, type = "l", lwd=2, ylab="BP (mmHg)")
      abline(v=c(dic, dia), col="grey", lty=3, lwd=2)
      mtext(c("Ed", "P3"), side = 3, at = c(dic,dia))
    }
    
    return(data.frame(dicrotic_notch = dic, dicrotic_peak = dia))
    
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
    
    ## sort f in case it's passed in backwards
    f <- sort(f)
    
    filt <- signal::butter(
      n = n,
      W = f * 2,
      type = "low",
      plane = "z"
    )
    
    ## remove mean
    yAvg <- mean(y)
    y <- y - yAvg
    
    ## pad the data to twice the max period
    pad <- max(p) * 2
    ny <- length(y)
    ## pad the data
    yPad <- c(y[pad:1], y, y[ny:(ny - pad)])
    ## run the filter
    yFilt <- signal::filtfilt(filt, yPad)
    ## unpad the filtered data
    yFilt <- yFilt[(pad + 1):(ny + pad)]
    ## return with mean added back in
    filt.sig <- yFilt + yAvg
    
    if(isTRUE(do.plot)){
      ## plot results
      plot(filt.sig,
           type = "l",
           lwd = 2)
    }
    
    ## return filtered signal
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
  
  # Create time
  time <- (0:(length(pw)-1)) / 200
  
  # Some additional calcs
  end <- length(pw)
  foot <- Tintersect(pw) - 1
  if(foot < 1) {foot <- 1}
  maxpi <- which.max(pw)
  notchdat <- dicrotic(pw)
  notch <- notchdat$dicrotic_notch
  notchpeak <- notchdat$dicrotic_peak
  
  # plot(pw,type="o"): abline(v=foot)
  
  # Find p1 from 4th derivative
  above <- d4[foot:end] > 0
  above[1:10] <- FALSE
  n1 <- which(diff(above) < 0)
  n2 <- n1[1]
  p1i <- n2+(foot-1)
  
  # plot(pw,type="l")
  # par(new=T)
  # plot(d4,type="l",col=3);abline(h=0,v=p1i,col=2)
  # par(new=T)
  # plot(d1,type='l')
  
  # Find p2 from 3rd derivative
  p2i <- which.min(d3[maxpi:(notch - 5)]) + maxpi
  
  # # Find p2 from 4th derivative (P3d is better method)
  # below <- d4[maxpi:notch] < 0
  # blw1 <- which(diff(below) < 0) + 1
  # blw2 <- blw1[1]
  # p2i <- blw2 + (maxpi)
  # 
  # plot(pw[maxpi:(notch)],type="l"); abline(v = c(p2i - maxpi, p2i.3rd - maxpi))
  # par(new = T)
  # plot(d3[maxpi:(notch)],type="o", col=2)
  # par(new = T)
  # plot(d4[maxpi:(notch)],type="o", col=4); abline(h=0,col=4)
  
  # Depending type of pressure waveform p1 or p2 will equal max p
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
  
  # Find inflection after p1 (local max of 2nd derivative)
  p1i2 <- 9999
  if(p1i < maxpi-3) {
    p1i2 <- which.max(d2[p1i:maxpi]) + (p1i - 1)
  }
  
  # Find p1 from 1st D (local min of 1st derivative. Kelly et al. 10.1161/01.CIR.80.6.1652)
  p1i1 <- 9999
  if(p1i < maxpi-3) {
    p1i1 <- which.min(d1[dpdt.max:p1i2]) + (dpdt.max - 1)
  }
  
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
    grid (NULL,NULL, lty = 3, col = "lightgrey") 
    legend("topright", type, bty = 'n')
    
    # mtext(
    #   c("Ft", "P1", "P2", "Ed", "P3"),
    #   side = 3,
    #   cex = .7,
    #   at = c(time[foot],
    #          time[p1i],
    #          time[p2i],
    #          time[notch],
    #          time[notchpeak])
    # )
    
    points(x = c(time[foot], 
                 time[p1i], 
                 #time[p1i1], 
                 time[p2i], 
                 time[notch], 
                 time[notchpeak]),
           y = c(pw[foot], 
                 pw[p1i], 
                 #pw[p1i1], 
                 pw[p2i], 
                 pw[notch], 
                 pw[notchpeak]),
           pch = "|",
           col = 2,
           lwd = 3,
           cex = 1.7)
    
    
    points(x = time[p1i1],
           y = pw[p1i1],
           pch = "|",
           col = 4,
           lwd = 3,
           cex = 1.7)
    
    
  }
  
  
  df <- data.frame(
    # Index of values
    MaxP.index = maxpi,
    Foot.index = foot,
    P1.index = p1i,
    P2.index = p2i,
    Ed.index = notch,
    P3.index = notchpeak,
    P1x.index = p1i1,
    DpDt.index = dpdt.max,
    # Values in unit seconds
    MaxP.sec = time[maxpi],
    Foot.sec = time[foot],
    P1.sec = time[p1i],
    P2.sec = time[p2i],
    Ed.sec = time[notch],
    P3.sec = time[notchpeak],
    P1x.sec = time[p1i1],
    DpDt.sec = time[dpdt.max],
    # Values in unit mmHg
    MaxP.mmhg = pw[maxpi],
    Foot.mmhg = pw[foot],
    P1.mmhg = pw[p1i],
    P2.mmhg = pw[p2i],
    Ed.mmhg = pw[notch],
    P3.mmhg = pw[notchpeak],
    P1x.mmhg = pw[p1i1],
    DpDt.mmhg = pw[dpdt.max],
    AP.mmHg = ap,
    AIX = aix,
    # Err check
    Type = type
  )
  
  
  return(df)
  
}
