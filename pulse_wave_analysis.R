# LOAD REQUIRED FUNCTIONS
install.packages("signal")

fsg721 <- function(x, smth = 7) {
  
  # 2nd order polynomial sg filter
  # windowing for smoothing = smth (default 7)
  sg <- signal::sgolay(p = 2, n = smth, m = 1)
  sig <- signal::filter(sg, x)
  #plot(sig)
  return(sig)
  
}

low.pass <- function(y, fq, do.plot = FALSE) {
  
  # Second order low pass filter
  # Removes high frequency components above fq
  # y = a numeric vector
  # fq = a numeric vector giving frequency or period of the filter.
  
  if (any(is.na(y))) stop("y contains NA")
  
  # n = a numeric value giving the order of the filter. 
  # Larger numbers create steeper fall off.
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
  
  filt <- signal::butter(n = n,
                         W = f * 2,
                         type = "low",
                         plane = "z")
  
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
  
  # End index
  end <- length(pw)
  
  # End index without potential perturbation at end diastole  
  end2 <- end * .9
  
  # Isolate notch area with 1st derivatives
  nni <- which.min(dp1)
  
  # Dicrotic notch from local dp2 max
  dic <- which.max(dp2[nni:end2]) + nni - 1
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = dic, h = 0)
  
  
  # FIND DICROTIC PEAK ------------------------------------------------------
  
  end3 <- ((end - dic) * .6) + dic # 60% of diastolic duration
  
  # Dicrotic peak from min of 2nd derivative
  # works better for subtle peaks
  if(sum(dp2[dic:end3] < 0) < 1) {
    dia <- 9999
  } else {
    dia <- which.min(dp2[dic:end3]) + dic - 1
  }
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = c(dic, dia), h = 0)
  
  
  # Dicrotic peak from 0 crossing of 1st derivative
  # works better for very definable peaks
  if (pw[dia] > pw[dic] & !is.na(pw[dia])) {
    hold <- RootSpline1(1:(end - nni + 1), dp1[nni:end], verbose = F)
    dia <- hold[2] + nni - 1
  }
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp1, type='o',col="grey")
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


pwa <- function(pw, filt = FALSE, plot = FALSE) {
  
  # Low pass waveform
  if (isTRUE(filt)) {
    pw <- low.pass(pw, 10, do.plot = F)
  }
  
  # Create derivatives
  d1 <- fsg721(pw)
  d2 <- fsg721(d1)
  d3 <- fsg721(d2)
  d4 <- fsg721(d3)
  
  # Some additional calcs
  end <- length(pw)
  foot <- Tintersect(pw)
  if(foot < 1) {foot <- 1}
  maxpi <- which.max(pw)
  notchdat <- dicrotic(pw)
  notch <- notchdat$dicrotic_notch
  notchpeak <- notchdat$dicrotic_peak
  time <- (0:(length(pw)-1)) / 200
  dpdt.max <- which.max(d1) # Find dp/dt max
  
  # plot(pw,type="o")
  # abline(v=c(foot, notch, notchpeak))
  
  
  # Find p1 up inflection ---------------------------------------------------------
  
  # TODO find a better marker than "dpdt.max+5"
  
  # Find inflection after p1 shoulder (local max of 2nd derivative per SphygmoCor)
  loc_hold <- which.min(d2[dpdt.max:maxpi]) + (dpdt.max - 1)
  p1i2 <- which.max(d2[loc_hold:maxpi]) + (loc_hold - 1)
  #p1i2 <- which.max(d2[which.min(d2):maxpi]) + (which.min(d2) - 1)
  
  # plot(pw,type="o")
  # abline(v=p1i2, col=3)
  # par(new=T)
  # plot(d2,type="o",col="grey"); abline(h=0, v=loc_hold)
  

  # Find p1 shoulder ---------------------------------------------------------
  
  # Find p1 from 0 crossing of 4th derivative, most common method
  
  inc <- pw[(foot+5):p1i2]
  
  # plot(inc)
  # par(new=T)
  # plot(d4[(foot+5):p1i2],
  #      type='l',
  #      col=2)
  # abline(h=0,col=2)
  # par(new=T)
  # plot(d3[(foot+5):p1i2],
  #      type='o',
  #      col=3)
  
  Xindex <- 1:length(pw)
  
  zero_cross <- RootSpline1(Xindex[(foot+5):p1i2], 
                            d4[(foot+5):p1i2], 
                            verbose = F) 
  
  p1i_4d <- zero_cross[d4[zero_cross] > 0]
  p1i_3d <- which.max(d3[(foot+5):p1i2]) + (foot+5) - 1
  
  # plot(pw)
  # abline(v=c(p1i_4d,p1i_3d), col=2:3)
  
  p1i <- p1i_4d[which.min(abs(p1i_4d - p1i_3d))]
  
  
  # -------------------------------------------------------------------------
  
  
  # Find p2 from 3rd derivative
  # when looking for p2 after sbp, the 3rd derivative method is more robust than
  # the 4th derivative method, per what ive experienced, though the derivative filter is important here.
  p2i <- which.min(d3[maxpi:(notch - 5)]) + maxpi - 1
  
  # plot(pw)
  # par(new=T)
  # plot(d3,type="l",xaxt='n',lwd=2,ylab="4th derivative")
  # abline(h=0,v=p2i,col=2,lwd=1.5)

  # Depending on type of pressure waveform p1 or p2 will approx equal max p
  # Find which is closest to max P
  distp1 <- abs(maxpi - p1i)
  distp2 <- abs(maxpi - p2i)
  
  # TODO find better solution to the if(length(distp1) < 1) problem below
  
  # major hackyness, eww
  if(length(distp1) < 1) {
    distp1 <- distp2 - 1
  }
  
  # Which is closer to max P, p1 | p2
  if(distp1 > distp2) {        # if p2 is closer then...
    p2i <- maxpi
  } else if(distp2 > distp1) { # if p1 is closer then...
    p1i <- maxpi
  }
  
  # plot(pw)
  # abline(v=c(p1i,p2i),col=2:3,lwd=1.5)
  
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
                 time[p2i], 
                 time[notch], 
                 time[notchpeak]),
           y = c(pw[foot], 
                 pw[p1i], 
                 pw[p2i], 
                 pw[notch], 
                 pw[notchpeak]),
           pch = "|",
           col = 2,
           lwd = 3,
           cex = 1.7)
    
    if(aix > 0) {
      
      points(x = time[p1i2],
             y = pw[p1i2],
             pch = "|",
             col = 5,
             lwd = 3,
             cex = 1.7)
      
    }
    
  }
  
  
  df <- data.frame(
    # Index of values
    MaxP_index = maxpi,
    Foot_index = foot,
    P1_index = p1i,
    P2_index = p2i,
    Ed_index = notch,
    P3_index = notchpeak,
    P1_inflect_index = p1i2,
    DpDt_index = dpdt.max,
    # Values in unit seconds
    MaxP_sec = time[maxpi],
    Foot_sec = time[foot],
    P1_sec = time[p1i],
    P2_sec = time[p2i],
    Ed_sec = time[notch],
    P3_sec = time[notchpeak],
    P1_inflect_sec = time[p1i2],
    DpDt_sec = time[dpdt.max],
    # Values in unit mmHg
    MaxP_mmhg = pw[maxpi],
    Foot_mmhg = pw[foot],
    P1_mmhg = pw[p1i],
    P2_mmhg = pw[p2i],
    Ed_mmhg = pw[notch],
    P3_mmhg = pw[notchpeak],
    P1_inflect_mmhg = pw[p1i2],
    DpDt_mmhg = pw[dpdt.max],
    AP_mmHg = ap,
    AIX = aix,
    Type = type
  )
  
  # round values in df
  num_cols <- unlist(lapply(df, is.numeric)) # Identify numeric cols
  df[num_cols] <-  round(df[num_cols], 3)    # round numeric cols
  
  # print values to console
  for(i in 1:length(df)){
    print(paste0(names(df[i]),": ", df[1,i]), quote = F)
  }
  
  return(df)
  
}
