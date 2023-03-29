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
    pw <- low.pass(pw, 0.10, do.plot = F)
  }
  
  # Create derivatives
  d1 <- fsg721(pw)
  d2 <- fsg721(d1)
  d3 <- fsg721(d2)
  d4 <- fsg721(d3)
  
  time <- (0:(length(pw)-1)) / 200
  
  # Some additional calcs
  foot <- Tintersect(pw)
  if(foot < 1) {foot <- 1}
  dpdt.max <- which.max(d1) # Find dp/dt max
  maxpi <- which.max(pw)
  notchdat <- dicrotic(pw)
  notch <- notchdat$dicrotic_notch
  notchpeak <- notchdat$dicrotic_peak
  end <- length(pw)
  
  
  # Find P1 concave ---------------------------------------------------------
  
  # Find inflection after convex p1 shoulder (local max of 2nd derivative per SphygmoCor)
  loc_hold <- which.min(d2[dpdt.max:maxpi]) + (dpdt.max - 1)
  p1i_alt <- which.max(d2[loc_hold:maxpi]) + (loc_hold - 1)
  #p1i_alt <- which.max(d2[(dpdt.max+2):maxpi]) + ((dpdt.max+2) - 1)
  
  # plot(pw,type="b")
  # abline(v=p1i_alt, col=3, lwd=2)
  # par(new=T)
  # plot(d2,type="o",col="grey")
  # abline(h=0, v=c(loc_hold, dpdt.max+2, maxpi), lty=c(1,3,1, 1))
  
  
  # Find P2 -----------------------------------------------------------------
  
  # Find p2 from 3rd derivative
  # when looking for p2 after sbp, the 3rd derivative method is more robust than
  # the 4th derivative method.
  p2i <- which.min(d3[maxpi:(notch - 5)]) + maxpi - 1
  
  # plot(pw)
  # par(new=T)
  # plot(d3, type="o", col="grey")
  # abline(h=0, v=c(maxpi, notch-5, p2i), col=c(1,1,1,2))
  
  
  # Determine waveform type -------------------------------------------------
  
  # In this version of pwa(), waveform type is determined from the value of d2 
  # at p1i_alt. If d2 is >= 0, then the early systolic p1 inflection was present.
  # A buffer zone was added following some miss classifications. 
  buff <- 0.15 * (abs(max(d2)) - abs(min(d2))) # 10% of max of d2
  
  if(d2[p1i_alt]+buff >= 0) {
    p2i <- maxpi
  } else {
    p1i_alt <- maxpi
    p1i <- maxpi
  }
  
  # plot(pw)
  # abline(v=c(p1i_alt, p2i), col=2:3, lwd=1.5)
  
  # Calculate augmentation index
  ap  <- (pw[p2i] - pw[p1i_alt])
  pp  <- (pw[maxpi] - pw[foot])
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
  
  # Find convex P1 ----------------------------------------------------------
  
  # Find p1 from 0 crossing of 4th derivative, most common in the literature
  # or using the first derivative
  
  if(type == "A" | type == "B") {
    
    p1_seg <- pw[(foot+5):p1i_alt]
    xidx   <- 1:end
    
    # use d1 for very definable p1 & d4 for subtle p1
    if(max(p1_seg) > pw[p1i_alt]) {
      
      zero_cross <- RootSpline1(xidx[(foot+5):p1i_alt], 
                                d1[(foot+5):p1i_alt], 
                                verbose = F) 
      # crossing from above to below
      p1i_4d <- zero_cross[d1[zero_cross] > 0]
      
    } 
    else {
      
      zero_cross <- RootSpline1(xidx[(foot+5):p1i_alt], 
                                d4[(foot+5):p1i_alt], 
                                verbose = F)
      # crossing from above to below
      p1i_4d <- zero_cross[d4[zero_cross] > 0]
      
    }
    
    
    # # as safety checks, find p1 using local max of 3rd d
    # p1i_3d <- which.max(d3[(foot + 5):p1i_alt]) + (foot + 5) - 1
    # 
    # # error trap if p1 is not found
    # if(length(p1i_4d) < 1) {
    #   # if p1 not found use p1 from 3rd D
    #   p1i <- p1i_3d
    # } else {
    #   # if multiple zero crossing in 4th (or 1st) d, find crossing closest to p1 from 3rd d 
    #   p1i <- p1i_4d[which.min(abs(p1i_4d - p1i_3d))]
    # }
    
    
    # error trap if p1 is not found
    if(length(p1i_4d) < 1) {
      # if p1 not found use p1 from 3rd D
      p1i <- which.max(d3[(foot + 5):p1i_alt]) + (foot + 5) - 1
    } else {
      # if found, take zero crossing closest to p1_alt
      p1i <- p1i_4d[which.min(abs(p1i_4d - p1i_alt))]
    }
    
  }
  
  
  # plot(p1_seg, pch=19,
  #      xlim=c(1,80), col=5)
  # 
  # par(new=T)
  # plot(d1[(foot+5):p1i_alt],
  #      xlim=c(1,80),
  #      type='l',
  #      col=1)
  # abline(h=0, col=1)
  # 
  # par(new=T)
  # plot(d3[(foot+5):p1i_alt],
  #      xlim=c(1,80),
  #      type='l',
  #      col=3)
  # abline(h=0, v=p1i_3d-foot-5, col=3)
  # 
  # par(new=T)
  # plot(d4[(foot+5):p1i_alt],
  #      xlim=c(1,80),
  #      type='l',
  #      col=4)
  # abline(h=0,col=4)
  
  
  # Plot --------------------------------------------------------------------
  
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
      
      points(x = time[p1i_alt],
             y = pw[p1i_alt],
             pch = "|",
             col = 5,
             lwd = 3,
             cex = 1.7)
      
    }
    
  }
  
  
  df <- data.frame(
    # Index of values
    sp_idx = maxpi,
    dp_idx = foot,
    p1_idx = p1i,
    p2_idx = p2i,
    es_idx = notch,
    p3_idx = notchpeak,
    p1_alt_idx = p1i_alt,
    dpdt_idx = dpdt.max,
    # Values in unit seconds
    sp_sec = time[maxpi],
    dp_sec = time[foot],
    p1_sec = time[p1i],
    p2_sec = time[p2i],
    es_sec = time[notch],
    p3_sec = time[notchpeak],
    p1_alt_sec = time[p1i_alt],
    dpdt_sec = time[dpdt.max],
    # Values in unit mmHg
    sp_mmhg = pw[maxpi],
    dp_mmhg = pw[foot],
    p1_mmhg = pw[p1i],
    p2_mmhg = pw[p2i],
    es_mmhg = pw[notch],
    p3_mmhg = pw[notchpeak],
    p1_alt_mmhg = pw[p1i_alt],
    dpdt_mmhg = pw[dpdt.max],
    ap_mmHg = ap,
    aix = aix,
    type
  )
  
  # round values in df
  num_cols <- unlist(lapply(df, is.numeric)) # Identify numeric cols
  df[num_cols] <-  round(df[num_cols], 3)    # round numeric cols
  
  # print values to console
  # for(i in 1:length(df)){
  #   print(paste0(names(df[i]),": ", df[1,i]), quote = F)
  # }
  
  return(df)
  
}
