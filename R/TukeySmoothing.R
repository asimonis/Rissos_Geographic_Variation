spNLFilter <- function(Signal, Type, ...) {
  if (!is.character(Type) || length(Type) != 1)
    stop("Type must be a single string")
  
  args <- list(...)
  
  # --- running median with shrinking windows (Matlab fallback behavior) ---
  median_shrink <- function(x, k) {
    if (k %% 2 == 0) stop("Median filter size must be odd")
    n <- length(x)
    if (n == 0) return(x)
    mid <- (k - 1) / 2
    out <- x
    
    # start
    for (i in seq_len(mid)) {
      out[i] <- stats::median(x[1:(i + mid)])
    }
    # middle
    for (i in (mid + 1):(n - mid)) {
      out[i] <- stats::median(x[(i - mid):(i + mid)])
    }
    # end
    for (i in (n - mid + 1):n) {
      out[i] <- stats::median(x[(i - mid):n])
    }
    out
  }
  
  # --- Hanning (Hann) convolution + trim (Matlab conv + phase dump) ---
  hanning_conv_trim <- function(x, filter_size) {
    if (filter_size %% 2 == 0) stop("Hanning window size must be odd")
    n <- length(x)
    if (n == 0) return(x)
    
    j <- 0:(filter_size - 1)
    w <- 0.5 - 0.5 * cos(2 * pi * j / (filter_size - 1))
    w <- w / sum(w)
    
    y_full <- as.numeric(stats::convolve(x, rev(w), type = "open"))
    offset <- floor(filter_size / 2)
    y_full[(offset + 1):(offset + n)]
  }
  
  t <- tolower(Type)
  
  if (t %in% c("median", "m")) {
    if (length(args) < 1) stop("Median requires FilterSize")
    return(median_shrink(as.numeric(Signal), as.integer(args[[1]])))
  }
  
  if (t %in% c("hanning", "hann", "h")) {
    if (length(args) < 1) stop("Hanning requires FilterSize")
    return(hanning_conv_trim(as.numeric(Signal), as.integer(args[[1]])))
  }
  
  if (t %in% c("tukey", "t")) {
    if (length(args) < 2) stop("Tukey requires Median and HanningSize")
    Median <- as.integer(args[[1]])
    HanningSize <- as.integer(args[[2]])
    
    # estimate
    EstSignal <- spNLFilter(Signal, "Median", Median)
    EstSignal <- spNLFilter(EstSignal, "Hanning", HanningSize)
    
    # residual
    Error <- as.numeric(Signal) - EstSignal
    Error <- spNLFilter(Error, "Median", Median)
    Error <- spNLFilter(Error, "Hanning", HanningSize)
    
    return(EstSignal + Error)
  }
  
  stop(sprintf('Unrecognized filter type "%s"', Type))
}



calculateAverageSpectra_tukey <- function(
    x, evNum = 1, calibration = NULL, wl = 512,
    channel = 1:2, filterfrom_khz = 0, filterto_khz = NULL,
    sr = NULL, snr = 0, norm = TRUE, plot = FALSE, noise = FALSE,
    decimate = 1, sort = FALSE, mode = "spec",
    title = NULL, ylim = NULL, flim = NULL,
    # ---- new args ----
    tukeySmooth = TRUE,
    tukeyMedian = 5L,
    tukeyHanningSize = 5L,
    smoothNoise = FALSE
) {
  if (!requireNamespace("PAMpal", quietly = TRUE)) {
    stop("Package 'PAMpal' is required.")
  }
  if (!exists("spNLFilter", mode = "function")) {
    stop("I can't find spNLFilter() in your session. Source/define it first.")
  }
  
  # --- helper: PAMpal's log-average behavior (mirrors doLogAvg usage) ---
  doLogAvg_local <- function(mat, log = TRUE) {
    if (is.null(mat) || length(mat) == 0) return(NULL)
    mat <- as.matrix(mat)
    if (!log) {
      return(apply(mat, 1, mean, na.rm = TRUE))
    }
    # spectra are in dB => average in linear domain then convert back
    20 * log10(apply(mat, 1, function(v) mean(10^(v / 20), na.rm = TRUE)))
  }
  
  # --- helper: smooth each click spectrum (columns are clicks in PAMpal) ---
  smooth_specMat_cols <- function(specMat, Median, HanningSize) {
    specMat <- as.matrix(specMat)
    dn <- dimnames(specMat)
    out <- apply(specMat, 2, function(spec) {
      spec <- as.numeric(spec)
      if (length(spec) < max(Median, HanningSize) || all(is.na(spec))) return(spec)
      spNLFilter(spec, "Tukey", Median, HanningSize)
    })
    out <- as.matrix(out)
    dimnames(out) <- dn
    out
  }
  
  # 1) Get RAW per-click spectra/noise from PAMpal, without SNR filtering, no plotting, no normalization
  raw <- PAMpal::calculateAverageSpectra(
    x = x, evNum = evNum, calibration = calibration, wl = wl,
    channel = channel, filterfrom_khz = filterfrom_khz, filterto_khz = filterto_khz,
    sr = sr, snr = 0, norm = FALSE, plot = FALSE, noise = TRUE,
    decimate = decimate, sort = FALSE, mode = mode,
    title = title, ylim = ylim, flim = flim
  )
  if (is.null(raw)) return(NULL)
  
  freq     <- raw$freq
  clickUID <- raw$UID
  specMat_raw  <- raw$allSpec   # freq x click
  noiseMat_raw <- raw$allNoise  # freq x click
  
  # 2) Compute snrKeep using RAW matrices (PAMpal logic)
  if (snr > 0) {
    snrVals <- rep(NA_real_, ncol(specMat_raw))
    for (i in seq_along(snrVals)) {
      minX <- switch(mode,
                     "spec" = 1,
                     "ceps" = min(which(specMat_raw[-(1:2), i] < 3 * stats::median(specMat_raw[-(1:2), i]))) + 2
      )
      wherePeak <- which.max(specMat_raw[minX:nrow(specMat_raw), i]) + minX - 1
      if (length(wherePeak) == 0) next
      if (is.na(noiseMat_raw[wherePeak, i])) {
        snrVals[i] <- Inf
      } else {
        snrVals[i] <- specMat_raw[wherePeak, i] - noiseMat_raw[wherePeak, i]
      }
    }
    
    if (!any(snrVals >= snr, na.rm = TRUE)) {
      warning("No clicks above SNR threshold.", call. = FALSE)
      return(NULL)
    }
    snrKeep <- snrVals >= snr
  } else {
    snrKeep <- rep(TRUE, ncol(specMat_raw))
  }
  
  # 3) Smooth spectra (and optionally noise) AFTER snrKeep is decided
  specMat_use <- specMat_raw
  noiseMat_use <- noiseMat_raw
  
  if (isTRUE(tukeySmooth)) {
    specMat_use <- smooth_specMat_cols(specMat_raw, tukeyMedian, tukeyHanningSize)
    if (isTRUE(smoothNoise)) {
      noiseMat_use <- smooth_specMat_cols(noiseMat_raw, tukeyMedian, tukeyHanningSize)
    }
  }
  
  # 4) Average (log-avg for spec, linear avg for ceps) on *kept clicks*
  avgNoise <- doLogAvg_local(noiseMat_use[, snrKeep, drop = FALSE], log = (mode == "spec"))
  avgSpec  <- doLogAvg_local(specMat_use[,  snrKeep, drop = FALSE], log = (mode == "spec"))
  
  if (isTRUE(norm) && !is.null(avgSpec)) {
    maxVal <- max(avgSpec, na.rm = TRUE)
    avgSpec  <- avgSpec - maxVal
    avgNoise <- avgNoise - maxVal
  }
  
  # Optional minimal plotting (average spectrum only)
  if (isTRUE(plot)) {
    xlab <- if (mode == "spec") "Frequency (kHz)" else "ICI (ms)"
    xf <- if (mode == "spec") freq / 1e3 else freq * 1e3
    graphics::plot(xf, avgSpec, type = "l", xlab = xlab, ylab = if (mode=="spec") "dB" else "Value",
                   main = if (is.null(title)) "Average (Tukey-smoothed) spectrum" else title)
    if (isTRUE(noise)) graphics::lines(xf, avgNoise, lty = 2)
    if (isTRUE(noise)) graphics::legend("topright", legend = c("Signal", "Noise"), lty = c(1,2), bty = "n")
  }
  
  # ---- concatenated spectrum plot payload ----
  concatPlot <- local({
    x_khz <- freq / 1000
    y_signal <- avgSpec
    y_noise <- if (!is.null(avgNoise)) avgNoise else NULL
    
    list(
      x_khz = x_khz,
      y_signal = y_signal,
      y_noise = y_noise,
      plot = function(main = "Concatenated spectrum", showNoise = TRUE, ...) {
        graphics::plot(
          x_khz, y_signal, type = "l",
          xlab = "Frequency (kHz)", ylab = "Level (dB)",
          main = main, ...
        )
        if (showNoise && !is.null(y_noise)) {
          graphics::lines(x_khz, y_noise, lty = 2)
          graphics::legend(
            "topright",
            legend = c("Signal", "Noise"),
            lty = c(1, 2),
            bty = "n"
          )
        }
      }
    )
  })
  
  
  invisible(list(
    freq = freq,
    UID = clickUID,
    avgSpec = avgSpec,
    allSpec = specMat_use,
    avgNoise = avgNoise,
    allNoise = noiseMat_use,
    snrKeep = snrKeep,  # extra: which clicks were kept by RAW SNR
    concatPlot=concatPlot
  ))
}

truncateAverageSpectra <- function(res, fmin_khz = 19, fmax_khz = 50) {
  if (is.null(res) || is.null(res$freq)) stop("Input doesn't look like calculateAverageSpectra() output")
  
  fmin_hz <- fmin_khz * 1000
  fmax_hz <- fmax_khz * 1000
  
  keep <- res$freq >= fmin_hz & res$freq <= fmax_hz
  if (!any(keep)) stop("No frequency bins remain after truncation.")
  
  res$freq <- res$freq[keep]
  
  # vectors (same length as freq)
  if (!is.null(res$avgSpec))  res$avgSpec  <- res$avgSpec[keep]
  if (!is.null(res$avgNoise)) res$avgNoise <- res$avgNoise[keep]
  
  # matrices (rows correspond to freq)
  if (!is.null(res$allSpec))  res$allSpec  <- res$allSpec[keep, , drop = FALSE]
  if (!is.null(res$allNoise)) res$allNoise <- res$allNoise[keep, , drop = FALSE]
  
  res
}


plotAverageSpectrum <- function(res, showNoise = TRUE, main = NULL,
                                xlab = "Frequency (kHz)", ylab = "Level (dB)",
                                colSpec = "black", colNoise = "grey50") {
  if (is.null(res$freq) || is.null(res$avgSpec))
    stop("Object doesn't look like calculateAverageSpectra output")
  
  x <- res$freq / 1000  # kHz
  
  graphics::plot(x, res$avgSpec, type = "l",
                 xlab = xlab, ylab = ylab,
                 main = main, col = colSpec)
  
  if (showNoise && !is.null(res$avgNoise)) {
    graphics::lines(x, res$avgNoise, lty = 2, col = colNoise)
    graphics::legend("topright",
                     legend = c("Signal", "Noise"),
                     lty = c(1, 2),
                     col = c(colSpec, colNoise),
                     bty = "n")
  }
}


plotCompareSpectra <- function(raw, smooth,
                               fmin_khz = 19, fmax_khz = 50,
                               main = "Raw vs Tukey-smoothed") {
  raw   <- truncateAverageSpectra(raw, fmin_khz, fmax_khz)
  smooth <- truncateAverageSpectra(smooth, fmin_khz, fmax_khz)
  
  x <- raw$freq / 1000
  
  graphics::plot(x, raw$avgSpec, type = "l", col = "grey60",
                 xlab = "Frequency (kHz)", ylab = "Level (dB)",
                 main = main)
  graphics::lines(x, smooth$avgSpec, col = "black", lwd = 2)
  
  graphics::legend("topright",
                   legend = c("Raw", "Tukey"),
                   col = c("grey60", "black"),
                   lwd = c(1, 2),
                   bty = "n")
}



#5 point Hann smoother
hann5 <- function(x) {
  w <- c(0.25, 0.5, 1, 0.5, 0.25)
  w <- w / sum(w)
  as.numeric(stats::filter(x, w, sides = 2))
}

#Plot a single click by UID, then step through clicks interactively (press Enter to advance)
plot_one_click <- function(raw, tukey, uid,
                           fmin_khz = 19, fmax_khz = 50,
                           tukeyMedian = 3, tukeyHanningSize = 5) {
  
  i_raw <- match(uid, raw$UID)
  i_tuk <- match(uid, tukey$UID)
  
  keep <- raw$freq >= fmin_khz*1000 & raw$freq <= fmax_khz*1000
  x_khz <- raw$freq[keep] / 1000
  
  y_raw <- raw$allSpec[keep, i_raw]
  y_tuk <- tukey$allSpec[keep, i_tuk]
  y_h5  <- hann5(y_raw)
  
  ylim <- range(c(y_raw, y_tuk, y_h5), na.rm = TRUE)
  
  plot(x_khz, y_raw, type="l", ylim=ylim,
       xlab="Frequency (kHz)", ylab="Level (dB)",
       main=paste("UID:", uid))
  
  lines(x_khz, y_tuk, lwd=2)
  lines(x_khz, y_h5, lty=2, lwd=2)
  
  legend("topright",
         legend=c("Raw", "Tukey (double)", "Hann-5"),
         lwd=c(1,2,2), lty=c(1,1,2), bty="n")
}



browse_clicks <- function(raw,
                          uids = raw$UID,
                          fmin_khz = 19, fmax_khz = 50,
                          tukeyMedian = 3L,
                          tukeyHanningSize = 5L) {
  
  # basic checks
  if (is.null(raw$freq) || is.null(raw$allSpec) || is.null(raw$UID))
    stop("raw must contain $freq, $allSpec, and $UID")
  
  if (!exists("spNLFilter", mode = "function"))
    stop("spNLFilter() not found. Source/define it first.")
  
  if (!exists("hann5", mode = "function"))
    stop("hann5() not found. Define it first.")
  
  # keep only valid uids
  uids <- intersect(uids, raw$UID)
  if (length(uids) == 0) stop("No UIDs to browse.")
  
  # frequency band
  keep_f <- raw$freq >= fmin_khz * 1000 & raw$freq <= fmax_khz * 1000
  if (!any(keep_f)) stop("No frequency bins in the requested range.")
  x_khz <- raw$freq[keep_f] / 1000
  
  for (uid in uids) {
    i <- match(uid, raw$UID)
    
    y_raw <- raw$allSpec[keep_f, i]
    y_tuk <- spNLFilter(y_raw, "Tukey", tukeyMedian, tukeyHanningSize)
    y_h5  <- hann5(y_raw)
    
    ylim <- range(c(y_raw, y_tuk, y_h5), na.rm = TRUE)
    
    graphics::plot(x_khz, y_raw, type = "l", ylim = ylim,
                   xlab = "Frequency (kHz)", ylab = "Level (dB)",
                   main = paste("UID:", uid))
    
    graphics::lines(x_khz, y_tuk, lwd = 2, lty = 1)
    graphics::lines(x_khz, y_h5,  lwd = 2, lty = 2)
    
    graphics::legend("topright",
                     legend = c("Raw", "Tukey (double)", "Hann-5"),
                     lty = c(1, 1, 2),
                     lwd = c(1, 2, 2),
                     bty = "n")
    
    ans <- readline("Enter = next, q = quit: ")
    if (tolower(ans) == "q") break
  }
  
  invisible(TRUE)
}

