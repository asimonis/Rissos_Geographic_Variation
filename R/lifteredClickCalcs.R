lifteredClickCalcs <- function (data, sr_hz = "auto", calibration = NULL,
                                filterfrom_khz = 10,
                                filterto_khz = 80, winLen_sec = 0.0025)
{
  result <- list()
  paramNames <- c("noiseLevel", "duration",
                  "peakTime", "peak", "peak2",
                  "peak3", "trough", "trough2",
                  "peakToPeak2", "peakToPeak3",
                  "peak2ToPeak3", "dBPP", "Q_10dB",
                  "fmin_10dB","fmax_10dB", "BW_10dB",
                  "centerkHz_10dB", "Q_3dB",
                  "fmin_3dB","fmax_3dB", "BW_3dB",
                  "centerkHz_3dB",
                  "wave")  # <-- NEW base name, becomes liftered_wave after rename
  
  ## ---- hpLifterClick helper -----------------
  hpLifterClick <- function(data,
                            calibration = NULL,
                            n_remove    = 6,
                            return_wave = TRUE,
                            chan        = 1,
                            ...) {
    
    if (is.null(data$wave)) {
      stop("data$wave is NULL – hpLifterClick expects data$wave to be a matrix.")
    }
    
    wave <- data$wave
    if (!is.matrix(wave)) wave <- as.matrix(wave)
    
    n_chan <- ncol(wave)
    if (n_chan < 1) stop("data$wave has zero channels.")
    if (chan < 1 || chan > n_chan) {
      stop("Requested chan = ", chan, " but data$wave has only ", n_chan, " channel(s).")
    }
    
    result <- data.frame(ChannelIdx = seq_len(n_chan))
    
    if (return_wave) {
      result$liftered_wave <- replicate(n_chan, NA, simplify = FALSE)
      
      x <- as.numeric(wave[, chan])
      n <- length(x)
      
      spec    <- fft(x)
      logSpec <- log(Mod(spec) + .Machine$double.eps)
      cep     <- Re(fft(logSpec, inverse = TRUE)) / n
      
      if (n_remove > 0 && n_remove < length(cep)) {
        cep[1:n_remove] <- 0
      }
      
      logSpec_hp <- Re(fft(cep))
      spec_hp    <- exp(logSpec_hp) * exp(1i * Arg(spec))
      x_hp       <- Re(fft(spec_hp, inverse = TRUE)) / n
      
      result$liftered_wave[[chan]] <- x_hp
    }
    
    return(result)
  }
  ## -------------------------------------------------------------------
  
  if (inherits(data, "Wave")) {
    data <- WaveMC(data)
    data <- list(wave = data@.Data, sr = data@samp.rate)
  }
  if (!is.matrix(data$wave)) {
    data$wave <- matrix(data$wave, ncol = 1)
  }
  
  neededVals <- c("wave")
  if (sr_hz == "auto") neededVals <- c(neededVals, "sr")
  missingVals <- neededVals[!(neededVals %in% names(data))]
  if (length(missingVals) > 0) {
    pamWarning("Values for ", paste(missingVals, collapse = ", "),
               " are missing.", "These are required for Click Calculations, please fix.")
    return(NULL)
  }
  
  for (chan in 1:ncol(data$wave)) {
    
    thisWave <- data$wave[, chan]
    
    # --- apply hp liftering to this channel
    lifter_res <- hpLifterClick(
      data = list(wave = as.matrix(thisWave)),
      chan = 1,
      return_wave = TRUE
    )
    thisWave <- lifter_res$liftered_wave[[1]]
    # ---------------------------------------------------------------
    
    if (all(thisWave == 0)) {
      blanks <- data.frame(matrix(NA, nrow = 1, ncol = length(paramNames)))
      colnames(blanks) <- paramNames
      blanks$wave <- I(list(NA))  # <-- ensure list-column exists
      result[[chan]] <- blanks
      next
    }
    
    if (sr_hz == "auto") sr <- data$sr else sr <- sr_hz
    
    if (inherits(filterfrom_khz, "Arma")) {
      thisWave <- signal::filtfilt(filterfrom_khz, thisWave)@.Data
    } else if (filterfrom_khz > 0) {
      if (!is.null(filterto_khz)) {
        to_hz <- filterto_khz * 1000
        if (to_hz >= sr/2) to_hz <- NULL
      } else {
        to_hz <- NULL
      }
      thisWave <- seewave::bwfilter(thisWave, f = sr, n = 4,
                                    from = filterfrom_khz * 1000,
                                    to   = to_hz,
                                    output = "sample")
    }
    
    peakTime <- (which.max(abs(thisWave)) - 1)/sr
    
    fftSize <- round(sr * winLen_sec, 0)
    fftSize <- fftSize + (fftSize%%2)
    
    thisWave <- PAMpal:::clipAroundPeak(thisWave, fftSize)
    
    # ---- NEW: keep the exact waveform used for spectral calcs
    liftered_wave_used <- thisWave
    
    thisTk <- seewave:::TKEO(thisWave, f = sr, M = 1, plot = FALSE)
    tkEnergy <- thisTk[1:length(thisWave), 2]
    tkDb <- 10 * log10(tkEnergy - min(tkEnergy, na.rm = TRUE))
    tkDb <- tkDb - max(tkDb, na.rm = TRUE)
    tkDb[!is.finite(tkDb)] <- NA
    
    noiseLevel <- median(tkDb, na.rm = TRUE)
    if (is.na(noiseLevel)) noiseLevel <- 0
    
    thisDf <- list(noiseLevel = noiseLevel)
    
    noiseThresh <- quantile(thisTk[, 2], probs = 0.4, na.rm = TRUE) * 100
    dur <- subset(thisTk, thisTk[, 2] >= noiseThresh)
    
    if (length(dur) == 0) dur <- 0 else dur <- 1e+06 * (max(dur[, 1]) - min(dur[, 1]))
    
    thisDf$duration <- dur
    thisDf$peakTime <- peakTime
    
    thisSpec <- seewave:::spec(thisWave, f = sr, wl = fftSize,
                               norm = FALSE,
                               correction = "amplitude",
                               plot = FALSE)
    
    if (any(thisSpec[, 2] == 0)) thisSpec[thisSpec[, 2] == 0, 2] <- 1e-13
    
    if (any(is.nan(thisSpec[, 2]))) {
      blanks <- data.frame(matrix(NA, nrow = 1, ncol = length(paramNames)))
      colnames(blanks) <- paramNames
      blanks$wave <- I(list(NA))
      result[[chan]] <- blanks
      next
    }
    
    relDb <- 20 * log10(thisSpec[, 2])
    if (any(!is.finite(relDb))) relDb[!is.finite(relDb)] <- NA
    
    freq <- seq(from = 0,
                by   = thisSpec[2, 1] - thisSpec[1, 1],
                length.out = nrow(thisSpec))
    
    if (!is.null(calibration)) {
      if (is.function(calibration)) {
        calFun <- calibration
      } else if (is.character(calibration)) {
        calFun <- findCalibration(calibration)
      }
      relDb <- relDb + calFun(freq * 1000)
    }
    
    calibratedClick <- cbind(freq, relDb)
    
    peakData <- lapply(peakTrough(calibratedClick), unname)
    thisDf <- c(thisDf, peakData)
    
    dBPP <- 20 * log10(max(thisWave) - min(thisWave))
    if (!is.null(calibration)) dBPP <- dBPP + calFun(peakData$peak * 1000)
    thisDf$dBPP <- dBPP
    
    dbBW10 <- PAMpal:::Qfast(calibratedClick, f = sr, level = -10, plot = FALSE)
    names(dbBW10) <- c("Q_10dB", "fmin_10dB", "fmax_10dB", "BW_10dB")
    dbBW10$centerkHz_10dB <- dbBW10$fmax_10dB - (dbBW10$BW_10dB/2)
    
    dbBW3 <- PAMpal:::Qfast(calibratedClick, f = sr, level = -3, plot = FALSE)
    names(dbBW3) <- c("Q_3dB", "fmin_3dB", "fmax_3dB", "BW_3dB")
    dbBW3$centerkHz_3dB <- dbBW3$fmax_3dB - (dbBW3$BW_3dB/2)
    
    thisDf <- c(thisDf, dbBW10, dbBW3)
    
    # ---- NEW: add waveform as a list-column under base name 'wave'
    thisDf$wave <- list(liftered_wave_used)
    
    result[[chan]] <- thisDf
  }
  
  result <- dplyr::bind_rows(result)
  
  # Prefix all names with liftered_
  names(result) <- paste0("liftered_", names(result))
  
  # Waveform column exactly: liftered_wave
  result
}


# ---- Plotting helpers left unchanged below ----

clickSpectraMatrix <- function(waves, sr, nfft = 512) {
  if (!is.list(waves) || length(waves) == 0) {
    stop("'waves' must be a non-empty list of numeric vectors.")
  }
  if (!all(vapply(waves, is.numeric, logical(1)))) {
    stop("All elements of 'waves' must be numeric vectors.")
  }
  
  Ls <- vapply(waves, length, integer(1))
  if (length(unique(Ls)) != 1) {
    stop("All waveforms in 'waves' must have the same length.")
  }
  
  freqs <- seq(0, sr/2, length.out = nfft/2 + 1)
  
  click_spec <- function(x) {
    if (length(x) < nfft) {
      x <- c(x, rep(0, nfft - length(x)))
    } else if (length(x) > nfft) {
      x <- x[seq_len(nfft)]
    }
    X <- fft(x)
    P <- (Mod(X)^2) / (nfft * sr)
    P[1:(nfft/2 + 1)]
  }
  
  spec_mat <- sapply(waves, click_spec)
  spec_db  <- 10 * log10(spec_mat + .Machine$double.eps)
  
  list(freqs = freqs, spec_db = spec_db)
}

plotClickSpectrogram <- function(
    waves,
    sr,
    nfft     = 512,
    flim     = NULL,
    log_freq = FALSE,
    main     = "Click spectra by click number"
) {
  cs <- clickSpectraMatrix(waves, sr, nfft)
  freqs   <- cs$freqs
  spec_db <- cs$spec_db
  
  if (!is.null(flim)) {
    if (length(flim) != 2) stop("flim must be c(fmin, fmax)")
    keep <- which(freqs >= flim[1] & freqs <= flim[2])
    freqs   <- freqs[keep]
    spec_db <- spec_db[keep, , drop = FALSE]
  }
  
  n_clicks <- ncol(spec_db)
  
  if (!log_freq) {
    image(
      x = 1:n_clicks,
      y = freqs,
      z = t(spec_db),
      xlab = "Click number",
      ylab = "Frequency (Hz)",
      main = main
    )
  } else {
    plot(1,
         type = "n",
         xlim = c(1, n_clicks),
         ylim = range(freqs),
         log  = "y",
         xlab = "Click number",
         ylab = "Frequency (Hz)",
         main = main)
    
    z <- t(spec_db)
    z_norm <- (z - min(z, na.rm = TRUE)) /
      (max(z, na.rm = TRUE) - min(z, na.rm = TRUE))
    
    rasterImage(
      as.raster(z_norm),
      xleft   = 1,
      xright  = n_clicks,
      ybottom = min(freqs),
      ytop    = max(freqs)
    )
  }
  
  invisible(list(freqs = freqs, S_click_db = spec_db))
}

plotMeanClickSpectrum <- function(
    waves,
    sr,
    nfft     = 512,
    flim     = NULL,
    log_freq = FALSE,
    main     = "Mean spectrum across clicks"
) {
  cs <- clickSpectraMatrix(waves, sr, nfft)
  freqs   <- cs$freqs
  spec_db <- cs$spec_db
  
  if (!is.null(flim)) {
    if (length(flim) != 2) stop("flim must be c(fmin, fmax)")
    keep <- which(freqs >= flim[1] & freqs <= flim[2])
    freqs   <- freqs[keep]
    spec_db <- spec_db[keep, , drop = FALSE]
  }
  
  mean_spec_db <- rowMeans(spec_db, na.rm = TRUE)
  freqs_khz <- freqs / 1000
  
  if (!log_freq) {
    plot(
      freqs_khz, mean_spec_db,
      type = "l",
      xlab = "Frequency (kHz)",
      ylab = "Mean spectrum (dB)",
      main = main
    )
  } else {
    plot(
      freqs_khz, mean_spec_db,
      type = "l",
      log  = "x",
      xlab = "Frequency (kHz)",
      ylab = "Mean spectrum (dB)",
      main = paste0(main, " (log freq)")
    )
  }
  
  invisible(list(freqs_khz = freqs_khz, mean_spec_db = mean_spec_db))
}

