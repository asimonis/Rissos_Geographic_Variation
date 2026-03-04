library(MSnbase)
library(xcms)

pick_peaks_one_click <- function(freq_hz, spec_db,
                                 fmin_hz = 19000, fmax_hz = 50000,
                                 fwhm_hz = 1500,     # expected peak width (tune!)
                                 snthresh = 5,
                                 max_peaks = 10) {
  
  keep <- freq_hz >= fmin_hz & freq_hz <= fmax_hz
  f <- freq_hz[keep]
  y_db <- spec_db[keep]
  
  # Prefer linear intensity for peak finding
  y_lin <- 10^(y_db / 20)
  
  chr <- MSnbase::Chromatogram(rtime = f, intensity = y_lin)
  
  xchr <- xcms::findChromPeaks(
    chr,
    param = xcms::MatchedFilterParam(
      fwhm = fwhm_hz,
      snthresh = snthresh,
      max = max_peaks
    )
  )
  
  # Extract detected peaks (rt is your frequency axis here)
  pk <- xcms::chromPeaks(xchr)
  if (nrow(pk) == 0) return(pk)
  
  # Add frequency (Hz) columns for convenience
  pk$freq_hz   <- pk$rt
  pk$freqmin_hz <- pk$rtmin
  pk$freqmax_hz <- pk$rtmax
  pk
}





pick_peaks_all_clicks <- function(raw,
                                  fmin_hz = 19000, fmax_hz = 50000,
                                  fwhm_hz = 1500,
                                  snthresh = 5,
                                  max_peaks = 10) {
  
  freq_hz <- raw$freq
  out <- vector("list", ncol(raw$allSpec))
  
  for (i in seq_len(ncol(raw$allSpec))) {
    spec_db <- raw$allSpec[, i]
    
    pk <- pick_peaks_one_click(
      freq_hz = freq_hz,
      spec_db = spec_db,
      fmin_hz = fmin_hz, fmax_hz = fmax_hz,
      fwhm_hz = fwhm_hz,
      snthresh = snthresh,
      max_peaks = max_peaks
    )
    
    if (nrow(pk) > 0) {
      pk$UID <- raw$UID[i]
      pk$click_index <- i
    }
    out[[i]] <- pk
  }
  
  do.call(rbind, out)
}
