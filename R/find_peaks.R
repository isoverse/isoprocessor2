#' Peak Detection for TFS 253 Plus IRMS Data
#'
#' Three-stage pipeline reverse-engineered from Qtegra's PeakFinder.dll:
#'
#' 1. \code{detect_peaks()} — slope-based state machine finds peak boundaries
#' 2. \code{integrate_peaks()} — background subtraction, area, center time,
#'    square-pulse classification
#' 3. \code{filter_peaks()} — post-detection filtering on height, width, etc.
#'
#' Each stage can be called independently or chained together via
#' \code{find_peaks()}, which runs all three with a single call.

# ==============================================================================
# Stage 1: Peak detection (state machine)
# ==============================================================================

#' Detect peak boundaries using slope-based state machine
#'
#' Runs the PeakFinder.dll detection algorithm: computes rolling-window
#' linear-regression slopes, then walks through the data with a 3-state
#' machine (BACKGROUND -> PEAKSTART -> PEAKTOP) to locate peak start, top,
#' and end positions.
#'
#' @param time Numeric vector. Time in seconds.
#' @param intensity Numeric vector. Intensity in volts (V).
#' @param slope_n Integer. Rolling slope regression window size
#'   (\code{SlopeNumberOfDataPoints}). Default 5.
#' @param start_slope Numeric. Slope threshold for peak start detection
#'   (\code{StartSlope}). Peak start triggers when
#'   \code{slope1 + slope2 > 2 * start_slope}. Default 0.001.
#' @param top_slope Numeric. Slope threshold for peak top zero-crossing
#'   (\code{TopSlope}). Default 5e-05.
#' @param end_slope Numeric. Slope threshold for peak end detection
#'   (\code{EndSlope}). Default 0.001.
#' @param max_peak_width Numeric. Maximum allowed peak width in seconds
#'   (\code{MaxPeakWidth}). Peaks exceeding this are force-ended. Default 180.
#' @param peak_height_pct Numeric. Percentage drop from top required before
#'   end-slope is checked (\code{PeakHeightPercentage}). Default 20.
#' @param eval_start Numeric or NULL. Evaluation start time (s)
#'   (\code{EvalStartTime}). Default NULL (use all data).
#' @param eval_end Numeric or NULL. Evaluation end time (s)
#'   (\code{EvalEndTime}). Default NULL (use all data).
#' @param bg_history Integer. Number of preceding points for background
#'   estimation at peak start (\code{BackgroundHistoryPoints}). Default 25.
#'
#' @return A tibble with one row per detected peak:
#'   \describe{
#'     \item{peak}{Integer. Sequential peak number.}
#'     \item{start_idx}{Integer. Index into the time/intensity vectors.}
#'     \item{top_idx}{Integer. Index of peak maximum.}
#'     \item{end_idx}{Integer. Index of peak end.}
#'     \item{start_time}{Numeric. Peak start time (s).}
#'     \item{top_time}{Numeric. Peak top time (s).}
#'     \item{end_time}{Numeric. Peak end time (s).}
#'     \item{start_intensity}{Numeric. Raw intensity at peak start (V).}
#'     \item{top_intensity}{Numeric. Raw intensity at peak top (V).}
#'     \item{end_intensity}{Numeric. Raw intensity at peak end (V).}
#'     \item{bg_at_start}{Numeric. Rough background estimate at peak start
#'       (mean of preceding \code{bg_history} points). Used internally by
#'       the \code{"linear"} background method; superseded by the
#'       \code{"individual"} method's smoothed-minimum calculation.}
#'     \item{is_shoulder}{Logical. TRUE if peak was split from a shoulder.}
#'   }
#'
#' @export
detect_peaks <- function(
  time,
  intensity,
  slope_n = 5L,
  start_slope = 0.001,
  top_slope = 5e-05,
  end_slope = 0.001,
  max_peak_width = 180,
  peak_height_pct = 20,
  eval_start = NULL,
  eval_end = NULL,
  bg_history = 25L
) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required")
  }

  stopifnot(
    is.numeric(time),
    is.numeric(intensity),
    length(time) == length(intensity),
    length(time) >= 3L
  )

  n <- length(time)

  # Apply evaluation time window
  if (!is.null(eval_start) || !is.null(eval_end)) {
    keep <- rep(TRUE, n)
    if (!is.null(eval_start)) {
      keep <- keep & (time >= eval_start)
    }
    if (!is.null(eval_end)) {
      keep <- keep & (time <= eval_end)
    }
    time <- time[keep]
    intensity <- intensity[keep]
    n <- length(time)
    if (n < 3L) return(.empty_detected())
  }

  # Rolling slope with asymmetric window matching DLL layout
  slopes <- .rolling_slope(time, intensity, as.integer(slope_n))
  valid_slopes <- which(!is.na(slopes))
  if (length(valid_slopes) < 3L) {
    return(.empty_detected())
  }

  # --- State machine ---
  # States: 0 = BACKGROUND, 1 = PEAKSTART, 2 = PEAKTOP
  #
  # DLL's slope buffer holds 3 consecutive slopes. GetCurrentSlopes
  # returns buffer[0] (oldest) and buffer[2] (newest) — 2 cycles apart.
  # We replicate this by pairing valid_slopes[k-2] with valid_slopes[k].
  state <- 0L
  peaks <- list()
  peak_count <- 0L

  start_idx <- NA_integer_
  top_idx <- NA_integer_
  top_val_y <- -Inf
  bg_at_start <- 0
  is_shoulder <- FALSE

  # Helper to record a detected peak
  record_peak <- function() {
    peak_count <<- peak_count + 1L
    peaks[[peak_count]] <<- data.frame(
      start_idx = start_idx,
      top_idx = top_idx,
      end_idx = ci,
      start_time = time[start_idx],
      top_time = time[top_idx],
      end_time = time[ci],
      start_intensity = intensity[start_idx],
      top_intensity = intensity[top_idx],
      end_intensity = intensity[ci],
      bg_at_start = bg_at_start,
      is_shoulder = is_shoulder,
      stringsAsFactors = FALSE
    )
  }

  for (k in seq(3L, length(valid_slopes))) {
    # slope1 = buffer[0] (2 cycles ago), slope2 = buffer[2] (current)
    idx_s1 <- valid_slopes[k - 2L]
    idx_s2 <- valid_slopes[k]
    slope1 <- slopes[idx_s1]
    slope2 <- slopes[idx_s2]

    # Target data point is at the position of the newest slope
    # (.rolling_slope already accounts for the -1 offset via targetIndex)
    ci <- idx_s2
    cx <- time[ci]
    cy <- intensity[ci]

    if (state == 0L) {
      # BACKGROUND: look for rising slope
      sum_slopes <- slope1 + slope2
      if (sum_slopes > -1 && sum_slopes > 2 * start_slope) {
        state <- 1L
        start_idx <- ci
        is_shoulder <- FALSE
        bg_range <- max(1L, ci - bg_history):max(1L, ci - 1L)
        bg_at_start <- mean(intensity[bg_range])
      }
    } else if (state == 1L) {
      # PEAKSTART: look for peak top (slope zero-crossing)
      elapsed <- cx - time[start_idx]
      if (elapsed > max_peak_width) {
        state <- 0L
        next
      }
      # FoundPeakTop: slope transitions from > top_slope to < top_slope
      if (slope2 < top_slope && slope1 > top_slope) {
        bg_subtracted <- cy - bg_at_start
        if (bg_subtracted < 0) {
          # Below background — not a real peak start, abort
          state <- 0L
          next
        }
        state <- 2L
        top_idx <- ci
        top_val_y <- bg_subtracted
      }
    } else if (state == 2L) {
      # PEAKTOP (descending): look for peak end
      elapsed <- cx - time[start_idx]

      # Update top if current bg-subtracted intensity is higher
      bg_sub_now <- cy - bg_at_start
      if (bg_sub_now > top_val_y) {
        top_val_y <- bg_sub_now
        top_idx <- ci
      }

      # Force-end if max width exceeded
      if (elapsed > max_peak_width) {
        record_peak()
        state <- 0L
        next
      }

      # BelowDropIntyCriteria: has intensity dropped enough from top?
      top_raw <- intensity[top_idx]
      drop_pct <- round(cy * 100 / top_raw + 0.5)
      below_drop <- (drop_pct <= (100 - peak_height_pct))

      if (below_drop) {
        sum_slopes <- slope1 + slope2
        if (sum_slopes <= -2 * end_slope) {
          next # slope still steeply negative — continue descending
        }

        # Check for shoulder: new rising inflection on the tail
        if (slope1 < top_slope && slope2 > top_slope) {
          record_peak()
          # Start new peak from this point
          state <- 1L
          start_idx <- ci
          bg_at_start <- intensity[ci]
          is_shoulder <- TRUE
          top_val_y <- -Inf
        } else {
          record_peak()
          state <- 0L
        }
      }
    }
  }

  if (peak_count == 0L) {
    return(.empty_detected())
  }

  out <- do.call(rbind, peaks)
  out$peak <- seq_len(nrow(out))
  out <- out[, c(
    "peak",
    "start_idx",
    "top_idx",
    "end_idx",
    "start_time",
    "top_time",
    "end_time",
    "start_intensity",
    "top_intensity",
    "end_intensity",
    "bg_at_start",
    "is_shoulder"
  )]
  tibble::as_tibble(out)
}


# ==============================================================================
# Stage 2: Peak integration & calculations
# ==============================================================================

#' Integrate detected peaks and classify shape
#'
#' Given raw data and peak boundaries from \code{detect_peaks()}, computes
#' background-subtracted area (trapezoidal), height, intensity-weighted
#' center time, and the square-pulse flag.
#'
#' @param time Numeric vector. Time in seconds (same as passed to
#'   \code{detect_peaks()}).
#' @param intensity Numeric vector. Intensity in volts (same as passed to
#'   \code{detect_peaks()}).
#' @param peaks Tibble returned by \code{detect_peaks()}.
#' @param bg_method Character. Background subtraction method:
#'   \code{"individual"} (default, matches DLL BackgroundMethod.INDIVIDUAL)
#'   smooths the pre-peak baseline with a weighted 5-point window, takes
#'   the minimum as a constant background;
#'   \code{"linear"} interpolates between start and end baseline values;
#'   \code{"rolling_min"} uses the minimum of start and end intensity.
#' @param bg_history Integer. Number of pre-peak data points for background
#'   estimation (\code{BackgroundHistoryPoints}). Default 25.
#' @param bg_smooth_points Integer. Number of smoothing points for the
#'   \code{"individual"} method (\code{NumSmoothingDataPoints}, must be 5).
#'   Default 5.
#' @param square_pulse_detection Logical. Enable square-pulse classification
#'   (\code{IsSquarePulseDetectionEnabled}). Default TRUE.
#' @param square_pulse_factor Numeric. Rectangularity threshold
#'   (\code{SquarePulseDetectionFactor}). Peaks with
#'   \code{area / (time_to_top * top_height) > factor} are flagged as
#'   square pulses (reference gas injections). Default 0.55.
#'
#' @return The input tibble with additional columns:
#'   \describe{
#'     \item{width}{Numeric. Peak width (s), end_time - start_time.}
#'     \item{height}{Numeric. Max background-subtracted intensity (V).}
#'     \item{area}{Numeric. Integrated peak area (V*s), trapezoidal.}
#'     \item{center_time}{Numeric. Intensity-weighted center time (s).}
#'     \item{background}{Numeric. Background level subtracted (V). For
#'       \code{"individual"} this is a constant (smoothed minimum of
#'       pre-peak baseline); for \code{"linear"} it is the value at peak
#'       start (the end-value can be derived from the linear model).}
#'     \item{n_points}{Integer. Number of data points in the peak.}
#'     \item{is_square_pulse}{Logical. TRUE if peak is a square pulse.}
#'   }
#'
#' @export
integrate_peaks <- function(
  time,
  intensity,
  peaks,
  bg_method = c("individual", "linear", "rolling_min"),
  bg_history = 25L,
  bg_smooth_points = 5L,
  square_pulse_detection = TRUE,
  square_pulse_factor = 0.55
) {
  bg_method <- match.arg(bg_method)

  if (nrow(peaks) == 0L) {
    peaks$width <- numeric(0)
    peaks$height <- numeric(0)
    peaks$area <- numeric(0)
    peaks$center_time <- numeric(0)
    peaks$background <- numeric(0)
    peaks$n_points <- integer(0)
    peaks$is_square_pulse <- logical(0)
    return(peaks)
  }

  widths <- numeric(nrow(peaks))
  heights <- numeric(nrow(peaks))
  areas <- numeric(nrow(peaks))
  center_ts <- numeric(nrow(peaks))
  bg_vals <- numeric(nrow(peaks))
  n_pts <- integer(nrow(peaks))
  is_square <- logical(nrow(peaks))

  for (i in seq_len(nrow(peaks))) {
    si <- peaks$start_idx[i]
    ti <- peaks$top_idx[i]
    ei <- peaks$end_idx[i]

    t_start <- time[si]
    t_top <- time[ti]
    t_end <- time[ei]

    peak_idx <- si:ei
    peak_time <- time[peak_idx]
    peak_inty <- intensity[peak_idx]
    np <- length(peak_idx)

    # Background subtraction
    if (bg_method == "individual") {
      # DLL's IndividualBackground.GetBGD:
      # 1. Take rolling window of (bg_history + bg_smooth_points) pre-peak points
      # 2. Smooth with 5-point Standard coefficients [0.16, 0.22, 0.24, 0.22, 0.16]
      # 3. Take minimum of all smoothed values as constant background
      bg_val <- .individual_background(
        intensity,
        si,
        bg_history,
        bg_smooth_points
      )
      bg <- pmin(bg_val, peak_inty) # DLL caps bg at signal level
    } else if (bg_method == "linear") {
      bg_val <- peaks$bg_at_start[i]
      if (np > 1L) {
        bg_end_val <- intensity[ei]
        frac <- (peak_time - t_start) / max(t_end - t_start, 1e-15)
        bg <- bg_val + frac * (bg_end_val - bg_val)
      } else {
        bg <- bg_val
      }
    } else {
      # rolling_min
      bg_val <- min(peaks$bg_at_start[i], intensity[ei])
      bg <- rep(bg_val, np)
    }

    corrected <- peak_inty - bg
    corrected[corrected < 0] <- 0

    # Trapezoidal: Area = sum(Y) * (peakWidth / numPoints)
    peak_width <- t_end - t_start
    sum_y <- sum(corrected)
    area <- sum_y * (peak_width / np)

    # Intensity-weighted center time
    sum_xy <- sum(peak_time * corrected)
    ct <- if (sum_y > 1e-15) sum_xy / sum_y else t_top

    # Height: max background-subtracted intensity
    height <- max(corrected)

    # Square-pulse classification (SquarePeakDetection from DLL)
    # DLL IL: ratio = RawArea / ((EndVal.X - StartVal.X) * TopVal.Y)
    #   - denominator uses full peak width (end - start), NOT time-to-top
    #   - TopVal.Y is RAW intensity (not background-subtracted)
    # If ratio > SquarePulseDetectionFactor => square pulse
    sq <- FALSE
    if (square_pulse_detection) {
      raw_top_y <- intensity[ti] # raw, not bg-subtracted
      if (peak_width > 0 && raw_top_y > 0) {
        ratio <- area / (peak_width * raw_top_y)
        sq <- (ratio > square_pulse_factor)
      }
    }

    widths[i] <- peak_width
    heights[i] <- height
    areas[i] <- area
    center_ts[i] <- ct
    bg_vals[i] <- bg_val
    n_pts[i] <- np
    is_square[i] <- sq
  }

  peaks$width <- widths
  peaks$height <- heights
  peaks$area <- areas
  peaks$center_time <- center_ts
  peaks$background <- bg_vals
  peaks$n_points <- n_pts
  peaks$is_square_pulse <- is_square

  peaks
}


# ==============================================================================
# Stage 3: Post-detection filtering
# ==============================================================================

#' Filter detected peaks by height, width, and other criteria
#'
#' Applies post-detection filters matching PeakFinder.dll's PeakFilter.
#' Removes peaks that fail the criteria and renumbers the remaining peaks.
#'
#' @param peaks Tibble returned by \code{integrate_peaks()}.
#' @param min_height Numeric. Minimum background-subtracted peak height (V)
#'   (\code{MinHeight}). Default 0.05.
#' @param min_height_post Numeric. Post-filter minimum height (V)
#'   (\code{MinHeightPostFilter}). Applied after integration. Default 0
#'   (disabled).
#' @param max_peak_width_post Numeric. Post-filter maximum peak width (s)
#'   (\code{MaxPeakWidthPostFilter}). Default 0 (disabled).
#' @param min_area Numeric. Minimum integrated area (V*s). Default 0
#'   (disabled).
#'
#' @return Filtered tibble with peaks renumbered.
#'
#' @export
filter_peaks <- function(
  peaks,
  min_height = 0.05,
  min_height_post = 0,
  max_peak_width_post = 0,
  min_area = 0
) {
  if (nrow(peaks) == 0L) {
    return(peaks)
  }

  keep <- rep(TRUE, nrow(peaks))

  # MinHeight — core filter (also applied in DLL state machine, but
  # with bg_at_start only; here we use the fully integrated height)
  if (min_height > 0) {
    keep <- keep & (peaks$height >= min_height)
  }

  # MinHeightPostFilter
  if (min_height_post > 0) {
    keep <- keep & (peaks$height >= min_height_post)
  }

  # MaxPeakWidthPostFilter
  if (max_peak_width_post > 0) {
    keep <- keep & (peaks$width <= max_peak_width_post)
  }

  # MinArea (convenience, not in DLL but commonly useful)
  if (min_area > 0) {
    keep <- keep & (peaks$area >= min_area)
  }

  out <- peaks[keep, , drop = FALSE]
  out$peak <- seq_len(nrow(out))
  out
}


# ==============================================================================
# Convenience wrapper: all three stages in one call
# ==============================================================================

#' Find, integrate, and filter peaks in one call
#'
#' Convenience wrapper that chains \code{detect_peaks()},
#' \code{integrate_peaks()}, and \code{filter_peaks()}.
#'
#' @inheritParams detect_peaks
#' @inheritParams integrate_peaks
#' @inheritParams filter_peaks
#'
#' @return A tibble with one row per detected peak. See
#'   \code{\link{integrate_peaks}} for column descriptions.
#'
#' @examples
#' \dontrun{
#'   library(dplyr)
#'   library(ggplot2)
#'
#'   dat <- read_measure_data("path/to/MeasureData.bin", long = TRUE)
#'   m44 <- dat |> filter(mass == 44)
#'   peaks <- find_peaks(m44$time_s, m44$intensity_V)
#'
#'   ggplot(m44, aes(time_s, intensity_V)) +
#'     geom_line() +
#'     geom_point(data = peaks, aes(x = top_time, y = top_intensity),
#'                colour = "red", size = 2)
#' }
#'
#' @export
find_peaks <- function(
  time,
  intensity,
  # detect_peaks params
  slope_n = 5L,
  start_slope = 0.001,
  top_slope = 5e-05,
  end_slope = 0.001,
  max_peak_width = 180,
  peak_height_pct = 20,
  eval_start = NULL,
  eval_end = NULL,
  bg_history = 25L,
  # integrate_peaks params
  bg_method = c("individual", "linear", "rolling_min"),
  bg_smooth_points = 5L,
  square_pulse_detection = TRUE,
  square_pulse_factor = 0.55,
  # filter_peaks params
  min_height = 0.05,
  min_height_post = 0,
  max_peak_width_post = 0,
  min_area = 0
) {
  bg_method <- match.arg(bg_method)

  peaks <- detect_peaks(
    time,
    intensity,
    slope_n = slope_n,
    start_slope = start_slope,
    top_slope = top_slope,
    end_slope = end_slope,
    max_peak_width = max_peak_width,
    peak_height_pct = peak_height_pct,
    eval_start = eval_start,
    eval_end = eval_end,
    bg_history = bg_history
  )

  peaks <- integrate_peaks(
    time,
    intensity,
    peaks,
    bg_method = bg_method,
    bg_history = bg_history,
    bg_smooth_points = bg_smooth_points,
    square_pulse_detection = square_pulse_detection,
    square_pulse_factor = square_pulse_factor
  )

  peaks <- filter_peaks(
    peaks,
    min_height = min_height,
    min_height_post = min_height_post,
    max_peak_width_post = max_peak_width_post,
    min_area = min_area
  )

  peaks
}


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Rolling linear-regression slope (matching PeakFinder.dll layout)
#'
#' The DLL's FindSlope uses an asymmetric regression window:
#'   additionalSpace = slope_n %/% 2
#'   bufferSize      = slope_n + additionalSpace
#'   targetIndex     = (bufferSize - additionalSpace) %/% 2 - 1 + additionalSpace
#' Regression runs over the last slope_n points of the buffer (indices
#' additionalSpace to bufferSize-1). The target (the data point the slope
#' is attributed to) sits at targetIndex, which is 1 position before the
#' regression center.
#'
#' @param x,y Numeric vectors (time, intensity).
#' @param slope_n Integer. SlopeNumberOfDataPoints.
#' @return Numeric vector of slopes, NA where the window doesn't fit.
#' @keywords internal
.rolling_slope <- function(x, y, slope_n) {
  n <- length(x)
  additional <- slope_n %/% 2L
  buffer_size <- slope_n + additional
  target_idx <- (buffer_size - additional) %/% 2L - 1L + additional

  # Offsets relative to the target point
  behind <- target_idx - additional # points before target in regression
  ahead <- (buffer_size - 1L) - target_idx # points after target in regression

  slopes <- rep(NA_real_, n)
  for (i in seq(behind + 1L, n - ahead)) {
    idx <- seq(i - behind, i + ahead)
    xi <- x[idx]
    yi <- y[idx]
    mx <- mean(xi)
    my <- mean(yi)
    denom <- sum((xi - mx)^2)
    if (denom > 0) {
      slopes[i] <- sum((xi - mx) * (yi - my)) / denom
    }
  }
  slopes
}

#' Compute INDIVIDUAL background level (matching DLL)
#'
#' IndividualBackground.GetBGD: takes the pre-peak baseline data, applies
#' 5-point Standard smoothing (coefficients 0.16, 0.22, 0.24, 0.22, 0.16),
#' and returns the minimum of all smoothed values as a constant background.
#'
#' @param intensity Numeric vector. Full intensity trace.
#' @param peak_start_idx Integer. Index of peak start.
#' @param bg_history Integer. BackgroundHistoryPoints (default 25).
#' @param bg_smooth_points Integer. NumSmoothingDataPoints (default 5).
#' @return Numeric scalar. Constant background level.
#' @keywords internal
.individual_background <- function(
  intensity,
  peak_start_idx,
  bg_history = 25L,
  bg_smooth_points = 5L
) {
  # Standard smoothing coefficients from StandardCoeff.Init in DLL
  smooth_coeffs <- c(0.16, 0.22, 0.24, 0.22, 0.16)

  # SK: attempt:
  window_end <- peak_start_idx + (length(smooth_coeffs) %/% 2)
  window_start <- peak_start_idx -
    (bg_history - 1L) -
    (length(smooth_coeffs) %/% 2)
  bg_window <- intensity[window_start:window_end]

  bg_smoothed <- bg_window |> stats::filter(filter = smooth_coeffs, sides = 2)
  return(min(bg_smoothed, na.rm = TRUE))

  # # Rolling window size = bg_history + bg_smooth_points (DLL: 25 + 5 = 30)
  # raw_window_size <- bg_history + bg_smooth_points

  # # Get pre-peak background data points
  # window_end <- max(1L, peak_start_idx + (length(smooth_coeffs) - 1) / 2)
  # window_start <- max(1L, peak_start_idx - raw_window_size)
  # bg_window <- intensity[window_start:window_end]
  # n_bg <- length(bg_window)

  # if (n_bg < bg_smooth_points) {
  #   # Not enough data — DLL fills window with copies of first point
  #   return(min(bg_window))
  # }

  # # If window is shorter than expected, pad with first value
  # # (DLL fills entire window with first point on initialization)
  # if (n_bg < raw_window_size) {
  #   bg_window <- c(rep(bg_window[1L], raw_window_size - n_bg), bg_window)
  #   n_bg <- raw_window_size
  # }

  # # Apply 5-point weighted smoothing and find minimum
  # # DLL smooths positions halfSmooth .. (n_bg - halfSmooth - 1)
  # half_smooth <- bg_smooth_points %/% 2L # = 2
  # n_smoothed <- n_bg - bg_smooth_points + 1L
  # min_val <- Inf

  # for (j in seq_len(n_smoothed)) {
  #   chunk <- bg_window[j:(j + bg_smooth_points - 1L)]
  #   smoothed <- sum(chunk * smooth_coeffs)
  #   if (smoothed < min_val) min_val <- smoothed
  # }

  # min_val
}


#' Parabolic (quadratic) fit to find sub-sample peak center
#'
#' Matches PeakFinder.dll's TimeShift.ParabolicFit: takes a window of data
#' points centered on the peak top, builds a Vandermonde-style normal equation
#' system, solves via Cramer's rule, and returns the parabola's vertex (the
#' sub-sample-resolution position of the peak maximum).
#'
#' @param time Numeric vector. Full time axis.
#' @param intensity Numeric vector. Full intensity trace.
#' @param top_idx Integer. Index of the peak top in the data vectors.
#' @param buffer_size Integer. Number of points for the fit (default 9;
#'   DLL uses 25 when integration_time < 0.2s).
#' @return Numeric scalar. The fitted peak center time.  Returns
#'   \code{time[top_idx]} if the fit fails.
#' @keywords internal
.parabolic_center <- function(time, intensity, top_idx, buffer_size = 9L) {
  n <- length(time)
  order <- 3L # quadratic: 3 coefficients (a, b, c)
  half <- buffer_size %/% 2L

  # Indices centered on peak top
  i_start <- top_idx - half
  i_end <- top_idx + (buffer_size - 1L) - half
  if (i_start < 1L || i_end > n) {
    return(time[top_idx])
  }

  idx <- i_start:i_end
  x <- time[idx]
  y <- intensity[idx]

  # Build smoothing sums: S[k] = sum(x^k), k = 0..(2*(order-1))
  # and right-hand side: R[k] = sum(y * x^k), k = 0..(order-1)
  smooth_size <- buffer_size %/% 2L + 1L
  S <- numeric(2L * (order - 1L) + 1L) # S[1]..S[5]
  R <- numeric(order)

  for (j in seq_along(idx)) {
    xi <- x[j]
    yi <- y[j]
    xpow <- 1.0
    for (k in seq_len(length(S))) {
      S[k] <- S[k] + xpow
      xpow <- xpow * xi
    }
    ypow <- yi
    for (k in seq_len(order)) {
      R[k] <- R[k] + ypow
      ypow <- ypow * xi
    }
  }

  # Build normal equation matrix (order x order)
  # M[i,j] = S[i+j-1]  (1-indexed), i.e. M[row,col] = sum(x^(row+col-2))
  M <- matrix(0, nrow = order, ncol = order)
  for (i in seq_len(order)) {
    for (j in seq_len(order)) {
      M[i, j] <- S[i + j - 1L]
    }
  }

  # Solve using Cramer's rule (matching DLL's Determinant method)
  det_M <- det(M)
  if (abs(det_M) < .Machine$double.eps) {
    return(time[top_idx])
  }

  coeffs <- numeric(order)
  for (k in seq_len(order)) {
    M_k <- M
    M_k[, k] <- R
    coeffs[k] <- det(M_k) / det_M
  }

  # Vertex of parabola ax^2 + bx + c:  center = -b / (2a)
  # coeffs[1] = c, coeffs[2] = b, coeffs[3] = a
  a <- coeffs[order] # x^2 coefficient
  b <- coeffs[order - 1L] # x^1 coefficient
  if (abs(a) < .Machine$double.eps) {
    return(time[top_idx])
  }

  center <- -b / (2.0 * a)
  center
}


#' Integrate peaks across multiple traces with time-shift correction
#'
#' Implements the full multi-trace integration pipeline from PeakFinder.dll:
#' \enumerate{
#'   \item Peaks are detected on the detection trace (e.g., mass 44).
#'   \item Each trace is integrated using the same peak boundaries.
#'   \item A parabolic fit finds the sub-sample peak center on each trace.
#'   \item Non-detection traces get a time-shift area correction based on
#'         the offset between their center and the detection trace's center.
#' }
#'
#' @param time Numeric vector. Shared time axis (seconds).
#' @param traces Named list of numeric vectors, one per trace. Names should
#'   be descriptive (e.g., \code{list("44" = v44, "45" = v45, "46" = v46)}).
#' @param detection_trace Character or integer. Name or index of the detection
#'   trace in \code{traces}. Default \code{1L} (first trace).
#' @param peaks Tibble from \code{detect_peaks()} on the detection trace.
#' @param bg_method Background method passed to \code{integrate_peaks()}.
#' @param bg_history Integer. BackgroundHistoryPoints. Default 25.
#' @param bg_smooth_points Integer. NumSmoothingDataPoints. Default 5.
#' @param square_pulse_detection Logical. Enable square-pulse classification.
#' @param square_pulse_factor Numeric. Threshold for square-pulse. Default 0.55.
#' @param perform_time_shift Logical. Apply time-shift correction to
#'   non-detection traces. Default TRUE.
#' @param max_time_shift Numeric. Maximum allowed time shift in seconds.
#'   Default 0.5.
#' @param parabolic_buffer Integer. Points for parabolic fit. Default 9.
#' @param gain_ratios Named numeric vector of amplifier gain ratios relative
#'   to the detection trace (e.g., \code{c("45" = 100, "46" = 333)}).
#'   Raw voltage areas for non-detection traces are divided by these values
#'   to produce areas comparable to the detection trace's scale. This
#'   compensates for different Faraday cup amplifier resistors. If \code{NULL}
#'   (default), areas are returned as raw voltage integrals.
#'
#' @return A tibble with one row per peak. Contains the detection trace columns
#'   (from \code{integrate_peaks()}) plus, for each additional trace, columns
#'   named \code{area_<trace>}, \code{background_<trace>}, and
#'   \code{time_shift_<trace>}.
#'
#' @export
integrate_peaks_multi <- function(
  time,
  traces,
  detection_trace = 1L,
  peaks = NULL,
  bg_method = c("individual", "linear", "rolling_min"),
  bg_history = 25L,
  bg_smooth_points = 5L,
  square_pulse_detection = TRUE,
  square_pulse_factor = 0.55,
  perform_time_shift = TRUE,
  max_time_shift = 0.5,
  parabolic_buffer = 9L,
  gain_ratios = NULL
) {
  bg_method <- match.arg(bg_method)

  # Resolve detection trace
  trace_names <- names(traces)
  if (is.null(trace_names)) {
    trace_names <- as.character(seq_along(traces))
    names(traces) <- trace_names
  }

  if (is.numeric(detection_trace)) {
    det_name <- trace_names[detection_trace]
  } else {
    det_name <- as.character(detection_trace)
  }
  det_intensity <- traces[[det_name]]

  # If no peaks supplied, detect on the detection trace
  if (is.null(peaks)) {
    peaks <- detect_peaks(time, det_intensity)
  }

  # Integrate detection trace (this gives the base result)
  result <- integrate_peaks(
    time,
    det_intensity,
    peaks,
    bg_method = bg_method,
    bg_history = bg_history,
    bg_smooth_points = bg_smooth_points,
    square_pulse_detection = square_pulse_detection,
    square_pulse_factor = square_pulse_factor
  )

  if (nrow(result) == 0L || length(traces) < 2L) {
    return(result)
  }

  # Compute parabolic center for detection trace (per peak)
  det_centers <- numeric(nrow(result))
  for (i in seq_len(nrow(result))) {
    det_centers[i] <- .parabolic_center(
      time,
      det_intensity,
      result$top_idx[i],
      parabolic_buffer
    )
  }

  # Compute dt (average seconds per point) for index-shift conversion
  dt <- if (length(time) > 1L) {
    (time[length(time)] - time[1L]) / (length(time) - 1L)
  } else {
    1.0
  }

  # Process each non-detection trace
  other_names <- setdiff(trace_names, det_name)

  for (tname in other_names) {
    tr_intensity <- traces[[tname]]

    # Integrate using the same peak boundaries (from detection trace)
    tr_integ <- integrate_peaks(
      time,
      tr_intensity,
      peaks,
      bg_method = bg_method,
      bg_history = bg_history,
      bg_smooth_points = bg_smooth_points,
      square_pulse_detection = FALSE,
      square_pulse_factor = square_pulse_factor
    )

    areas <- tr_integ$area
    backgrounds <- tr_integ$background
    time_shifts <- numeric(nrow(result))

    if (perform_time_shift) {
      for (i in seq_len(nrow(result))) {
        # Skip square pulses — DLL skips time shift for square peaks
        if (result$is_square_pulse[i]) {
          next
        }
        if (areas[i] <= 0) {
          next
        }

        # Parabolic center for this trace
        tr_center <- .parabolic_center(
          time,
          tr_intensity,
          result$top_idx[i],
          parabolic_buffer
        )

        # Time shift = difference between this trace's center and detection
        shift <- abs(tr_center) - abs(det_centers[i])

        # Clamp to max_time_shift
        shift_sign <- sign(shift)
        if (abs(shift) > max_time_shift) {
          shift <- shift_sign * max_time_shift
        }

        time_shifts[i] <- shift

        # Skip if shift is effectively zero
        if (abs(shift) < .Machine$double.eps) {
          next
        }

        # Convert shift to integer index offset
        idx_offset <- as.integer(shift / dt)

        # Standard path (non-extended): linear area correction
        # fractional_shift = shift - idx_offset * dt
        # area_correction  = (endY - startY) * fractional_shift
        si <- result$start_idx[i]
        ei <- result$end_idx[i]
        frac_shift <- shift - idx_offset * dt

        start_y <- tr_intensity[si]
        end_y <- tr_intensity[ei]
        height_diff <- end_y - start_y

        areas[i] <- areas[i] + height_diff * frac_shift
      }
    }

    # Apply gain ratio (amplifier resistor correction)
    if (!is.null(gain_ratios) && tname %in% names(gain_ratios)) {
      gr <- gain_ratios[[tname]]
      areas <- areas / gr
      backgrounds <- backgrounds / gr
    }

    # Add columns to result
    result[[paste0("area_", tname)]] <- areas
    result[[paste0("background_", tname)]] <- backgrounds
    result[[paste0("time_shift_", tname)]] <- time_shifts
  }

  result
}


#' Empty tibble for detect_peaks
#' @keywords internal
.empty_detected <- function() {
  tibble::tibble(
    peak = integer(0),
    start_idx = integer(0),
    top_idx = integer(0),
    end_idx = integer(0),
    start_time = numeric(0),
    top_time = numeric(0),
    end_time = numeric(0),
    start_intensity = numeric(0),
    top_intensity = numeric(0),
    end_intensity = numeric(0),
    bg_at_start = numeric(0),
    is_shoulder = logical(0)
  )
}


# ==============================================================================
# Convenience plotting
# ==============================================================================

#' Plot detected peaks on a chromatogram
#'
#' @param time Numeric vector. Time in seconds.
#' @param intensity Numeric vector. Intensity in volts.
#' @param peaks Tibble from \code{integrate_peaks()} or \code{find_peaks()}.
#' @param show_background Logical. If TRUE, draw background lines. Default TRUE.
#'
#' @return A ggplot2 object.
#' @export
plot_peaks <- function(time, intensity, peaks, show_background = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  df <- data.frame(time = time, intensity = intensity)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = intensity)) +
    ggplot2::geom_line(colour = "grey30", linewidth = 0.4) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Time (s)",
      y = "Intensity (V)",
      title = paste(nrow(peaks), "peaks detected")
    )

  if (nrow(peaks) > 0L) {
    p <- p +
      ggplot2::geom_rect(
        data = peaks,
        mapping = ggplot2::aes(
          xmin = start_time,
          xmax = end_time,
          ymin = -Inf,
          ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "steelblue",
        alpha = 0.08
      )

    p <- p +
      ggplot2::geom_point(
        data = peaks,
        mapping = ggplot2::aes(x = top_time, y = top_intensity),
        inherit.aes = FALSE,
        colour = "red",
        size = 1.5
      )

    p <- p +
      ggplot2::geom_text(
        data = peaks,
        mapping = ggplot2::aes(x = top_time, y = top_intensity, label = peak),
        inherit.aes = FALSE,
        vjust = -0.8,
        size = 3,
        colour = "red"
      )

    if (show_background && "background" %in% names(peaks)) {
      p <- p +
        ggplot2::geom_segment(
          data = peaks,
          mapping = ggplot2::aes(
            x = start_time,
            xend = end_time,
            y = background,
            yend = background
          ),
          inherit.aes = FALSE,
          colour = "forestgreen",
          linewidth = 0.5,
          linetype = "dashed"
        )
    }
  }

  p
}


# ==============================================================================
# Self-test
# ==============================================================================

if (identical(environment(), globalenv())) {
  cat("=== find_peaks.R self-test ===\n\n")

  # --- Synthetic IRMS data: 5 Gaussian peaks + 2 reference gas pulses ---
  # Mimics a typical GasBench CO2 measurement sequence:
  #   - 2 reference gas square pulses (instant rise/fall, 20s wide)
  #   - 5 Gaussian sample peaks from GC separation
  set.seed(42)
  dt <- 0.25 # 4 Hz (typical 253 Plus)
  t_max <- 800
  time <- seq(0, t_max, by = dt)
  n <- length(time)

  baseline <- 0.5 # 0.5 V baseline
  noise_sd <- 0.002
  intensity <- baseline + rnorm(n, sd = noise_sd)

  # 5 Gaussian peaks
  gauss_centers <- c(200, 300, 420, 540, 660)
  gauss_heights <- c(2.0, 1.5, 3.0, 1.0, 2.5)
  gauss_sigmas <- c(4, 5, 6, 3, 5)
  for (j in seq_along(gauss_centers)) {
    intensity <- intensity +
      gauss_heights[j] *
        exp(-0.5 * ((time - gauss_centers[j]) / gauss_sigmas[j])^2)
  }

  # 2 reference gas square pulses: instant rise, flat top, instant fall
  # These have time_to_top ≈ 0, so ratio = area / (tiny × height) → very large
  ref_centers <- c(50, 120)
  ref_height <- 2.0
  ref_half_w <- 10 # 20s wide total
  for (rc in ref_centers) {
    mask <- time >= (rc - ref_half_w) & time <= (rc + ref_half_w)
    intensity[mask] <- intensity[mask] + ref_height
  }

  cat(sprintf("Synthetic data: %d points, %.0fs\n", n, t_max))
  cat(sprintf(
    "  %d Gaussian sample peaks + %d square reference pulses\n",
    length(gauss_centers),
    length(ref_centers)
  ))

  # --- Stage 1: detect ---
  cat("\n--- Stage 1: detect_peaks() ---\n")
  det <- detect_peaks(
    time,
    intensity,
    slope_n = 7,
    start_slope = 0.005,
    top_slope = 0.001,
    end_slope = 0.005,
    max_peak_width = 60
  )
  cat(sprintf("  Detected %d raw peaks\n", nrow(det)))
  for (i in seq_len(nrow(det))) {
    p <- det[i, ]
    cat(sprintf(
      "  Peak %d: start=%.1f  top=%.1f  end=%.1fs  width=%.1fs\n",
      p$peak,
      p$start_time,
      p$top_time,
      p$end_time,
      p$end_time - p$start_time
    ))
  }

  # --- Stage 2: integrate ---
  cat("\n--- Stage 2: integrate_peaks() ---\n")
  integ <- integrate_peaks(
    time,
    intensity,
    det,
    square_pulse_detection = TRUE,
    square_pulse_factor = 0.55
  )
  for (i in seq_len(nrow(integ))) {
    p <- integ[i, ]
    tt <- p$top_time - p$start_time # time to top
    cat(sprintf(
      "  Peak %d: t=%.1f-%.1f-%.1fs  h=%.3fV  area=%.2f  t_to_top=%.1fs  sq=%s\n",
      p$peak,
      p$start_time,
      p$top_time,
      p$end_time,
      p$height,
      p$area,
      tt,
      p$is_square_pulse
    ))
  }

  # --- Stage 3: filter ---
  cat("\n--- Stage 3: filter_peaks() ---\n")
  filt <- filter_peaks(integ, min_height = 0.1)
  cat(sprintf("  After min_height=0.1V: %d peaks remain\n", nrow(filt)))

  # --- Verification ---
  cat("\n--- Verification ---\n")
  all_true <- c(ref_centers, gauss_centers)
  matched <- 0L
  for (i in seq_len(nrow(filt))) {
    diffs <- abs(filt$top_time[i] - all_true)
    best <- which.min(diffs)
    ok <- diffs[best] < 15
    matched <- matched + ok
    cat(sprintf(
      "  Peak %d: top=%.1fs  h=%.3fV  w=%.1fs  nearest_true=%.0fs (d=%.1f) %s\n",
      i,
      filt$top_time[i],
      filt$height[i],
      filt$width[i],
      all_true[best],
      diffs[best],
      if (ok) "OK" else "MISS"
    ))
  }
  cat(sprintf("\nMatched: %d/%d\n", matched, length(all_true)))

  if (
    matched == length(all_true) &&
      nrow(filt) == length(all_true)
  ) {
    cat("*** ALL PEAKS CORRECTLY DETECTED ***\n")
  }

  # --- Convenience wrapper test ---
  cat("\n--- find_peaks() wrapper ---\n")
  all_peaks <- find_peaks(
    time,
    intensity,
    slope_n = 7,
    start_slope = 0.005,
    top_slope = 0.001,
    end_slope = 0.005,
    min_height = 0.1,
    max_peak_width = 60
  )
  cat(sprintf(
    "  find_peaks() returned %d peaks (%d square, %d Gaussian)\n",
    nrow(all_peaks),
    sum(all_peaks$is_square_pulse),
    sum(!all_peaks$is_square_pulse)
  ))

  # --- Real IRMS data comparison with Qtegra ---
  real_csv <- "/Users/seko0922/Dropbox/Tools/software/R/isoreader2/tmp/qtegra/R_from_claude/raw_data_example.csv"
  if (file.exists(real_csv)) {
    cat("\n--- Real IRMS data: comparison with Qtegra ---\n")
    real <- read.csv(real_csv)
    real_t <- real$time.sec
    real_v <- real$voltage

    # Qtegra reference values (start_slope=0.02, end_slope=0.04)
    qteg <- data.frame(
      peak = 1:6,
      ret_time = c(33.4, 83.3, 133.5, 219.4, 284.5, 349.6),
      start = c(13.2, 63.2, 113.1, 216.5, 281.6, 346.7),
      end = c(35.7, 85.8, 135.8, 226.3, 291.2, 356.3),
      area = c(184.205, 184.295, 186.333, 20.200, 19.016, 17.935)
    )

    # Run with matching parameters
    real_peaks <- find_peaks(
      real_t,
      real_v,
      start_slope = 0.02,
      end_slope = 0.04,
      min_height = 0
    )

    cat(sprintf(
      "  Detected %d peaks (Qtegra: %d)\n",
      nrow(real_peaks),
      nrow(qteg)
    ))

    # Compare each background method
    for (method in c("individual", "linear")) {
      real_int <- integrate_peaks(
        real_t,
        real_v,
        detect_peaks(real_t, real_v, start_slope = 0.02, end_slope = 0.04),
        bg_method = method
      )
      cat(sprintf("\n  Background method: %s\n", method))
      cat(sprintf(
        "  %4s  %8s %8s  %8s %8s  %10s %10s  %7s\n",
        "Peak",
        "RetOurs",
        "RetQteg",
        "StartO",
        "StartQ",
        "AreaOurs",
        "AreaQteg",
        "Diff%"
      ))

      for (j in seq_len(min(nrow(real_int), nrow(qteg)))) {
        diff_pct <- (real_int$area[j] - qteg$area[j]) / qteg$area[j] * 100
        cat(sprintf(
          "  %4d  %8.1f %8.1f  %8.1f %8.1f  %10.3f %10.3f  %+6.2f%%\n",
          j,
          real_int$top_time[j],
          qteg$ret_time[j],
          real_int$start_time[j],
          qteg$start[j],
          real_int$area[j],
          qteg$area[j],
          diff_pct
        ))
      }
    }
  }

  # --- Multi-trace test with real IRMS data ---
  multi_csv <- "/Users/seko0922/Dropbox/Tools/software/R/isoreader2/tmp/qtegra/R_from_claude/raw_data_example_all_traces.csv"
  if (file.exists(multi_csv)) {
    cat("\n--- Multi-trace integration: comparison with Qtegra ---\n")
    multi <- read.csv(multi_csv)
    n_per_trace <- nrow(multi) / 3L

    t_all <- multi$time.sec[1:n_per_trace]
    v44 <- multi$voltage[1:n_per_trace]
    v45 <- multi$voltage[(n_per_trace + 1):(2 * n_per_trace)]
    v46 <- multi$voltage[(2 * n_per_trace + 1):(3 * n_per_trace)]

    traces <- list("44" = v44, "45" = v45, "46" = v46)

    # TFS 253 Plus amplifier gain ratios (relative to mass 44):
    #   mass 44: 3e8 Ω  (reference, gain = 1)
    #   mass 45: 3e10 Ω (gain = 100)
    #   mass 46: 1e11 Ω (gain ≈ 333)
    gains <- c("45" = 100, "46" = 333)

    # Detect on trace 44
    det44 <- detect_peaks(t_all, v44, start_slope = 0.02, end_slope = 0.04)

    # Multi-trace integration with gain correction
    multi_result <- integrate_peaks_multi(
      t_all,
      traces,
      detection_trace = "44",
      peaks = det44,
      bg_method = "individual",
      gain_ratios = gains
    )

    # Qtegra reference values (from screenshot)
    qteg_44 <- c(184.205, 184.295, 186.333, 20.200, 19.016, 17.935)
    qteg_45 <- c(2.132, 2.133, 2.156, 0.231, 0.217, 0.205)
    qteg_46 <- c(0.789, 0.789, 0.800, 0.086, 0.081, 0.077)

    cat("  Gain ratios: 45/44 = 100, 46/44 = 333\n\n")
    cat(sprintf(
      "  %4s  %7s  %9s %9s %8s  %9s %9s %8s  %9s %9s %8s  %8s %8s\n",
      "Peak",
      "square",
      "44_ours",
      "44_qteg",
      "diff%",
      "45_ours",
      "45_qteg",
      "diff%",
      "46_ours",
      "46_qteg",
      "diff%",
      "ts_45",
      "ts_46"
    ))

    for (j in seq_len(min(nrow(multi_result), 6L))) {
      d44 <- (multi_result$area[j] - qteg_44[j]) / qteg_44[j] * 100
      d45 <- (multi_result$area_45[j] - qteg_45[j]) / qteg_45[j] * 100
      d46 <- (multi_result$area_46[j] - qteg_46[j]) / qteg_46[j] * 100
      sq <- if (multi_result$is_square_pulse[j]) "YES" else "no"
      cat(sprintf(
        "  %4d  %7s  %9.3f %9.3f %+7.2f%%  %9.3f %9.3f %+7.2f%%  %9.3f %9.3f %+7.2f%%  %+8.4f %+8.4f\n",
        j,
        sq,
        multi_result$area[j],
        qteg_44[j],
        d44,
        multi_result$area_45[j],
        qteg_45[j],
        d45,
        multi_result$area_46[j],
        qteg_46[j],
        d46,
        multi_result$time_shift_45[j],
        multi_result$time_shift_46[j]
      ))
    }
  }

  cat("\n=== Self-test complete ===\n")
}
