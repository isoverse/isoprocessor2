# background detector class ========

#' Background detectors
#'
#' Background detectors describe how a peak's background (the baseline that is
#' subtracted from the signal) is determined. They are supplied to a peak
#' detector (e.g. [ip_slope_based_peak_detector()]) and run from within
#' [ip_detect_peaks()] after peak detection. Use the generic `ip_bgrd_detector()`
#' to build a custom one, or the built-in `ip_no_bgrd_detector()`,
#' `ip_individual_bgrd_detector()` and `ip_isodat_default_background()`.
#'
#' `ip_bgrd_detector()` creates an `ip_bgrd_detector` object: unlike an
#' [ip_peak_detector] (a tibble), it is a simple single object (a list) with
#' three elements:
#'
#' - `type` (character): the kind of background detector (shown by [print()]).
#' - `detect` (function): a function `function(traces, peak_traces)` that returns
#'   `peak_traces` with an added background column `bgrd.<unit>` (in the same unit
#'   as the `intensity.<unit>` column, e.g. `bgrd.mV`).
#' - `details` (character): an optional plain-text summary of the detector's
#'   parameters (shown by [print()]).
#'
#' @param type the kind of background detector (a single non-empty string).
#' @param detect a function `function(traces, peak_traces)` returning
#'   `peak_traces` with an added `bgrd.<unit>` column.
#' @param details an optional plain-text description of the detector (a single
#'   string), shown by [print()]. Defaults to `NA_character_`.
#' @return an `ip_bgrd_detector` object.
#' @examples
#' # a background detector that sets the background to zero everywhere
#' ip_no_bgrd_detector()
#' @export
ip_bgrd_detector <- function(type, detect, details = NA_character_) {
  check_arg(
    type,
    !missing(type) && is_scalar_character(type) && !is.na(type) && nzchar(type),
    "must be a non-empty string"
  )
  check_arg(
    detect,
    !missing(detect) && is_function(detect),
    "must be a function"
  )
  check_arg(details, is_scalar_character(details), "must be a single string")
  structure(
    list(type = type, detect = detect, details = details),
    class = "ip_bgrd_detector"
  )
}

# is x an ip_bgrd_detector?
is_bgrd_detector <- function(x) {
  is(x, "ip_bgrd_detector")
}

# name of the background column matching the intensity column in `peak_traces`
# (e.g. intensity.mV -> bgrd.mV); errors if there is no intensity column
bg_column_name <- function(peak_traces, .call = caller_env()) {
  intensity_col <- grep("^intensity\\.", names(peak_traces), value = TRUE)[1L]
  if (is.na(intensity_col)) {
    cli_abort(
      "the peak traces must include an {.field intensity.*} column for background detection",
      call = .call
    )
  }
  sub("^intensity\\.", "bgrd.", intensity_col)
}

# the inline representation of a background detector ("<type> detector <details>")
# used both by print() and when summarizing it inside a peak detector
format_bgrd_detector <- function(x) {
  if (!is.na(x$details) && nzchar(x$details)) {
    format_inline("{.field {x$type}} detector {x$details}")
  } else {
    format_inline("{.field {x$type}} detector")
  }
}

#' @export
print.ip_bgrd_detector <- function(x, ...) {
  cli_text(format_bgrd_detector(x))
  return(invisible(x))
}

#' @exportS3Method knitr::knit_print
knit_print.ip_bgrd_detector <- function(x, ...) {
  print(x, ...)
}

# built-in background detectors ========

#' @rdname ip_bgrd_detector
#' @description
#' `ip_no_bgrd_detector()` performs no background subtraction: its `detect` sets
#' the background to `0` for every data point (`bgrd.<unit> = 0`).
#'
#' @examples
#' ip_no_bgrd_detector()
#' @export
ip_no_bgrd_detector <- function() {
  ip_bgrd_detector(
    type = "no background",
    detect = function(traces, peak_traces) {
      peak_traces[[bg_column_name(peak_traces)]] <- 0
      peak_traces
    }
  )
}

#' @rdname ip_bgrd_detector
#' @description
#' `ip_individual_bgrd_detector()` determines an individual background for each
#' peak (and each mass) from the baseline leading up to the peak. For every peak
#' it takes the data points in the chosen history window ending at (and
#' including) the peak's first point, optionally smooths them, and applies `func`
#' (the minimum by default) to get a single background value per peak. The
#' history window is given either as a number of seconds (`history.s`) or of data
#' points (`history.pts`); supply exactly one. `smooth_coefficients` is an
#' odd-length vector of (typically normalized) weights for a centered moving
#' average (the default `1` means no smoothing); the extra points needed on each
#' side to smooth the edge of the window are included automatically.
#'
#' @param history.s,history.pts the background window ending at each peak's first
#'   point, as a number of seconds (`history.s`) or of data points
#'   (`history.pts`). Supply exactly one.
#' @param smooth_coefficients an odd-length numeric vector of weights for a
#'   centered moving-average smoothing of the background window. The default `1`
#'   applies no smoothing.
#' @param func a function applied to the (smoothed) background window to produce
#'   the background value, e.g. [min] (the default) or [max].
#' @examples
#' ip_individual_bgrd_detector(history.s = 5)
#' ip_individual_bgrd_detector(history.pts = 25, func = min)
#' @export
ip_individual_bgrd_detector <- function(
  history.s = NA_real_,
  history.pts = NA_real_,
  smooth_coefficients = 1,
  func = min
) {
  # exactly one of history.s / history.pts must be supplied
  has_s <- length(history.s) == 1L && !is.na(history.s)
  has_pts <- length(history.pts) == 1L && !is.na(history.pts)
  if (has_s && has_pts) {
    cli_abort(
      "supply only one of {.field history.s} or {.field history.pts}, not both"
    )
  }
  if (!has_s && !has_pts) {
    cli_abort("must supply one of {.field history.s} or {.field history.pts}")
  }
  if (has_s) {
    check_arg(
      history.s,
      is.numeric(history.s) && length(history.s) == 1L && history.s > 0,
      "must be a single positive number"
    )
  } else {
    check_arg(
      history.pts,
      is_scalar_integerish(history.pts) && history.pts > 0,
      "must be a single positive integer"
    )
    history.pts <- as.integer(history.pts)
  }
  check_arg(
    smooth_coefficients,
    is.numeric(smooth_coefficients) && length(smooth_coefficients) %% 2L == 1L,
    "must be a numeric vector of odd length (1 = no smoothing)"
  )
  check_arg(func, is_function(func), "must be a function")

  # plain-text summary, e.g. "min of last 5 seconds before the peak with 5-point
  # smoothing" (the type "individual background" is added by print())
  history_label <- if (has_s) {
    format_inline("{.strong {history.s} seconds}")
  } else {
    format_inline("{.strong {history.pts} points}")
  }
  smoothing_label <- if (length(smooth_coefficients) > 1L) {
    format_inline(
      "with {.strong {length(smooth_coefficients)}-point} smoothing"
    )
  } else {
    "without smoothing"
  }
  details <- format_inline(
    "({.emph {.strong {describe_function(func)}} intensity of the last {history_label} before the peak {smoothing_label}})"
  )

  ip_bgrd_detector(
    type = "individual background",
    detect = function(traces, peak_traces) {
      detect_individual_background(
        traces,
        peak_traces,
        history.s = history.s,
        history.pts = history.pts,
        smooth_coefficients = smooth_coefficients,
        func = func
      )
    },
    details = details
  )
}

#' @rdname ip_bgrd_detector
#' @description
#' `ip_isodat_default_background()` is the Isodat-default individual background:
#' the minimum of the smoothed baseline over the last `history.pts` data points
#' up to each peak (5-point smoothing).
#'
#' @examples
#' ip_isodat_default_background()
#' @export
ip_isodat_default_background <- function(
  history.pts = 25,
  smooth_coefficients = c(0.16, 0.22, 0.24, 0.22, 0.16)
) {
  ip_individual_bgrd_detector(
    history.pts = history.pts,
    smooth_coefficients = smooth_coefficients,
    func = min
  )
}

# describe a (background) function for the detector summary
describe_function <- function(f) {
  if (identical(f, min)) {
    "min"
  } else if (identical(f, max)) {
    "max"
  } else if (identical(f, mean)) {
    "mean"
  } else {
    "function"
  }
}

# worker for the individual background detector. for each (mass, peak) in
# `peak_traces` it computes one background value from that mass's full trace (in
# `traces`): the data points in the history window ending at (and including) the
# peak's first point are smoothed (centered moving average) and reduced with
# `func`. returns `peak_traces` with the background column (bgrd.<unit>) added.
detect_individual_background <- function(
  traces,
  peak_traces,
  history.s,
  history.pts,
  smooth_coefficients,
  func
) {
  intensity_col <- grep("^intensity\\.", names(peak_traces), value = TRUE)[1L]
  if (is.na(intensity_col)) {
    cli_abort(
      "the peak traces must include an {.field intensity.*} column for background detection"
    )
  }
  bg_col <- sub("^intensity\\.", "bgrd.", intensity_col)
  if (nrow(peak_traces) == 0L) {
    peak_traces[[bg_col]] <- numeric(0)
    return(peak_traces)
  }

  # the first (earliest) point of each peak, per mass
  peak_starts <- peak_traces |>
    dplyr::summarize(
      .by = c("mass", "peak"),
      start_tp = .data$tp[which.min(.data$time.s)],
      start_time = min(.data$time.s)
    )

  # smooth each mass's full trace once; the centered average automatically pulls
  # in the points on each side needed to smooth the edges of the history window
  smoothed_by_mass <- unique(peak_starts$mass) |>
    rlang::set_names() |>
    purrr::map(function(m) {
      mt <- traces |>
        dplyr::filter(as.character(.data$mass) == as.character(m)) |>
        dplyr::arrange(.data$time.s)
      list(
        tp = mt$tp,
        time.s = mt$time.s,
        smoothed = as.numeric(
          stats::filter(mt[[intensity_col]], smooth_coefficients, sides = 2)
        )
      )
    })

  use_pts <- !is.na(history.pts)
  peak_starts[[bg_col]] <- purrr::pmap_dbl(
    peak_starts[c("mass", "start_tp", "start_time")],
    function(mass, start_tp, start_time) {
      s <- smoothed_by_mass[[as.character(mass)]]
      start_pos <- match(start_tp, s$tp)
      if (is.na(start_pos)) {
        return(NA_real_)
      }
      # the history window ends at (includes) the peak's first point
      core <- if (use_pts) {
        seq.int(start_pos - history.pts + 1L, start_pos)
      } else {
        which(s$time.s >= start_time - history.s & s$time.s <= start_time)
      }
      core <- core[core >= 1L & core <= length(s$tp)]
      vals <- s$smoothed[core]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0L) {
        return(NA_real_)
      }
      func(vals)
    }
  )

  dplyr::left_join(
    peak_traces,
    peak_starts[c("mass", "peak", bg_col)],
    by = c("mass", "peak")
  )
}
