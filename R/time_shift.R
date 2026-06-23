# time shift detector class ========

#' Time shift detectors
#'
#' Time shift detectors describe how the per-mass time offset of a peak (relative
#' to its detection mass) is determined. They are supplied to [ip_detect_peaks()]
#' and run after the peaks (and their rectangularity) have been determined, adding
#' a `time_shift.s` column to the `peaks` dataset. Use the generic
#' `ip_time_shift_detector()` to build a custom one, or the built-in
#' `ip_no_time_shift_detector()` and `ip_parabolic_time_shift_detector()`.
#'
#' `ip_time_shift_detector()` creates an `ip_time_shift_detector` object: like an
#' [ip_bgrd_detector], it is a simple single object (a list) with three elements:
#'
#' - `type` (character): the kind of time shift detector (shown by [print()]).
#' - `detect` (function): a function `function(traces, peaks)` that returns
#'   `peaks` with an added `time_shift.s` column.
#' - `details` (character): an optional plain-text summary of the detector's
#'   parameters (shown by [print()]).
#'
#' @param type the kind of time shift detector (a single non-empty string).
#' @param detect a function `function(traces, peaks)` returning `peaks` with an
#'   added `time_shift.s` column.
#' @param details an optional plain-text description of the detector (a single
#'   string), shown by [print()]. Defaults to `NA_character_`.
#' @return an `ip_time_shift_detector` object.
#' @examples
#' ip_no_time_shift_detector()
#' ip_parabolic_time_shift_detector(apex_window = 9)
#' @export
ip_time_shift_detector <- function(type, detect, details = NA_character_) {
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
    class = "ip_time_shift_detector"
  )
}

# is x an ip_time_shift_detector?
is_time_shift_detector <- function(x) {
  is(x, "ip_time_shift_detector")
}

# the inline representation of a time shift detector ("<type> detector <details>",
# mirroring format_bgrd_detector) used both by print() and when summarizing it
# inside a peak workflow
format_time_shift_detector <- function(x) {
  if (!is.na(x$details) && nzchar(x$details)) {
    format_inline("{.field {x$type}} detector {x$details}")
  } else {
    format_inline("{.field {x$type}} detector")
  }
}

#' @export
print.ip_time_shift_detector <- function(x, ...) {
  cli_text(format_time_shift_detector(x))
  return(invisible(x))
}

#' @exportS3Method knitr::knit_print
knit_print.ip_time_shift_detector <- function(x, ...) {
  print(x, ...)
}

# built-in time shift detectors ========

#' @rdname ip_time_shift_detector
#' @description
#' `ip_no_time_shift_detector()` applies no time shift: its `detect` sets
#' `time_shift.s = 0` for every peak.
#'
#' @examples
#' ip_no_time_shift_detector()
#' @export
ip_no_time_shift_detector <- function() {
  ip_time_shift_detector(
    type = "no time shift",
    detect = function(traces, peaks) {
      peaks[["time_shift.s"]] <- if (nrow(peaks) == 0L) numeric(0) else 0
      peaks
    }
  )
}

#' @rdname ip_time_shift_detector
#' @description
#' `ip_parabolic_time_shift_detector()` estimates each (non-rectangular) peak's
#' sub-sample apex time on every mass trace with a parabolic fit over `apex_window`
#' points around the peak's apex index, and reports `time_shift.s` as the offset
#' of each mass's apex from the detection mass's apex (the detection mass itself
#' therefore gets `0`). Rectangular peaks always get `time_shift.s = 0`.
#'
#' @param apex_window the number of points (an odd integer >= 3, default `9`) used
#'   for the parabolic apex fit around each peak's apex.
#' @examples
#' ip_parabolic_time_shift_detector(apex_window = 9)
#' @export
ip_parabolic_time_shift_detector <- function(apex_window = 9L) {
  check_arg(
    apex_window,
    is_scalar_integerish(apex_window) &&
      apex_window >= 3L &&
      apex_window %% 2L == 1L,
    "must be a single odd integer >= 3"
  )
  apex_window <- as.integer(apex_window)
  ip_time_shift_detector(
    type = "parabolic time shift",
    detect = function(traces, peaks) {
      detect_parabolic_time_shift(traces, peaks, apex_window = apex_window)
    },
    details = format_inline(
      "({.strong {apex_window} point} parabolic apex window)"
    )
  )
}

# worker for the parabolic time shift detector. for each non-rectangular peak it
# finds the sub-sample apex time on each mass trace (a parabolic fit of width
# `apex_window` around the apex index, via `.parabolic_center`) and reports
# `time_shift.s` = mass apex - detection-mass apex. rectangular peaks (and the
# detection mass itself) get 0.
detect_parabolic_time_shift <- function(traces, peaks, apex_window) {
  if (nrow(peaks) == 0L) {
    peaks[["time_shift.s"]] <- numeric(0)
    return(peaks)
  }
  intensity_col <- grep("^intensity\\.", names(traces), value = TRUE)[1L]
  if (is.na(intensity_col)) {
    cli_abort(
      "the {.field traces} must include an {.field intensity.*} column for parabolic time shift detection"
    )
  }
  if (!"tp" %in% names(traces)) {
    cli_abort(
      "the {.field traces} must include a {.field tp} column for parabolic time shift detection"
    )
  }
  for (col in c("apex.idx", "detection_mass", "rectangular")) {
    if (!col %in% names(peaks)) {
      cli_abort(
        "the {.field peaks} must include a {.field {col}} column for parabolic time shift detection"
      )
    }
  }

  # split the (time-ordered) traces by mass trace for quick per-peak lookup
  trace_keys <- intersect(
    c("uidx", "analysis", "species", "mass"),
    names(traces)
  )
  traces <- dplyr::arrange(traces, .data$time.s)
  trace_key <- if (length(trace_keys) > 0L) {
    do.call(paste, c(traces[trace_keys], sep = "\r"))
  } else {
    rep("", nrow(traces))
  }
  trace_list <- split(traces, trace_key)
  peak_key <- if (length(trace_keys) > 0L) {
    do.call(paste, c(peaks[trace_keys], sep = "\r"))
  } else {
    rep("", nrow(peaks))
  }

  # parabolic apex time for each peaks row (only the non-rectangular ones need it)
  peaks$.apex_center <- vapply(
    seq_len(nrow(peaks)),
    function(i) {
      if (isTRUE(peaks$rectangular[i])) {
        return(NA_real_)
      }
      tr <- trace_list[[peak_key[i]]]
      if (is.null(tr)) {
        return(NA_real_)
      }
      pos <- match(peaks$apex.idx[i], tr$tp)
      if (is.na(pos)) {
        return(NA_real_)
      }
      .parabolic_center(tr$time.s, tr[[intensity_col]], pos, apex_window)
    },
    numeric(1)
  )

  # the detection-mass apex center per peak, joined back onto every mass
  peak_keys <- intersect(c("uidx", "analysis", "species", "peak"), names(peaks))
  det_center <- dplyr::filter(peaks, .data$detection_mass)
  det_center <- det_center[c(peak_keys, ".apex_center")]
  names(det_center)[length(det_center)] <- ".det_center"
  peaks <- dplyr::left_join(peaks, det_center, by = peak_keys)

  peaks[["time_shift.s"]] <- dplyr::if_else(
    peaks$rectangular | is.na(peaks$.apex_center) | is.na(peaks$.det_center),
    0,
    peaks$.apex_center - peaks$.det_center
  )
  peaks$.apex_center <- NULL
  peaks$.det_center <- NULL
  peaks
}
