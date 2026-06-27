# slope-based peak detector ========

#' @rdname ip_peak_detector
#' @description
#' `ip_slope_based_peak_detector()` creates a detector with `type =
#' "slope-based"` that locates peak boundaries from rolling-window regression
#' slopes of the intensity trace (a background -> peak-start -> peak-top state
#' machine). The slope thresholds are supplied in the unit matching the intensity
#' trace: either mV/s (voltage) or pA/s (current); for each of the start/top/end
#' thresholds supply exactly one of the `.mV_s` / `.pA_s` arguments, all using the
#' same unit. Its detection parameters are deliberately required (no defaults) so
#' each detector is explicit.
#'
#' @param detection_mass which mass/trace to run detection on: a function applied
#'   to the available masses to pick one (e.g. `min` or `max`), a specific mass
#'   number, or a mass label (string).
#' @param start_slope.mV_s,start_slope.pA_s slope threshold (in mV/s or pA/s) for
#'   detecting the rising edge / start of a peak. Supply exactly one.
#' @param top_slope.mV_s,top_slope.pA_s slope threshold (in mV/s or pA/s) for
#'   detecting the peak top (slope zero-crossing). Supply exactly one, in the same
#'   unit as `start_slope`.
#' @param end_slope.mV_s,end_slope.pA_s slope threshold (in mV/s or pA/s) for
#'   detecting the falling edge / end of a peak. Supply exactly one, in the same
#'   unit as `start_slope`.
#' @param slope_window the rolling-regression window size (odd integer >= 3) used
#'   to compute slopes (the regression window for the rolling slope).
#' @param slope_window_shift how far to shift the rolling-regression window
#'   relative to the target point (0 = centered, positive = forward).
#' @param max_peak_width.s the maximum allowed peak width in seconds; peaks
#'   exceeding this are force-ended.
#' @param peak_resolution.pct the peak resolution (in %, matching the
#'   Isodat/Qtegra "Peak Resolution"): the minimum drop from a peak's top - as a
#'   percentage of its (background-subtracted) height - required before the peak
#'   may end or before a shoulder is split off as a separate peak. The end/split
#'   logic only engages once the signal has dropped to at or below
#'   `(100 - peak_resolution.pct)%` of the peak height, so a *higher* resolution
#'   demands a deeper valley between two maxima for them to be counted as two
#'   peaks (and otherwise merges them into one).
#' @param min_height.mV,min_height.pA an optional minimum peak height (in mV or
#'   pA; supply at most one, matching the slope unit's signal kind). Peaks whose
#'   detection-mass height above the peak's own start point is below this are
#'   dropped and the rest renumbered. (The true background, determined later by
#'   [ip_detect_peaks()] after the time shift, is not yet available at detection.)
#' @param bgrd_detector an [ip_bgrd_detector] (e.g. [ip_no_bgrd_detector()])
#'   bundled with the detector (stored in its `bgrd` column). It is not run during
#'   detection but by [ip_detect_peaks()], after the time shift, to add a
#'   `bgrd.<unit>` background column to the detected peak traces.
#' @examples
#' ip_slope_based_peak_detector(
#'   "CO2",
#'   detection_mass = min,
#'   start_slope.mV_s = 0.2,
#'   top_slope.mV_s = 0.05,
#'   end_slope.mV_s = 0.4,
#'   slope_window = 5L,
#'   slope_window_shift = 1L,
#'   max_peak_width.s = 180,
#'   peak_resolution.pct = 50,
#'   bgrd_detector = ip_no_bgrd_detector()
#' )
#' @export
ip_slope_based_peak_detector <- function(
  species = NA_character_,
  start.s = 60 * start.min,
  stop.s = 60 * stop.min,
  start.min = 0,
  stop.min = Inf,
  detection_mass,
  start_slope.mV_s = NA_real_,
  start_slope.pA_s = NA_real_,
  top_slope.mV_s = NA_real_,
  top_slope.pA_s = NA_real_,
  end_slope.mV_s = NA_real_,
  end_slope.pA_s = NA_real_,
  slope_window,
  slope_window_shift,
  max_peak_width.s,
  peak_resolution.pct,
  min_height.mV = NA_real_,
  min_height.pA = NA_real_,
  bgrd_detector
) {
  # detection mass selector: a function (e.g. min/max), a mass number, or a label
  check_arg(
    detection_mass,
    !missing(detection_mass) &&
      (is_function(detection_mass) ||
        (length(detection_mass) == 1L &&
          !is.na(detection_mass) &&
          (is.numeric(detection_mass) || is.character(detection_mass)))),
    "must be a function, a single number, or a single string"
  )

  # resolve each slope threshold to a single value + unit, enforcing that
  # exactly one of the .mV_s / .pA_s pair is supplied
  start <- resolve_slope_unit(start_slope.mV_s, start_slope.pA_s, "start_slope")
  top <- resolve_slope_unit(top_slope.mV_s, top_slope.pA_s, "top_slope")
  end <- resolve_slope_unit(end_slope.mV_s, end_slope.pA_s, "end_slope")

  # all three thresholds must use the same unit
  if (length(unique(c(start$unit, top$unit, end$unit))) > 1L) {
    cli_abort(c(
      "{.field start_slope}, {.field top_slope} and {.field end_slope} must all use the same unit (mV/s or pA/s), not a mix",
      "x" = "got start = {.emph {start$unit}}, top = {.emph {top$unit}}, end = {.emph {end$unit}}"
    ))
  }
  unit <- start$unit
  start_slope <- start$value
  top_slope <- top$value
  end_slope <- end$value

  # validate the remaining type-specific parameters (all required, no defaults)
  check_arg(
    slope_window,
    !missing(slope_window) &&
      is_scalar_integerish(slope_window) &&
      slope_window >= 3L &&
      slope_window %% 2L == 1L,
    "must be a single odd integer >= 3"
  )
  slope_window <- as.integer(slope_window)
  check_arg(
    slope_window_shift,
    !missing(slope_window_shift) &&
      is_scalar_integerish(slope_window_shift) &&
      abs(slope_window_shift) <= slope_window %/% 2L,
    "must be a single integer with magnitude <= slope_window %/% 2"
  )
  slope_window_shift <- as.integer(slope_window_shift)
  check_arg(
    max_peak_width.s,
    !missing(max_peak_width.s) &&
      is.numeric(max_peak_width.s) &&
      length(max_peak_width.s) == 1L &&
      !is.na(max_peak_width.s) &&
      max_peak_width.s > 0,
    "must be a single positive number"
  )
  check_arg(
    peak_resolution.pct,
    !missing(peak_resolution.pct) &&
      is.numeric(peak_resolution.pct) &&
      length(peak_resolution.pct) == 1L &&
      !is.na(peak_resolution.pct) &&
      peak_resolution.pct >= 0 &&
      peak_resolution.pct <= 100,
    "must be a single number between 0 and 100"
  )
  check_arg(
    bgrd_detector,
    !missing(bgrd_detector) && is_bgrd_detector(bgrd_detector),
    format_inline(
      "must be an {.cls ip_bgrd_detector}, e.g. {.emph ip_individual_bgrd_detector()}"
    )
  )
  # optional min height (at most one of mV/pA); NULL means no height filtering
  min_height <- resolve_min_height(min_height.mV, min_height.pA)

  # detection step — the resolved parameters are captured via this closure's
  # enclosing environment and handed to the slope-based detection worker. note
  # that the background detector is NOT run here: it is bundled with the peak
  # detector (the `bgrd` column) and run by ip_detect_peaks() after the time shift.
  detect <- function(traces) {
    detect_slope_based_peaks(
      traces,
      detection_mass = detection_mass,
      start_slope = start_slope,
      top_slope = top_slope,
      end_slope = end_slope,
      unit = unit,
      slope_window = slope_window,
      slope_window_shift = slope_window_shift,
      max_peak_width.s = max_peak_width.s,
      peak_resolution.pct = peak_resolution.pct,
      min_height = min_height
    )
  }

  # summary of the supplied parameters, with the values emphasized so they are
  # easy to spot in the printed detector
  unit_label <- if (unit == "mV_s") "mV/s" else "pA/s"
  min_height_label <- if (!is.null(min_height)) {
    format_inline(
      "{col_silver('min height')} = {.strong {min_height$value} {min_height$unit}}; "
    )
  } else {
    ""
  }
  details <- format_inline(
    "{col_silver('detection mass')} = {.strong {describe_detection_mass(detection_mass)}}; ",
    "{col_silver('start/top/end slope')} = {.strong {start_slope}/{top_slope}/{end_slope} {unit_label}}; ",
    "{col_silver('slope window')} = {.strong {slope_window}} ({col_silver('shift')} {.strong {slope_window_shift}}); ",
    "{col_silver('max width')} = {.strong {max_peak_width.s} s}; ",
    "{col_silver('peak resolution')} = {.strong {peak_resolution.pct}%}; ",
    "{min_height_label}",
    "{col_silver('background')} = {format_bgrd_detector(bgrd_detector)}"
  )

  ip_peak_detector(
    type = "slope-based",
    species = species,
    start.s = start.s,
    stop.s = stop.s,
    detect = detect,
    bgrd = bgrd_detector,
    details = details
  )
}

# resolve a slope threshold given its mV/s and pA/s alternatives: exactly one
# must be supplied (non-NA). returns list(value=, unit=) where unit is "mV_s" or
# "pA_s"; errors if both or neither are supplied, or a supplied value is invalid.
resolve_slope_unit <- function(mV_s, pA_s, arg, .call = caller_env()) {
  validate <- function(z, suffix) {
    # the NA_real_ default means "not supplied"
    if (length(z) == 1L && is.na(z)) {
      return(FALSE)
    }
    if (!is.numeric(z) || length(z) != 1L) {
      cli_abort("{.field {arg}.{suffix}} must be a single number", call = .call)
    }
    TRUE
  }
  has_mV <- validate(mV_s, "mV_s")
  has_pA <- validate(pA_s, "pA_s")
  if (has_mV && has_pA) {
    cli_abort(
      "supply only one of {.field {arg}.mV_s} or {.field {arg}.pA_s}, not both",
      call = .call
    )
  }
  if (!has_mV && !has_pA) {
    cli_abort(
      "must supply one of {.field {arg}.mV_s} or {.field {arg}.pA_s}",
      call = .call
    )
  }
  if (has_mV) {
    list(value = as.double(mV_s), unit = "mV_s")
  } else {
    list(value = as.double(pA_s), unit = "pA_s")
  }
}

# resolve the optional min height given its mV / pA alternatives: at most one may
# be supplied. returns NULL if neither, else list(value=, unit="mV"/"pA").
resolve_min_height <- function(
  min_height.mV,
  min_height.pA,
  .call = caller_env()
) {
  validate <- function(z, suffix) {
    if (length(z) == 1L && is.na(z)) {
      return(FALSE) # the NA_real_ default means "not supplied"
    }
    if (!is.numeric(z) || length(z) != 1L || z <= 0) {
      cli_abort(
        "{.field {paste0('min_height.', suffix)}} must be a single positive number",
        call = .call
      )
    }
    TRUE
  }
  has_mV <- validate(min_height.mV, "mV")
  has_pA <- validate(min_height.pA, "pA")
  if (has_mV && has_pA) {
    cli_abort(
      "supply only one of {.field min_height.mV} or {.field min_height.pA}, not both",
      call = .call
    )
  }
  if (!has_mV && !has_pA) {
    return(NULL)
  }
  if (has_mV) {
    list(value = as.double(min_height.mV), unit = "mV")
  } else {
    list(value = as.double(min_height.pA), unit = "pA")
  }
}

# describe the detection_mass selector for the detector summary
describe_detection_mass <- function(x) {
  if (is.function(x)) {
    if (identical(x, min)) {
      return("min")
    }
    if (identical(x, max)) {
      return("max")
    }
    return("function")
  }
  as.character(x)
}

#' @rdname ip_peak_detector
#' @description
#' `ip_isodat_default_detector()` is a convenience wrapper around
#' `ip_slope_based_peak_detector()` preconfigured with the slope-based settings
#' Isodat uses for continuous-flow data: mV/s slope thresholds of 0.2 (start),
#' 0.05 (top) and 0.4 (end), the minimum mass as the detection trace, a 5-point
#' slope window (shifted forward by 1), a 180 s maximum peak width, a 50% peak
#' resolution (the Isodat default), a 50 mV minimum height, and the Isodat
#' individual background ([ip_isodat_default_background()]).
#'
#' @examples
#' ip_isodat_default_detector("CO2")
#' @export
ip_isodat_default_detector <- function(
  species = NA_character_,
  start.s = 60 * start.min,
  stop.s = 60 * stop.min,
  start.min = 0,
  stop.min = Inf,
  detection_mass = min,
  start_slope.mV_s = 0.2,
  top_slope.mV_s = 0.05,
  end_slope.mV_s = 0.4,
  max_peak_width.s = 180,
  peak_resolution.pct = 50,
  min_height.mV = 50,
  bgrd_detector = ip_isodat_default_background()
) {
  ip_slope_based_peak_detector(
    species = species,
    start.s = start.s,
    stop.s = stop.s,
    detection_mass = detection_mass,
    start_slope.mV_s = start_slope.mV_s,
    top_slope.mV_s = top_slope.mV_s,
    end_slope.mV_s = end_slope.mV_s,
    slope_window = 5L,
    slope_window_shift = 1L,
    max_peak_width.s = max_peak_width.s,
    peak_resolution.pct = peak_resolution.pct,
    min_height.mV = min_height.mV,
    bgrd_detector = bgrd_detector
  )
}

# detection helpers ========

# calculate the rolling linear-regression (OLS) slope of `y` against `x` over a
# window of `window_size` points, using slider::slide2_dbl (no explicit loop).
# `window_size` must be an odd integer >= 3. `window_shift` shifts the regression
# window relative to the target point: 0 keeps it centered, positive values shift
# it forward (more points after the target), negative values shift it backward.
# for window_size = 2k+1 the window spans `k - window_shift` points before and
# `k + window_shift` points after the target, so `abs(window_shift)` must be <= k.
# `window_shift = 1` reproduces the original slope_n behaviour (the regression
# effectively centered on the next point). returns a numeric vector the same
# length as the inputs, with NA wherever the full window does not fit (near the
# edges). `x` is assumed to be the (strictly increasing) time axis, so its window
# variance is always > 0.
calculate_rolling_slope <- function(x, y, window_size, window_shift) {
  check_arg(x, is.numeric(x), "must be numeric")
  check_arg(y, is.numeric(y), "must be numeric")
  check_arg(
    window_size,
    is_scalar_integerish(window_size) &&
      window_size >= 3L &&
      window_size %% 2L == 1L,
    "must be a single odd integer >= 3"
  )
  half <- as.integer(window_size) %/% 2L
  check_arg(
    window_shift,
    !missing(window_shift) && is_scalar_integerish(window_shift),
    "must be a single integer"
  )
  if (abs(window_shift) > half) {
    cli_abort(c(
      "{.field window_shift} ({window_shift}) is too large for {.field window_size} ({window_size})",
      "i" = "{.field window_shift} must satisfy {.code abs(window_shift) <= window_size %/% 2} (here, <= {half})"
    ))
  }
  if (length(x) != length(y)) {
    cli_abort(
      "{.field x} (length {length(x)}) and {.field y} (length {length(y)}) must have the same length"
    )
  }

  # window shifted relative to the target point: for window_size = 2k+1 the
  # window spans (k - window_shift) points before and (k + window_shift) after
  window_shift <- as.integer(window_shift)
  before <- half - window_shift
  after <- half + window_shift

  slider::slide2_dbl(
    x,
    y,
    .f = function(xi, yi) {
      mx <- mean(xi)
      sum((xi - mx) * (yi - mean(yi))) / sum((xi - mx)^2)
    },
    .before = before,
    .after = after,
    .complete = TRUE
  )
}

# slope-based detection ========

# Parse a signal unit like "mV", "V", "nA", "pA" into its base ("V" for voltage,
# "A" for current) and the metric-prefix scale factor relative to that base (e.g.
# "mV" -> base "V", factor 1e-3). Returns NULL for anything that is not a
# metric-prefixed volt or ampere (e.g. "cps").
parse_signal_unit <- function(unit_str) {
  base <- if (endsWith(unit_str, "V")) {
    "V"
  } else if (endsWith(unit_str, "A")) {
    "A"
  } else {
    return(NULL)
  }
  prefix <- substr(unit_str, 1L, nchar(unit_str) - 1L)
  idx <- match(prefix, names(.metric_prefixes))
  if (is.na(idx)) {
    return(NULL)
  }
  list(base = base, factor = .metric_prefixes[[idx]])
}

# Slope-based peak detection worker.
#
# This is the engine behind `ip_slope_based_peak_detector()`'s `detect` step. It
# is given a `traces` tibble for a *single species* (already restricted to the
# detector's time window by `ip_detect_peaks()`) and returns the trace rows that
# belong to detected peaks, each tagged with an integer `peak` column (1, 2, 3,
# ...) and a logical `detection_mass` column. Returns `NULL` if nothing is
# detected. Peaks are *found* on the detection mass only, but the rows *returned*
# cover **all** masses (the detection mass plus the others) at the detected time
# points, so downstream steps can integrate every mass over the same peaks.
#
# Steps:
#   1. Pick the detection mass. If `detection_mass` is a function it is applied to
#      the `mass` column to choose a mass (e.g. `min`/`max`); otherwise the mass
#      is matched by `as.character()`.
#   2. Detection runs independently per analysis (grouping by whichever of `uidx`
#      / `analysis` are present) so separate time series are not mixed.
#   3. Within each group `find_slope_peak_ranges()` returns the start/end index of
#      every peak in the detection trace. Each peak's time points (`tp`) are then
#      joined back to *all* masses in that group, so every mass contributes its
#      rows for that peak; a `tp` shared by two adjacent peaks (a start/end
#      boundary) is emitted once per peak, for every mass.
#
# The slope thresholds are given in mV/s or pA/s (`unit`). The trace's intensity
# column may use any metric prefix of the *same* kind of signal (volt or ampere):
# the thresholds are converted to the intensity's unit (e.g. mV/s -> V/s, pA/s ->
# nA/s). A mismatch in kind (e.g. mV/s thresholds against an ampere intensity, or
# a non-volt/ampere unit such as cps) is an error.
#
# @param traces tibble of traces for one species (must contain `mass`, `time.s`,
#   and an `intensity.*` column of the same signal kind as `unit`).
# @param detection_mass function / number / string selecting the detection mass.
# @param start_slope,top_slope,end_slope slope thresholds (in the `unit`'s
#   per-second unit) for the rising edge, peak top and falling edge.
# @param unit `"mV_s"` or `"pA_s"`; the signal kind the thresholds are expressed
#   in (voltage or current).
# @param slope_window,slope_window_shift passed to `calculate_rolling_slope()`.
# @param max_peak_width.s force-end peaks wider than this (s).
# @param peak_resolution.pct peak resolution (%): the end/shoulder-split logic is
#   only checked once the signal has dropped to <= (100 - peak_resolution.pct)% of
#   the peak top (i.e. dropped by >= peak_resolution.pct% from the top).
# @param min_height NULL, or list(value=, unit="mV"/"pA"): peaks whose
#   detection-mass height above the peak's own start point (max(intensity) minus
#   the intensity at the first point) is below this are dropped and the remaining
#   peaks renumbered. (The true background is not yet known at this stage - it is
#   determined later, by ip_detect_peaks(), after the time shift.) The number
#   dropped is attached as the "filter_info" attribute of the result.
# @return a tibble of all masses' rows belonging to detected peaks, with an
#   integer `peak` column and a logical `detection_mass` column (placed right
#   after `mass`), but no background column, or `NULL` if no peaks were found.
detect_slope_based_peaks <- function(
  traces,
  detection_mass,
  start_slope,
  top_slope,
  end_slope,
  unit,
  slope_window,
  slope_window_shift,
  max_peak_width.s,
  peak_resolution.pct,
  min_height = NULL
) {
  if (nrow(traces) == 0L) {
    return(NULL)
  }

  # find the intensity column and reconcile its unit with the slope unit: the
  # signal kind must match (voltage vs current), but a different metric prefix is
  # fine - the slope thresholds are simply converted to the intensity's unit
  intensity_cols <- grep("^intensity\\.", names(traces), value = TRUE)
  if (length(intensity_cols) == 0L) {
    cli_abort("the traces must include an {.field intensity.*} column")
  }
  intensity_col <- intensity_cols[1L]
  intensity_unit <- sub("^intensity\\.", "", intensity_col)
  slope_unit <- if (unit == "mV_s") "mV" else "pA"
  slope_signal <- parse_signal_unit(slope_unit)
  intensity_signal <- parse_signal_unit(intensity_unit)

  describe_signal <- function(parsed) {
    if (is.null(parsed)) {
      "neither a voltage nor a current"
    } else if (parsed$base == "V") {
      "a voltage"
    } else {
      "a current"
    }
  }
  if (is.null(intensity_signal) || intensity_signal$base != slope_signal$base) {
    cli_abort(c(
      "slope thresholds in {.emph {slope_unit}/s} are not compatible with the trace intensity in {.emph {intensity_unit}}",
      "x" = "{.emph {slope_unit}/s} measures {describe_signal(slope_signal)}, but {.field {intensity_col}} is {describe_signal(intensity_signal)}"
    ))
  }

  # same signal kind, possibly different prefix -> convert thresholds to the
  # intensity unit (e.g. mV/s -> V/s, or pA/s -> nA/s)
  if (intensity_unit != slope_unit) {
    conversion <- slope_signal$factor / intensity_signal$factor
    start_slope <- start_slope * conversion
    top_slope <- top_slope * conversion
    end_slope <- end_slope * conversion
  }

  # reconcile the optional min height the same way (same signal kind required,
  # different prefix converted to the intensity unit)
  min_height_value <- NULL
  if (!is.null(min_height)) {
    mh_signal <- parse_signal_unit(min_height$unit)
    if (is.null(intensity_signal) || intensity_signal$base != mh_signal$base) {
      cli_abort(c(
        "the min height in {.emph {min_height$unit}} is not compatible with the trace intensity in {.emph {intensity_unit}}",
        "x" = "{.emph {min_height$unit}} measures {describe_signal(mh_signal)}, but {.field {intensity_col}} is {describe_signal(intensity_signal)}"
      ))
    }
    min_height_value <- if (intensity_unit != min_height$unit) {
      min_height$value * (mh_signal$factor / intensity_signal$factor)
    } else {
      min_height$value
    }
  }

  if (!"mass" %in% names(traces)) {
    cli_abort(
      "the traces must include a {.field mass} column for slope-based detection"
    )
  }

  # 1. pick the detection mass
  detection_mass_value <- if (is.function(detection_mass)) {
    detection_mass(traces$mass)
  } else {
    detection_mass
  }
  if (!any(as.character(traces$mass) == as.character(detection_mass_value))) {
    return(NULL)
  }

  # 2. detect per analysis (separate time series must not be concatenated), but
  #    keep all masses in each analysis group so they can be returned together
  group_cols <- intersect(c("uidx", "analysis"), names(traces))
  groups <- if (length(group_cols) > 0L) {
    traces |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::group_split()
  } else {
    list(traces)
  }

  # 3. find peaks on the detection trace, then join each peak's time points (tp)
  #    back to all masses in the group
  results <- list()
  peak_offset <- 0L
  for (group in groups) {
    detection_trace <- group |>
      dplyr::filter(
        as.character(.data$mass) == as.character(!!detection_mass_value)
      ) |>
      dplyr::arrange(.data$time.s)
    if (nrow(detection_trace) == 0L) {
      next
    }
    ranges <- find_slope_peak_ranges(
      time = detection_trace$time.s,
      intensity = detection_trace[[intensity_col]],
      start_slope = start_slope,
      top_slope = top_slope,
      end_slope = end_slope,
      slope_window = slope_window,
      slope_window_shift = slope_window_shift,
      max_peak_width.s = max_peak_width.s,
      peak_resolution.pct = peak_resolution.pct
    )
    if (nrow(ranges) == 0L) {
      next
    }

    # each peak's detection-trace time points; a tp shared by two adjacent peaks
    # appears once per peak (so overlapping peaks stay fully represented)
    tp_by_peak <- purrr::map2(
      ranges$start,
      ranges$end,
      ~ detection_trace$tp[seq.int(.x, .y)]
    )
    membership <- tibble(
      tp = unlist(tp_by_peak),
      peak = rep(seq_len(nrow(ranges)) + peak_offset, lengths(tp_by_peak))
    )
    # pull all masses in this group at those time points and flag the detection
    # mass (placed right after the mass column); the background is NOT added here
    group_peak_traces <- dplyr::inner_join(
      group,
      membership,
      by = "tp",
      relationship = "many-to-many"
    ) |>
      dplyr::mutate(
        detection_mass = as.character(.data$mass) ==
          as.character(!!detection_mass_value)
      ) |>
      dplyr::relocate("detection_mass", .after = "mass")

    results[[length(results) + 1L]] <- group_peak_traces
    peak_offset <- peak_offset + nrow(ranges)
  }

  if (length(results) == 0L) {
    return(NULL)
  }

  combined <- dplyr::bind_rows(results) |>
    dplyr::arrange(.data$peak, .data$mass, .data$time.s)

  # drop peaks whose detection-mass height above the peak's own start point is
  # below the min height, then renumber the remaining peaks 1..n. (the true
  # background is determined later, by ip_detect_peaks(), so it cannot be used
  # here; the height above the start point is the rough estimate the state machine
  # itself uses.)
  filtered_out <- 0L
  if (!is.null(min_height_value)) {
    peak_height <- combined |>
      dplyr::filter(.data$detection_mass) |>
      dplyr::summarize(
        .by = "peak",
        height = max(.data[[intensity_col]]) -
          .data[[intensity_col]][which.min(.data$time.s)]
      )
    keep <- peak_height$peak[peak_height$height >= min_height_value]
    filtered_out <- length(unique(combined$peak)) - length(keep)
    combined <- combined |>
      dplyr::filter(.data$peak %in% keep) |>
      dplyr::mutate(peak = match(.data$peak, sort(unique(.data$peak))))
  }

  if (filtered_out > 0L) {
    attr(combined, "filter_info") <- format_inline(
      "{filtered_out} additional below the {col_silver('min height')} threshold filtered out"
    )
  }
  combined
}

# Slope-based peak-boundary state machine.
#
# Given a single, ascending `(time, intensity)` series, walk it once with a
# 3-state machine to locate peak boundaries and return their index ranges as a
# tibble with `start` and `end` columns (indices into `time`/`intensity`).
#
# The walk uses the rolling regression slope (see `calculate_rolling_slope()`).
# Following the original Isodat-style detector, two slopes two valid steps apart
# (`slope1` = older, `slope2` = current) are combined: their sum is compared to
# the start/end thresholds and their individual crossing of `top_slope` marks the
# peak top / shoulders. The states are:
#
#   * BACKGROUND - waiting for a rising edge: when `slope1 + slope2 > 2 *
#     start_slope` a peak starts. Its *first data point* is taken as the peak's
#     background reference (no separate background history is used).
#   * RISING - climbing to the top: when the slope falls back through `top_slope`
#     (`slope1 > top_slope` and `slope2 < top_slope`) the top is reached. If that
#     candidate top is below the peak's first point (it never actually rose) the
#     peak is discarded.
#   * FALLING - descending from the top: the top is updated while the signal keeps
#     climbing; once the (background-subtracted) signal has dropped to at or below
#     `(100 - peak_resolution.pct)%` of the top height (i.e. a drop of at least the
#     peak resolution) and the slope is no longer steeply negative (`slope1 +
#     slope2 > -2 * end_slope`) the peak ends. If instead the slope turns back up
#     through `top_slope` a shoulder is split off: the current peak ends and a new
#     peak begins at the same point. A higher `peak_resolution.pct` therefore
#     demands a deeper valley before two maxima are split into two peaks.
#
# Peaks wider than `max_peak_width.s` are force-ended. A peak still open when the
# series ends is not recorded.
#
# @return tibble with integer `start` and `end` columns (one row per peak).
find_slope_peak_ranges <- function(
  time,
  intensity,
  start_slope,
  top_slope,
  end_slope,
  slope_window,
  slope_window_shift,
  max_peak_width.s,
  peak_resolution.pct
) {
  empty <- tibble(start = integer(0), end = integer(0))
  n <- length(time)
  if (n < 3L) {
    return(empty)
  }

  # rolling slope; only positions with a defined slope take part in the walk
  slopes <- calculate_rolling_slope(
    time,
    intensity,
    window_size = slope_window,
    window_shift = slope_window_shift
  )
  valid <- which(!is.na(slopes))
  if (length(valid) < 3L) {
    return(empty)
  }

  # the fraction of the top height the signal must drop to before the end/shoulder
  # check engages: a peak resolution of R% requires the signal to fall to at or
  # below (100 - R)% of the top (i.e. a drop of at least R%), matching the Isodat /
  # Qtegra "Peak Resolution"
  drop_frac <- (100 - peak_resolution.pct) / 100

  state <- 0L # 0 = background, 1 = rising, 2 = falling
  start_idx <- NA_integer_
  top_height <- -Inf
  bg <- NA_real_ # background reference = peak's first data point
  starts <- integer(0)
  ends <- integer(0)

  for (k in seq.int(3L, length(valid))) {
    slope1 <- slopes[valid[k - 2L]] # older slope
    slope2 <- slopes[valid[k]] # current slope
    ci <- valid[k] # current data point
    cy <- intensity[ci]

    if (state == 0L) {
      # BACKGROUND: rising edge starts a peak
      if (slope1 + slope2 > 2 * start_slope) {
        state <- 1L
        start_idx <- ci
        bg <- intensity[ci]
      }
    } else if (state == 1L) {
      # RISING: look for the top (slope dropping back through top_slope)
      if (time[ci] - time[start_idx] > max_peak_width.s) {
        state <- 0L
        next
      }
      if (slope1 > top_slope && slope2 < top_slope) {
        if (cy - bg < 0) {
          # never rose above its own start -> not a real peak
          state <- 0L
          next
        }
        state <- 2L
        top_height <- cy - bg
      }
    } else {
      # FALLING: track the top, then look for the end
      if (cy - bg > top_height) {
        top_height <- cy - bg
      }
      if (time[ci] - time[start_idx] > max_peak_width.s) {
        starts <- c(starts, start_idx)
        ends <- c(ends, ci)
        state <- 0L
        next
      }
      if (cy - bg <= drop_frac * top_height) {
        if (slope1 + slope2 <= -2 * end_slope) {
          next # still steeply descending
        }
        if (slope1 < top_slope && slope2 > top_slope) {
          # shoulder: end this peak and start the next at the same point
          starts <- c(starts, start_idx)
          ends <- c(ends, ci)
          state <- 1L
          start_idx <- ci
          bg <- intensity[ci]
          top_height <- -Inf
        } else {
          starts <- c(starts, start_idx)
          ends <- c(ends, ci)
          state <- 0L
        }
      }
    }
  }

  tibble(start = starts, end = ends)
}
