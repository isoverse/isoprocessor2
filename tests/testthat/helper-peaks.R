# shared test helpers for peak detection (auto-loaded by testthat before tests)

# minimal fake ir_aggregated_data (a named list of tibbles with the class); adds
# a `tp` time-point column (always present in real isoreader2 data) if missing
fake_aggregated <- function(traces) {
  if (!"tp" %in% names(traces)) {
    traces$tp <- seq_len(nrow(traces))
  }
  structure(list(traces = traces), class = "ir_aggregated_data")
}

# capture cli output, returning list(result, msgs); a wide width keeps each
# reported message on a single line so it can be matched in full; `...` is passed
# on to ip_detect_peaks()
run_detect <- function(agg, det, ...) {
  result <- NULL
  msgs <- withr::with_options(
    list(cli.width = 1000),
    cli::cli_fmt(result <- ip_detect_peaks(agg, det, ...))
  )
  list(result = result, msgs = cli::ansi_strip(msgs))
}

# a detect function that marks every row as belonging to `n_peaks` peaks and adds
# the columns a complete detector result must include (peak + detection_mass). it
# does NOT add a background column - that is the job of the detector's bundled
# background detector, run later by ip_detect_peaks()
mark_peaks <- function(n_peaks = 1L) {
  function(traces) {
    if (nrow(traces) == 0L) {
      return(NULL)
    }
    out <- dplyr::mutate(
      traces,
      peak = rep(seq_len(n_peaks), length.out = dplyr::n())
    )
    out$detection_mass <- TRUE
    out
  }
}

# build a slope-based detector with all required args, overriding via ...
make_slope_detector <- function(...) {
  args <- utils::modifyList(
    list(
      species = "CO2",
      detection_mass = min,
      start_slope.mV_s = 1,
      top_slope.mV_s = 0.05,
      end_slope.mV_s = 1,
      slope_window = 5L,
      slope_window_shift = 1L,
      max_peak_width.s = 180,
      peak_end_max_height.pct = 20,
      bgrd_detector = ip_no_bgrd_detector()
    ),
    list(...)
  )
  do.call(ip_slope_based_peak_detector, args)
}

# the summary is built with cli markup; strip styling before matching content
detector_details <- function(d) cli::ansi_strip(d$details[[1]])

# synthetic single-mass trace with gaussian peaks at the given centers/heights
synth_trace <- function(
  centers,
  heights,
  sigma = 12,
  mass = "44",
  t = seq(0, 600)
) {
  y <- rowSums(mapply(
    function(c, h) h * exp(-((t - c)^2) / (2 * sigma^2)),
    centers,
    heights
  ))
  tibble::tibble(
    species = "CO2",
    mass = mass,
    tp = seq_along(t),
    time.s = t,
    intensity.mV = y
  )
}

# the isodat default detection parameters, for direct worker calls. the worker no
# longer runs a background detector (that happens later, in ip_detect_peaks()), so
# there is no bgrd_detector here
isodat_args <- list(
  start_slope = 0.2,
  top_slope = 0.05,
  end_slope = 0.4,
  unit = "mV_s",
  slope_window = 5L,
  slope_window_shift = 1L,
  max_peak_width.s = 180,
  peak_end_max_height.pct = 50
)

run_slope_detect <- function(traces, detection_mass = min, ...) {
  do.call(
    detect_slope_based_peaks,
    utils::modifyList(
      c(list(traces = traces, detection_mass = detection_mass), isodat_args),
      list(...)
    )
  )
}
