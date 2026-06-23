# Extracted from test-peak-detection.R:540

# prequel ----------------------------------------------------------------------
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
      peak_height_pct = 20
    ),
    list(...)
  )
  do.call(ip_slope_based_peak_detector, args)
}
detector_details <- function(d) cli::ansi_strip(d$details[[1]])
fake_aggregated <- function(traces) {
  structure(list(traces = traces), class = "ir_aggregated_data")
}
run_detect <- function(agg, det) {
  result <- NULL
  msgs <- cli::cli_fmt(result <- ip_detect_peaks(agg, det))
  list(result = result, msgs = cli::ansi_strip(msgs))
}
mark_peaks <- function(n_peaks = 1L) {
  function(traces) {
    if (nrow(traces) == 0L) {
      return(NULL)
    }
    dplyr::mutate(traces, peak = rep(seq_len(n_peaks), length.out = dplyr::n()))
  }
}

# test -------------------------------------------------------------------------
agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(0, 10, 20),
    intensity.mV = c(1, 2, 3)
  ))
det <- ip_peak_detector(
    "test",
    "CO2",
    detect = mark_peaks(2L),
    details = "a test detector"
  )
msgs <- run_detect(agg, det)$msgs
expect_true(any(grepl("detected 2 CO2 peaks", msgs)))
expect_true(any(grepl("with a test detector", msgs)))
