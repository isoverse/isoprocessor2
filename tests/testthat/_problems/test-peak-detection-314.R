# Extracted from test-peak-detection.R:314

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

# test -------------------------------------------------------------------------
d <- make_slope_detector()
expect_error(d$detect[[1]](tibble::tibble()), "not yet implemented")
