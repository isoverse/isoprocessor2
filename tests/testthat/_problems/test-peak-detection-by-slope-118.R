# Extracted from test-peak-detection-by-slope.R:118

# test -------------------------------------------------------------------------
d <- ip_isodat_default_detector(
    c("CO2", "N2"),
    start.min = c(0, 5),
    stop.min = c(5, Inf),
    peak_end_max_height.pct = 30
  )
expect_equal(nrow(d), 4L)
expect_match(detector_details(d), "end max height = 30%")
