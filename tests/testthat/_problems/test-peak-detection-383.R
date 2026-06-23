# Extracted from test-peak-detection.R:383

# test -------------------------------------------------------------------------
agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(20, NA, 0, 10),
    intensity.mV = c(2, 9, 1, NA)
  ))
det <- ip_peak_detector("test", "CO2", detect = function(t) dplyr::mutate(t, peak = 1L))
pt <- run_detect(agg, det)$result$peak_traces
