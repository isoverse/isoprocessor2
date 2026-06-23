# Extracted from test-peak-detection.R:596

# test -------------------------------------------------------------------------
agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(0, 10),
    intensity.mV = c(1, 2)
  ))
det <- ip_peak_detector("test", "CO2", detect = function(t) NULL)
out <- run_detect(agg, det)
expect_equal(nrow(out$result$peak_traces), 0L)
expect_false(any(grepl("detected", out$msgs)))
