# Extracted from test-peak-detection.R:406

# test -------------------------------------------------------------------------
agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(0, 10, 20, 30),
    intensity.mV = c(1, 2, 3, 4)
  ))
seen <- NULL
det <- ip_peak_detector(
    "test",
    "CO2",
    start.s = 10,
    stop.s = 30,
    detect = function(t) {
      seen <<- t$time.s
      dplyr::mutate(t, peak = 1L)
    }
  )
run_detect(agg, det)
