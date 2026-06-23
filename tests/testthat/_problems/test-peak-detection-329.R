# Extracted from test-peak-detection.R:329

# test -------------------------------------------------------------------------
pt <- tibble::tibble(
    species = "CO2",
    mass = "44",
    peak = 1L,
    detection_mass = TRUE,
    tp = 1:5,
    time.s = 0:4,
    intensity.mV = c(2, 7, 12, 7, 2),
    bgrd.mV = 2
  )
pk <- summarize_peak_traces(pt)
expect_true("area.mVs" %in% names(pk))
expect_equal(pk[["area.mVs"]], 20)
expect_true("bgrd.mVs" %in% names(pk))
expect_equal(pk[["bgrd.mVs"]], 8)
