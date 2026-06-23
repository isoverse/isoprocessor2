# Extracted from test-peak-detection-by-slope.R:277

# test -------------------------------------------------------------------------
mV <- synth_trace(c(150, 400), c(200, 150))
res_mV <- run_slope_detect(mV)
expect_equal(max(res_mV$peak), 2L)
in_V <- tibble::tibble(
    species = mV$species,
    mass = mV$mass,
    time.s = mV$time.s,
    intensity.V = mV$intensity.mV / 1000
  )
res_V <- run_slope_detect(in_V)
