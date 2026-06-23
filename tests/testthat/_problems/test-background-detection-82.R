# Extracted from test-background-detection.R:82

# test -------------------------------------------------------------------------
d <- ip_individual_bgrd_detector(history.s = 5)
expect_s3_class(d, "ip_bgrd_detector")
expect_identical(d$type, "individual background")
expect_match(cli::ansi_strip(d$details), "min of last 5 seconds before the peak")
expect_match(cli::ansi_strip(d$details), "without smoothing")
d2 <- ip_individual_bgrd_detector(
    history.pts = 25,
    smooth_coefficients = c(0.25, 0.5, 0.25),
    func = max
  )
expect_match(cli::ansi_strip(d2$details), "max of last 25 points before the peak")
