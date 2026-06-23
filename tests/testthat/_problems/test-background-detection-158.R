# Extracted from test-background-detection.R:158

# test -------------------------------------------------------------------------
d <- ip_isodat_default_background()
expect_s3_class(d, "ip_bgrd_detector")
expect_identical(d$type, "individual background")
expect_match(cli::ansi_strip(d$details), "min of last 5 seconds before the peak")
