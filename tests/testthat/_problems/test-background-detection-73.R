# Extracted from test-background-detection.R:73

# test -------------------------------------------------------------------------
d <- ip_individual_bgrd_detector(history.s = 5)
expect_s3_class(d, "ip_bgrd_detector")
expect_identical(d$type, "individual background")
expect_match(cli::ansi_strip(d$details), "min of last 5 seconds before the peak")
