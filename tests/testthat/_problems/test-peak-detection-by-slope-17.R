# Extracted from test-peak-detection-by-slope.R:17

# test -------------------------------------------------------------------------
d <- make_slope_detector()
expect_s3_class(d, "ip_peak_detector")
expect_true(all(d$type == "slope-based"))
expect_equal(d$species, "CO2")
expect_true(is.function(d$detect[[1]]))
expect_match(detector_details(d), "detection mass = min")
expect_match(detector_details(d), "mV/s")
expect_match(detector_details(d), "slope window = 5")
expect_match(detector_details(d), "max width = 180 s")
expect_match(detector_details(d), "end max height = 20%")
