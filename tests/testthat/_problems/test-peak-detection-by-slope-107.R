# Extracted from test-peak-detection-by-slope.R:107

# test -------------------------------------------------------------------------
d <- ip_isodat_default_detector("CO2")
expect_s3_class(d, "ip_peak_detector")
expect_true(all(d$type == "slope-based"))
expect_match(detector_details(d), "detection mass = min")
expect_match(detector_details(d), "0.2/0.05/0.4 mV/s")
expect_match(detector_details(d), "slope window = 5")
expect_match(detector_details(d), "max width = 180 s")
expect_match(detector_details(d), "end max height = 50%")
