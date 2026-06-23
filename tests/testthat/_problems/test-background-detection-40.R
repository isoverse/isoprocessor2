# Extracted from test-background-detection.R:40

# test -------------------------------------------------------------------------
out_no <- cli::ansi_strip(cli::cli_fmt(print(ip_no_bgrd_detector())))
expect_match(paste(out_no, collapse = " "), "no background detector")
expect_false(any(grepl("with", out_no)))
d <- ip_bgrd_detector(
    "custom",
    function(traces, peak_traces) peak_traces,
    details = "p = 5"
  )
out <- cli::ansi_strip(cli::cli_fmt(print(d)))
expect_match(paste(out, collapse = " "), "custom detector with p = 5")
