# tests for background detectors (background.R)

# ip_bgrd_detector() ----

test_that("ip_bgrd_detector() builds the expected object", {
  d <- ip_bgrd_detector(
    "custom",
    detect = function(traces, peak_traces) peak_traces,
    details = "p = 5"
  )
  expect_s3_class(d, "ip_bgrd_detector")
  expect_type(d, "list")
  expect_identical(d$type, "custom")
  expect_true(is.function(d$detect))
  expect_identical(d$details, "p = 5")
})

test_that("ip_bgrd_detector() validates its arguments", {
  expect_error(ip_bgrd_detector(), "must be a non-empty string")
  expect_error(ip_bgrd_detector(""), "must be a non-empty string")
  expect_error(ip_bgrd_detector("t"), "must be a function")
  expect_error(
    ip_bgrd_detector(
      "t",
      function(traces, peak_traces) peak_traces,
      details = 5
    ),
    "must be a single string"
  )
})

test_that("print.ip_bgrd_detector() prints the type and details", {
  out_no <- cli::ansi_strip(cli::cli_fmt(print(ip_no_bgrd_detector())))
  expect_match(paste(out_no, collapse = " "), "no background detector")
  # no details -> no "with"
  expect_false(any(grepl("with", out_no)))

  d <- ip_bgrd_detector(
    "custom",
    function(traces, peak_traces) peak_traces,
    details = "p = 5"
  )
  out <- cli::ansi_strip(cli::cli_fmt(print(d)))
  expect_match(paste(out, collapse = " "), "custom detector p = 5")
})

# ip_no_bgrd_detector() ----

test_that("ip_no_bgrd_detector() sets the background to zero in the intensity unit", {
  pt <- tibble::tibble(
    species = "CO2",
    mass = "44",
    time.s = c(0, 10),
    intensity.mV = c(1, 5),
    peak = 1L
  )
  res <- ip_no_bgrd_detector()$detect(traces = NULL, peak_traces = pt)
  expect_true("bgrd.mV" %in% names(res))
  expect_equal(res$bgrd.mV, c(0, 0))

  # the bg column follows the intensity unit
  pt_V <- tibble::tibble(
    species = "CO2",
    time.s = 0,
    intensity.V = 0.1,
    peak = 1L
  )
  expect_true("bgrd.V" %in% names(ip_no_bgrd_detector()$detect(NULL, pt_V)))
})

# ip_individual_bgrd_detector() ----

test_that("ip_individual_bgrd_detector() builds the detector and summarizes its parameters", {
  d <- ip_individual_bgrd_detector(history.s = 5)
  expect_s3_class(d, "ip_bgrd_detector")
  expect_identical(d$type, "individual background")
  expect_match(
    cli::ansi_strip(d$details),
    "min intensity of the last 5 seconds before the peak"
  )
  expect_match(cli::ansi_strip(d$details), "without smoothing")

  # points + smoothing + a different function
  d2 <- ip_individual_bgrd_detector(
    history.pts = 25,
    smooth_coefficients = c(0.25, 0.5, 0.25),
    func = max
  )
  expect_match(
    cli::ansi_strip(d2$details),
    "max intensity of the last 25 points before the peak"
  )
  expect_match(cli::ansi_strip(d2$details), "with 3-point smoothing")
})

test_that("ip_individual_bgrd_detector() validates its arguments", {
  # exactly one of history.s / history.pts
  expect_error(ip_individual_bgrd_detector(), "must supply one")
  expect_error(
    ip_individual_bgrd_detector(history.s = 5, history.pts = 25),
    "only one"
  )
  # positive history
  expect_error(ip_individual_bgrd_detector(history.s = -1), "positive number")
  expect_error(ip_individual_bgrd_detector(history.pts = 0), "positive integer")
  # smooth_coefficients must be odd length
  expect_error(
    ip_individual_bgrd_detector(
      history.s = 5,
      smooth_coefficients = c(0.5, 0.5)
    ),
    "odd length"
  )
  # func must be a function
  expect_error(
    ip_individual_bgrd_detector(history.s = 5, func = 5),
    "must be a function"
  )
})

test_that("ip_individual_bgrd_detector() detect computes a per-peak background", {
  # flat baseline of 5, then a step up to 50 (the "peak")
  traces <- tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:50,
    time.s = 0:49,
    intensity.mV = c(rep(5, 20), rep(50, 30))
  )
  peak_traces <- traces |>
    dplyr::filter(.data$tp >= 21) |>
    dplyr::mutate(peak = 1L, detection_mass = TRUE)

  # 10 points before the start (all baseline 5) -> background = 5
  res_pts <- ip_individual_bgrd_detector(history.pts = 10)$detect(
    traces,
    peak_traces
  )
  expect_true("bgrd.mV" %in% names(res_pts))
  expect_true(all(res_pts$bgrd.mV == 5))

  # 10 seconds before the start -> same baseline
  res_s <- ip_individual_bgrd_detector(history.s = 10)$detect(
    traces,
    peak_traces
  )
  expect_true(all(res_s$bgrd.mV == 5))

  # smoothing a flat baseline leaves it unchanged
  res_sm <- ip_individual_bgrd_detector(
    history.pts = 10,
    smooth_coefficients = c(0.25, 0.5, 0.25)
  )$detect(traces, peak_traces)
  expect_equal(unique(res_sm$bgrd.mV), 5)
})

test_that("ip_individual_bgrd_detector() history window includes the peak's first point", {
  # rising baseline so the start point differs from the point before it
  traces <- tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:6,
    time.s = 0:5,
    intensity.mV = c(10, 20, 30, 40, 50, 60)
  )
  peak_traces <- traces |>
    dplyr::filter(.data$tp >= 4) |> # peak starts at tp 4 (intensity 40)
    dplyr::mutate(peak = 1L, detection_mass = TRUE)

  # 1 point of history -> just the peak's first point (40), not the one before
  expect_equal(
    unique(
      ip_individual_bgrd_detector(history.pts = 1)$detect(
        traces,
        peak_traces
      )$bgrd.mV
    ),
    40
  )
  # 0.5 s of history -> only the start point (time 3) qualifies
  expect_equal(
    unique(
      ip_individual_bgrd_detector(history.s = 0.5)$detect(
        traces,
        peak_traces
      )$bgrd.mV
    ),
    40
  )
})

test_that("ip_individual_bgrd_detector() excludes tf_added points and applies bgrd to all rows", {
  # flat baseline of 5, then a step up to 50 (the "peak")
  traces <- tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:30,
    time.s = 0:29,
    intensity.mV = c(rep(5, 15), rep(50, 15))
  )
  # the peak window = tp 16-30; simulate a time shift: the first original point is
  # removed (kept, flagged), and an interpolated boundary point (tp = NA) is added
  peak_traces <- traces |>
    dplyr::filter(.data$tp >= 16) |>
    dplyr::mutate(
      peak = 1L,
      detection_mass = TRUE,
      tf_removed = .data$tp == 16L,
      tf_added = FALSE
    )
  edge <- peak_traces[1, ]
  edge$tp <- NA_integer_
  edge$time.s <- 15.5
  edge$intensity.mV <- 999 # an arbitrary value that must not affect the background
  edge$tf_removed <- FALSE
  edge$tf_added <- TRUE
  pt <- dplyr::bind_rows(peak_traces, edge)

  res <- ip_individual_bgrd_detector(history.pts = 10)$detect(traces, pt)
  # the background is the baseline (5), unaffected by the tp = NA added point, and
  # is present on every row - including the added (tp = NA) and the removed one
  expect_false(anyNA(res$bgrd.mV))
  expect_true(all(res$bgrd.mV == 5))
})

test_that("ip_individual_bgrd_detector() computes background per mass", {
  # two masses with different baselines
  traces <- dplyr::bind_rows(
    tibble::tibble(
      species = "CO2",
      mass = "44",
      tp = 1:30,
      time.s = 0:29,
      intensity.mV = c(rep(5, 15), rep(50, 15))
    ),
    tibble::tibble(
      species = "CO2",
      mass = "45",
      tp = 1:30,
      time.s = 0:29,
      intensity.mV = c(rep(2, 15), rep(20, 15))
    )
  )
  peak_traces <- traces |>
    dplyr::filter(.data$tp >= 16) |>
    dplyr::mutate(peak = 1L, detection_mass = .data$mass == "44")
  res <- ip_individual_bgrd_detector(history.pts = 10)$detect(
    traces,
    peak_traces
  )
  expect_equal(unique(res$bgrd.mV[res$mass == "44"]), 5)
  expect_equal(unique(res$bgrd.mV[res$mass == "45"]), 2)
})

# ip_isodat_default_background() ----

test_that("ip_isodat_default_background() uses the Isodat settings", {
  d <- ip_isodat_default_background()
  expect_s3_class(d, "ip_bgrd_detector")
  expect_identical(d$type, "individual background")
  expect_match(
    cli::ansi_strip(d$details),
    "min intensity of the last 25 points before the peak"
  )
  expect_match(cli::ansi_strip(d$details), "with 5-point smoothing")
})
