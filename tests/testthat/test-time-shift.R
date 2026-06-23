# strip ANSI for matching the printed representation
ts_format <- function(x) {
  cli::ansi_strip(paste(cli::cli_fmt(print(x)), collapse = " "))
}

# two-mass traces where mass 45's gaussian is shifted +0.3 s from mass 44
shift_traces <- function() {
  t <- seq(95, 105, by = 0.1)
  g <- function(center) 100 * exp(-0.5 * ((t - center) / 1)^2)
  dplyr::bind_rows(
    tibble::tibble(
      species = "CO2",
      mass = "44",
      tp = seq_along(t),
      time.s = t,
      intensity.mV = g(100.0)
    ),
    tibble::tibble(
      species = "CO2",
      mass = "45",
      tp = seq_along(t),
      time.s = t,
      intensity.mV = g(100.3)
    )
  )
}

# the matching peaks tibble (apex at t = 100.0 -> tp 51, t = 100.3 -> tp 54)
shift_peaks <- function(rectangular = FALSE) {
  tibble::tibble(
    species = "CO2",
    mass = c("44", "45"),
    peak = 1L,
    detection_mass = c(TRUE, FALSE),
    apex.idx = c(51L, 54L),
    rectangular = rectangular
  )
}

test_that("ip_time_shift_detector() builds and validates the object", {
  d <- ip_time_shift_detector(
    "custom",
    function(traces, peaks) peaks,
    details = "(custom)"
  )
  expect_s3_class(d, "ip_time_shift_detector")
  expect_true(is_time_shift_detector(d))
  expect_false(is_time_shift_detector(list()))
  expect_equal(d$type, "custom")
  expect_true(is.function(d$detect))

  expect_error(
    ip_time_shift_detector("", function(traces, peaks) peaks),
    "non-empty"
  )
  expect_error(
    ip_time_shift_detector("x", "not a function"),
    "must be a function"
  )
  expect_error(
    ip_time_shift_detector("x", function(t, p) p, details = 1),
    "single string"
  )
})

test_that("time shift detectors print like background detectors", {
  expect_match(ts_format(ip_no_time_shift_detector()), "no time shift detector")
  expect_match(
    ts_format(ip_parabolic_time_shift_detector(apex_window = 9)),
    "parabolic time shift detector .*9 point.*parabolic apex window"
  )
})

test_that("ip_no_time_shift_detector() sets time_shift.s = 0", {
  res <- ip_no_time_shift_detector()$detect(shift_traces(), shift_peaks())
  expect_equal(res$time_shift.s, c(0, 0))
  # handles an empty peaks tibble
  empty <- ip_no_time_shift_detector()$detect(
    shift_traces(),
    shift_peaks()[0, ]
  )
  expect_equal(nrow(empty), 0L)
  expect_true("time_shift.s" %in% names(empty))
})

test_that("ip_parabolic_time_shift_detector() validates apex_window", {
  expect_error(ip_parabolic_time_shift_detector(2), "odd integer >= 3")
  expect_error(ip_parabolic_time_shift_detector(8), "odd integer >= 3") # even
  expect_no_error(ip_parabolic_time_shift_detector(3))
  expect_no_error(ip_parabolic_time_shift_detector(9))
})

test_that("ip_parabolic_time_shift_detector() measures the per-mass apex offset", {
  res <- ip_parabolic_time_shift_detector(apex_window = 9)$detect(
    shift_traces(),
    shift_peaks()
  )
  # detection mass -> 0, shifted mass -> ~ +0.3 s
  expect_equal(res$time_shift.s[res$mass == "44"], 0)
  expect_equal(res$time_shift.s[res$mass == "45"], 0.3, tolerance = 1e-2)
})

test_that("ip_parabolic_time_shift_detector() forces rectangular peaks to 0", {
  res <- ip_parabolic_time_shift_detector()$detect(
    shift_traces(),
    shift_peaks(rectangular = TRUE)
  )
  expect_equal(res$time_shift.s, c(0, 0))
})

test_that("ip_parabolic_time_shift_detector() requires the needed columns", {
  d <- ip_parabolic_time_shift_detector()
  pk <- shift_peaks()
  expect_error(
    d$detect(shift_traces(), dplyr::select(pk, -"apex.idx")),
    "apex.idx"
  )
  expect_error(
    d$detect(shift_traces(), dplyr::select(pk, -"detection_mass")),
    "detection_mass"
  )
  expect_error(
    d$detect(shift_traces(), dplyr::select(pk, -"rectangular")),
    "rectangular"
  )
})

test_that("ip_detect_peaks() applies the time shift detector and validates it", {
  agg <- structure(list(traces = shift_traces()), class = "ir_aggregated_data")
  det <- ip_peak_detector(
    "test",
    "CO2",
    detect = function(t) {
      t$peak <- 1L
      t$detection_mass <- t$mass == "44"
      t$bgrd.mV <- 0
      t
    }
  )

  # default parabolic detector -> mass 45 shifted relative to detection mass 44
  peaks <- run_detect(agg, det)$result$peaks
  expect_true("time_shift.s" %in% names(peaks))
  expect_equal(peaks$time_shift.s[peaks$mass == "44"], 0)
  expect_equal(peaks$time_shift.s[peaks$mass == "45"], 0.3, tolerance = 1e-2)

  # overriding with the no-shift detector gives 0 everywhere
  peaks0 <- run_detect(
    agg,
    det,
    time_shift = ip_no_time_shift_detector()
  )$result$peaks
  expect_equal(peaks0$time_shift.s, c(0, 0))

  # invalid time_shift argument errors
  expect_error(
    ip_detect_peaks(agg, det, time_shift = "nope"),
    "ip_time_shift_detector"
  )
})

# apply_peak_time_shift() ----

# a mass trace with 0.2 s spacing and points on both sides of the peak window so
# the shift has room to pull points in / interpolate at the edges
ts_full_trace <- function() {
  tibble::tibble(
    species = "CO2",
    mass = "45",
    tp = 8:16,
    time.s = seq(1.6, 3.2, by = 0.2),
    intensity.mV = c(0, 1, 2, 4, 6, 8, 10, 8, 5)
  )
}

test_that("apply_peak_time_shift() shifts the window and interpolates fractional ends", {
  full <- ts_full_trace()
  # the detected peak window = tp 10-13 (time 2.0-2.6)
  peak_traces <- dplyr::filter(full, .data$tp %in% 10:13) |>
    dplyr::mutate(peak = 1L, detection_mass = FALSE)
  # +0.5 s = two whole 0.2 s points + a 0.1 s fractional remainder -> window
  # moves to [2.5, 3.1]
  peaks <- tibble::tibble(
    species = "CO2",
    mass = "45",
    peak = 1L,
    time_shift.s = 0.5
  )
  out <- apply_peak_time_shift(full, peak_traces, peaks) |>
    dplyr::arrange(.data$time.s)

  # all four original points are kept; the two interpolated boundary points
  # (tp = NA) sit at 2.5 and 3.1
  expect_equal(out$time.s, c(2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.1))
  expect_equal(out$tp, c(10L, 11L, 12L, NA, 13L, 14L, 15L, NA))
  # interpolated: 2.5 between (2.4, 6) & (2.6, 8) = 7; 3.1 between (3.0, 8) & (3.2, 5) = 6.5
  expect_equal(out$intensity.mV, c(2, 4, 6, 7, 8, 10, 8, 6.5))
  # the originals before the shifted window are flagged removed; the brought-in
  # real points and both interpolated boundaries are flagged added
  expect_equal(
    out$tf_removed,
    c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  )
  expect_equal(
    out$tf_added,
    c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE)
  )
  # the active window (kept points) is [2.5, 3.1]
  expect_equal(out$time.s[!out$tf_removed], c(2.5, 2.6, 2.8, 3.0, 3.1))
})

test_that("apply_peak_time_shift() adds no interpolated points for a whole-point shift", {
  full <- ts_full_trace()
  peak_traces <- dplyr::filter(full, .data$tp %in% 10:13) |>
    dplyr::mutate(peak = 1L, detection_mass = FALSE)
  # +0.4 s = exactly two 0.2 s points (no fractional remainder) -> window [2.4, 3.0]
  peaks <- tibble::tibble(
    species = "CO2",
    mass = "45",
    peak = 1L,
    time_shift.s = 0.4
  )
  out <- apply_peak_time_shift(full, peak_traces, peaks) |>
    dplyr::arrange(.data$time.s)
  expect_equal(out$time.s, c(2.0, 2.2, 2.4, 2.6, 2.8, 3.0))
  expect_false(anyNA(out$tp)) # no interpolated boundary points
  expect_equal(out$tf_removed, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(out$tf_added, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))
})

test_that("apply_peak_time_shift() leaves zero-shift (detection mass) traces untouched", {
  full <- ts_full_trace()
  peak_traces <- dplyr::filter(full, .data$tp %in% 10:13) |>
    dplyr::mutate(peak = 1L, detection_mass = TRUE)
  peaks <- tibble::tibble(
    species = "CO2",
    mass = "45",
    peak = 1L,
    time_shift.s = 0
  )
  out <- apply_peak_time_shift(full, peak_traces, peaks)
  expect_equal(nrow(out), nrow(peak_traces))
  expect_false(any(out$tf_removed))
  expect_false(any(out$tf_added))
  expect_equal(out$intensity.mV, peak_traces$intensity.mV)
})

test_that("ip_detect_peaks() shifts each mass's window to its own peak (not aligned)", {
  # mass 45 peaks 0.3 s later than the detection mass 44 (see shift_traces())
  agg <- structure(list(traces = shift_traces()), class = "ir_aggregated_data")
  det <- ip_peak_detector(
    "test",
    "CO2",
    detect = function(t) {
      t$peak <- 1L
      t$detection_mass <- t$mass == "44"
      t
    }
  )
  out <- run_detect(agg, det)$result
  pk <- out$peaks
  # the shift is recorded (mass 45 ~ +0.3 s, the detection mass 44 = 0)
  expect_equal(pk$time_shift.s[pk$mass == "44"], 0)
  expect_equal(pk$time_shift.s[pk$mass == "45"], 0.3, tolerance = 1e-2)
  # the windows are deliberately NOT aligned: mass 45 sits ~0.3 s later, at its
  # own true apex (~100.3 s), while the detection mass apex is at 100.0 s
  expect_equal(pk$apex.s[pk$mass == "44"], 100.0, tolerance = 0.05)
  expect_equal(pk$apex.s[pk$mass == "45"], 100.3, tolerance = 0.05)

  # the peak traces keep all points with the time-shift flags; the shift removed
  # some of mass 45's leading points, the detection mass is untouched
  pt <- out$peak_traces
  expect_true(all(c("tf_removed", "tf_added") %in% names(pt)))
  expect_true(any(pt$tf_removed[pt$mass == "45"]))
  expect_false(any(pt$tf_removed[pt$mass == "44"]))
  expect_false(any(pt$tf_added[pt$mass == "44"]))

  # the removed points keep their measured values and still carry a background
  rem45 <- dplyr::filter(pt, .data$mass == "45", .data$tf_removed)
  orig45 <- dplyr::filter(shift_traces(), .data$mass == "45")
  expect_equal(
    rem45$intensity.mV,
    orig45$intensity.mV[match(rem45$tp, orig45$tp)]
  )
  expect_false(anyNA(pt$bgrd.mV)) # every row (incl. removed/added) has a background
})
