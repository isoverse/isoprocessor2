# tests for the slope-based peak detector (peak_detection_by_slope.R)
# shared helpers live in helper-peaks.R

# ip_slope_based_peak_detector() ----

test_that("ip_slope_based_peak_detector() builds a slope-based ip_peak_detector", {
  d <- make_slope_detector()
  expect_s3_class(d, "ip_peak_detector")
  expect_true(all(d$type == "slope-based"))
  expect_equal(d$species, "CO2")
  expect_true(is.function(d$detect[[1]]))
  # details summarize the parameters, including the detection mass and unit
  expect_match(detector_details(d), "detection mass = min")
  expect_match(detector_details(d), "mV/s")
  expect_match(detector_details(d), "slope window = 5")
  expect_match(detector_details(d), "max width = 180 s")
  expect_match(detector_details(d), "peak resolution = 80%")
})

test_that("ip_slope_based_peak_detector() requires the detection parameters", {
  # detection_mass is required
  expect_error(
    ip_slope_based_peak_detector(
      "CO2",
      start_slope.mV_s = 1,
      top_slope.mV_s = 0.05,
      end_slope.mV_s = 1,
      slope_window = 5L,
      slope_window_shift = 1L,
      max_peak_width.s = 180,
      peak_resolution.pct = 80
    ),
    "function, a single number, or a single string"
  )
  # the numeric detection params are required too (no defaults)
  expect_error(
    ip_slope_based_peak_detector(
      "CO2",
      detection_mass = min,
      start_slope.mV_s = 1,
      top_slope.mV_s = 0.05,
      end_slope.mV_s = 1
    ),
    "odd integer >= 3"
  )
})

test_that("ip_slope_based_peak_detector() accepts pA/s units and passes detector args through", {
  d <- make_slope_detector(
    species = c("CO2", "N2"),
    start.min = c(0, 5),
    stop.min = c(5, Inf),
    start_slope.mV_s = NA_real_,
    top_slope.mV_s = NA_real_,
    end_slope.mV_s = NA_real_,
    start_slope.pA_s = 10,
    top_slope.pA_s = 1,
    end_slope.pA_s = 10
  )
  # 2 species x 2 intervals
  expect_equal(nrow(d), 4L)
  expect_match(detector_details(d), "pA/s")
})

test_that("ip_slope_based_peak_detector() enforces one consistent slope unit", {
  # both units for one threshold
  expect_error(make_slope_detector(start_slope.pA_s = 1), "not both")
  # neither unit for one threshold
  expect_error(
    make_slope_detector(start_slope.mV_s = NA_real_),
    "must supply one"
  )
  # mixed units across thresholds
  expect_error(
    make_slope_detector(start_slope.mV_s = NA_real_, start_slope.pA_s = 1),
    "same unit"
  )
})

test_that("ip_slope_based_peak_detector() validates type-specific numeric params", {
  expect_error(make_slope_detector(slope_window = 4), "odd integer >= 3")
  expect_error(
    make_slope_detector(slope_window = 5, slope_window_shift = 3),
    "magnitude"
  )
  expect_error(make_slope_detector(max_peak_width.s = -1), "positive")
  expect_error(
    make_slope_detector(peak_resolution.pct = 150),
    "between 0 and 100"
  )
})

test_that("ip_slope_based_peak_detector() accepts numeric/text detection_mass", {
  expect_match(
    detector_details(make_slope_detector(detection_mass = 44)),
    "detection mass = 44"
  )
  expect_match(
    detector_details(make_slope_detector(detection_mass = "44")),
    "detection mass = 44"
  )
  expect_match(
    detector_details(make_slope_detector(detection_mass = max)),
    "detection mass = max"
  )
  expect_error(
    make_slope_detector(detection_mass = c(1, 2)),
    "function, a single number"
  )
})

test_that("ip_slope_based_peak_detector() detect runs the slope-based worker", {
  d <- make_slope_detector()
  # empty traces -> NULL (nothing to detect)
  expect_null(d$detect[[1]](tibble::tibble()))
})

# ip_isodat_default_detector() ----

test_that("ip_isodat_default_detector() applies the Isodat defaults", {
  d <- ip_isodat_default_detector("CO2")
  expect_s3_class(d, "ip_peak_detector")
  expect_true(all(d$type == "slope-based"))
  expect_match(detector_details(d), "detection mass = min")
  expect_match(detector_details(d), "0.2/0.05/0.4 mV/s")
  expect_match(detector_details(d), "slope window = 5")
  expect_match(detector_details(d), "max width = 180 s")
  expect_match(detector_details(d), "peak resolution = 50%")
})

test_that("ip_isodat_default_detector() forwards species/start/stop and overrides", {
  d <- ip_isodat_default_detector(
    c("CO2", "N2"),
    start.min = c(0, 5),
    stop.min = c(5, Inf),
    peak_resolution.pct = 30
  )
  expect_equal(nrow(d), 4L)
  expect_match(detector_details(d), "peak resolution = 30%")
})

# calculate_rolling_slope() ----

test_that("calculate_rolling_slope() recovers a constant slope (asymmetric window)", {
  x <- 1:10
  y <- 2 * x + 3
  s <- calculate_rolling_slope(x, y, window_size = 3, window_shift = 1)
  expect_length(s, 10L)
  # window_size = 3, shift 1 -> before = 0, after = 2; valid for i in 1..(n-2)
  expect_equal(s[1:8], rep(2, 8))
  expect_true(all(is.na(s[9:10])))

  # non-integer slope, larger window: window_size = 5 -> before = 1, after = 3
  y2 <- -0.5 * x + 1
  s2 <- calculate_rolling_slope(x, y2, window_size = 5, window_shift = 1)
  expect_equal(s2[2:7], rep(-0.5, 6))
  expect_true(all(is.na(s2[c(1, 8, 9, 10)])))
})

test_that("calculate_rolling_slope() matches the original .rolling_slope", {
  set.seed(42)
  y <- cumsum(rnorm(50))
  x <- seq_along(y)
  # window_shift = 1 reproduces the original slope_n behaviour
  for (w in c(3L, 5L, 7L)) {
    expect_equal(
      calculate_rolling_slope(x, y, window_size = w, window_shift = 1),
      .rolling_slope(x, y, w)
    )
  }
})

test_that("calculate_rolling_slope() requires window_shift", {
  expect_error(
    calculate_rolling_slope(1:5, 1:5, window_size = 3),
    "must be a single integer"
  )
})

test_that("calculate_rolling_slope() window_shift repositions the window", {
  x <- 1:10
  y <- 2 * x + 3 # constant slope -> recovered wherever the window fits

  # centered (shift 0): window_size 3 -> before 1, after 1; NA at both ends
  s0 <- calculate_rolling_slope(x, y, window_size = 3, window_shift = 0)
  expect_true(is.na(s0[1]) && is.na(s0[10]))
  expect_equal(s0[2:9], rep(2, 8))

  # forward (shift 1, the default): before 0, after 2; NA at the end
  s1 <- calculate_rolling_slope(x, y, window_size = 3, window_shift = 1)
  expect_equal(s1[1:8], rep(2, 8))
  expect_true(all(is.na(s1[9:10])))

  # backward (shift -1): before 2, after 0; NA at the start
  sm <- calculate_rolling_slope(x, y, window_size = 3, window_shift = -1)
  expect_true(all(is.na(sm[1:2])))
  expect_equal(sm[3:10], rep(2, 8))

  # shift larger than half the window errors
  expect_error(
    calculate_rolling_slope(x, y, window_size = 3, window_shift = 2),
    "too large"
  )
  expect_error(
    calculate_rolling_slope(x, y, window_size = 3, window_shift = 1.5),
    "must be a single integer"
  )
})

test_that("calculate_rolling_slope() requires a single odd integer >= 3", {
  expect_error(
    calculate_rolling_slope(1:5, 1:5, window_size = 1, window_shift = 1),
    "odd integer >= 3"
  )
  expect_error(
    calculate_rolling_slope(1:5, 1:5, window_size = 2, window_shift = 1),
    "odd integer >= 3"
  )
  expect_error(
    calculate_rolling_slope(1:5, 1:5, window_size = 4, window_shift = 1),
    "odd integer >= 3"
  )
  expect_error(
    calculate_rolling_slope(1:5, 1:5, window_size = 3.5, window_shift = 1),
    "odd integer >= 3"
  )
  expect_error(
    calculate_rolling_slope(1:5, 1:4, window_size = 3, window_shift = 1),
    "same length"
  )
  expect_error(
    calculate_rolling_slope("a", 1:3, window_size = 3, window_shift = 1),
    "must be numeric"
  )
})

# slope-based detection ----

test_that("detect_slope_based_peaks() finds well-separated peaks at the right places", {
  res <- run_slope_detect(synth_trace(c(150, 400), c(200, 150)))
  expect_s3_class(res, "tbl_df")
  expect_true("peak" %in% names(res))
  expect_equal(max(res$peak), 2L)
  # the apex of each detected peak matches the gaussian centers
  apex <- res |>
    dplyr::group_by(.data$peak) |>
    dplyr::summarize(
      t = .data$time.s[which.max(.data$intensity.mV)],
      .groups = "drop"
    )
  expect_equal(apex$t, c(150, 400))
})

test_that("detect_slope_based_peaks() returns NULL for a flat trace", {
  flat <- tibble::tibble(
    species = "CO2",
    mass = "44",
    time.s = 0:600,
    intensity.mV = 5
  )
  expect_null(run_slope_detect(flat))
})

test_that("detect_slope_based_peaks() selects the detection mass", {
  traces <- dplyr::bind_rows(
    synth_trace(150, 200, mass = "44"),
    synth_trace(150, 100, mass = "45")
  )
  # the mass flagged with detection_mass = TRUE is the one detection ran on
  detected <- function(res) unique(res$mass[res$detection_mass])
  # function picks the minimum / maximum mass
  expect_equal(detected(run_slope_detect(traces, detection_mass = min)), "44")
  expect_equal(detected(run_slope_detect(traces, detection_mass = max)), "45")
  # numeric / string select by as.character match
  expect_equal(detected(run_slope_detect(traces, detection_mass = 45)), "45")
  expect_equal(detected(run_slope_detect(traces, detection_mass = "45")), "45")
  # a mass not present -> NULL
  expect_null(run_slope_detect(traces, detection_mass = 46))
})

test_that("detect_slope_based_peaks() returns all masses flagged by detection_mass", {
  traces <- dplyr::bind_rows(
    synth_trace(c(150, 400), c(200, 150), mass = "44"),
    synth_trace(c(150, 400), c(80, 60), mass = "45")
  )
  res <- run_slope_detect(traces, detection_mass = min)
  # both masses are returned over the detected peaks
  expect_setequal(unique(res$mass), c("44", "45"))
  # detection_mass marks the detection mass (44) only
  expect_true(all(res$detection_mass[res$mass == "44"]))
  expect_false(any(res$detection_mass[res$mass == "45"]))
  # each mass contributes the same time points per peak
  counts <- res |> dplyr::count(.data$mass, .data$peak)
  expect_equal(
    counts$n[counts$mass == "44"],
    counts$n[counts$mass == "45"]
  )
})

test_that("detect_slope_based_peaks() repeats a boundary point shared by adjacent peaks", {
  # two overlapping peaks produce a shoulder split that shares one boundary point
  res <- run_slope_detect(synth_trace(
    c(150, 220),
    c(200, 200),
    sigma = 15,
    t = 0:400
  ))
  expect_equal(max(res$peak), 2L)
  shared <- res |>
    dplyr::distinct(.data$peak, .data$time.s) |>
    dplyr::count(.data$time.s) |>
    dplyr::filter(.data$n > 1L)
  expect_equal(nrow(shared), 1L)
})

test_that("detect_slope_based_peaks() converts slope thresholds within a signal family", {
  mV <- synth_trace(c(150, 400), c(200, 150)) # intensity.mV
  res_mV <- run_slope_detect(mV)
  expect_equal(max(res_mV$peak), 2L)

  # same signal expressed in V; mV/s thresholds get converted to V/s, so the
  # detection is identical
  in_V <- tibble::tibble(
    species = mV$species,
    mass = mV$mass,
    tp = mV$tp,
    time.s = mV$time.s,
    intensity.V = mV$intensity.mV / 1000
  )
  res_V <- run_slope_detect(in_V)
  expect_identical(res_V$time.s, res_mV$time.s)
  expect_identical(res_V$peak, res_mV$peak)

  # current family: data in nA with pA/s thresholds (converted to nA/s)
  in_nA <- tibble::tibble(
    species = mV$species,
    mass = mV$mass,
    tp = mV$tp,
    time.s = mV$time.s,
    intensity.nA = mV$intensity.mV
  )
  expect_equal(max(run_slope_detect(in_nA, unit = "pA_s")$peak), 2L)
})

test_that("detect_slope_based_peaks() errors on incompatible signal kinds", {
  mV <- synth_trace(150, 200) # intensity.mV (voltage)
  make_trace <- function(col) {
    out <- tibble::tibble(species = "CO2", mass = "44", time.s = mV$time.s)
    out[[col]] <- mV$intensity.mV
    out
  }
  # current thresholds vs voltage intensity
  expect_error(run_slope_detect(mV, unit = "pA_s"), "not compatible")
  # voltage thresholds vs current intensity
  expect_error(
    run_slope_detect(make_trace("intensity.nA"), unit = "mV_s"),
    "not compatible"
  )
  # neither a voltage nor a current (cps)
  expect_error(
    run_slope_detect(make_trace("intensity.cps"), unit = "mV_s"),
    "not compatible"
  )
})

test_that("ip_slope_based_peak_detector() detect plugs into ip_detect_peaks", {
  traces <- dplyr::bind_rows(
    synth_trace(c(150, 400), c(200, 150), mass = "44"),
    synth_trace(c(150, 400), c(100, 75), mass = "45")
  )
  agg <- structure(list(traces = traces), class = "ir_aggregated_data")
  out <- run_detect(agg, ip_isodat_default_detector("CO2"))$result
  expect_equal(max(out$peak_traces$peak), 2L)
  # all masses are returned, with detection_mass flagging the detection mass (44)
  expect_setequal(unique(out$peak_traces$mass), c("44", "45"))
  expect_equal(
    unique(out$peak_traces$mass[out$peak_traces$detection_mass]),
    "44"
  )
  # the peaks summary has one row per mass per peak (2 masses x 2 peaks)
  expect_equal(nrow(out$peaks), 4L)
  expect_true("detection_mass" %in% names(out$peaks))
  # trace columns (e.g. tp) carry through to peak_traces
  expect_true("tp" %in% names(out$peak_traces))
})

# background detector integration ----

test_that("ip_slope_based_peak_detector() requires an ip_bgrd_detector", {
  # missing bgrd_detector
  expect_error(
    ip_slope_based_peak_detector(
      "CO2",
      detection_mass = min,
      start_slope.mV_s = 1,
      top_slope.mV_s = 0.05,
      end_slope.mV_s = 1,
      slope_window = 5L,
      slope_window_shift = 1L,
      max_peak_width.s = 180,
      peak_resolution.pct = 80
    ),
    "ip_bgrd_detector"
  )
  # wrong type
  expect_error(make_slope_detector(bgrd_detector = "nope"), "ip_bgrd_detector")
})

test_that("slope detector worker does not add a background column", {
  # the worker returns the raw peak traces (no bgrd.<unit>): the background is
  # added later, by ip_detect_peaks(), from the detector's bundled bgrd detector
  res <- run_slope_detect(synth_trace(c(150, 400), c(200, 150)))
  expect_false(any(grepl("^bgrd\\.", names(res))))

  # the background detector's own print form still shows up in the detector summary
  d <- make_slope_detector(bgrd_detector = ip_no_bgrd_detector())
  expect_match(detector_details(d), "background = no background detector")
})

test_that("ip_detect_peaks() runs the detector's bundled background detector", {
  agg <- structure(
    list(traces = synth_trace(c(150, 400), c(200, 150))),
    class = "ir_aggregated_data"
  )
  # the no-background detector -> bgrd.mV column of zeros in the peak traces
  out <- run_detect(
    agg,
    ip_isodat_default_detector("CO2", bgrd_detector = ip_no_bgrd_detector())
  )$result
  expect_true("bgrd.mV" %in% names(out$peak_traces))
  expect_true(all(out$peak_traces$bgrd.mV == 0))

  # a background detector that fails to add a bg column is an error
  bad <- ip_isodat_default_detector(
    "CO2",
    bgrd_detector = ip_bgrd_detector("bad", function(traces, peak_traces) {
      peak_traces
    })
  )
  expect_error(ip_detect_peaks(agg, bad), "must return the peak traces with a")
})

# min height filtering ----

test_that("ip_slope_based_peak_detector() validates min_height", {
  # at most one of mV/pA
  expect_error(
    make_slope_detector(min_height.mV = 50, min_height.pA = 50),
    "only one"
  )
  # must be positive
  expect_error(make_slope_detector(min_height.mV = -5), "positive")
  # one is fine
  expect_no_error(make_slope_detector(min_height.mV = 50))
  expect_match(
    detector_details(make_slope_detector(min_height.mV = 50)),
    "min height = 50 mV"
  )
})

test_that("slope detector filters peaks below min_height and renumbers", {
  t <- seq(0, 800)
  # three peaks of 200, 30 and 150 mV (the 30 mV one is below a 50 mV threshold)
  y <- 200 *
    exp(-((t - 150)^2) / (2 * 12^2)) +
    30 * exp(-((t - 400)^2) / (2 * 12^2)) +
    150 * exp(-((t - 650)^2) / (2 * 12^2))
  agg <- structure(
    list(
      traces = tibble::tibble(
        species = "CO2",
        mass = "44",
        tp = seq_along(t),
        time.s = t,
        intensity.mV = y
      )
    ),
    class = "ir_aggregated_data"
  )

  out <- run_detect(agg, ip_isodat_default_detector("CO2", min_height.mV = 50))
  # 2 kept, renumbered to 1, 2 (the middle 30 mV peak removed)
  expect_equal(max(out$result$peaks$peak), 2L)
  expect_setequal(out$result$peaks$peak, c(1L, 2L))
  expect_equal(out$result$peaks$apex.s, c(150, 650))
  # the message reports the filtered peak
  expect_match(
    paste(out$msgs, collapse = " "),
    "1 additional below the min height threshold filtered out"
  )

  # a low threshold keeps all three
  out_all <- run_detect(
    agg,
    ip_isodat_default_detector("CO2", min_height.mV = 1)
  )
  expect_equal(max(out_all$result$peaks$peak), 3L)
})

test_that("slope detector min_height converts prefixes and checks signal class", {
  t <- seq(0, 400)
  y <- 100 * exp(-((t - 200)^2) / (2 * 12^2)) # 100 mV peak (baseline 0)
  agg_V <- structure(
    list(
      traces = tibble::tibble(
        species = "CO2",
        mass = "44",
        tp = seq_along(t),
        time.s = t,
        intensity.V = y / 1000
      )
    ),
    class = "ir_aggregated_data"
  )

  # min_height.mV = 50 -> 0.05 V; the 0.1 V peak is kept
  keep <- run_detect(
    agg_V,
    ip_isodat_default_detector("CO2", min_height.mV = 50)
  )
  expect_equal(max(keep$result$peaks$peak), 1L)

  # min_height.mV = 150 -> 0.15 V; the 0.1 V peak is filtered out
  filt <- run_detect(
    agg_V,
    ip_isodat_default_detector("CO2", min_height.mV = 150)
  )
  expect_equal(nrow(filt$result$peaks), 0L)

  # a current-unit min height against voltage intensity is an error
  bad <- make_slope_detector(min_height.mV = NA_real_, min_height.pA = 50)
  expect_error(ip_detect_peaks(agg_V, bad), "not compatible")
})
