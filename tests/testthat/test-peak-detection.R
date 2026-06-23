# ip_peak_detector() construction ----

test_that("ip_peak_detector() builds the expected structure", {
  d <- ip_peak_detector("peaks")
  expect_s3_class(d, "ip_peak_detector")
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 1L)
  # type is the first column; bgrd sits between detect and details
  expect_identical(
    names(d),
    c("type", "species", "start.s", "stop.s", "detect", "bgrd", "details")
  )
  # column types
  expect_type(d$type, "character")
  expect_type(d$species, "character")
  expect_type(d$start.s, "double")
  expect_type(d$stop.s, "double")
  expect_type(d$detect, "list")
  expect_type(d$bgrd, "list")
  expect_type(d$details, "character")
  expect_true(is.function(d$detect[[1]]))
  # the bgrd column holds an ip_bgrd_detector (defaults to no background)
  expect_s3_class(d$bgrd[[1]], "ip_bgrd_detector")
  # values
  expect_identical(d$type, "peaks")
  expect_true(is.na(d$species))
  expect_identical(d$start.s, 0)
  expect_identical(d$stop.s, Inf)
  expect_true(is.na(d$details))
})

test_that("ip_peak_detector() requires a single non-empty type", {
  expect_error(ip_peak_detector(), "must be a")
  expect_error(ip_peak_detector(c("a", "b")), "must be a")
  expect_error(ip_peak_detector(NA_character_), "must be a")
  expect_error(ip_peak_detector(""), "must be a")
})

test_that("ip_peak_detector() converts minutes to seconds, seconds take precedence", {
  expect_identical(
    ip_peak_detector("d", start.min = 1, stop.min = 2)$start.s,
    60
  )
  expect_identical(
    ip_peak_detector("d", start.min = 1, stop.min = 2)$stop.s,
    120
  )
  # seconds win over minutes when both are supplied
  d <- ip_peak_detector(
    "d",
    start.s = 10,
    stop.s = 20,
    start.min = 99,
    stop.min = 99
  )
  expect_identical(d$start.s, 10)
  expect_identical(d$stop.s, 20)
})

test_that("ip_peak_detector() expands vectors of species and intervals", {
  # multiple species, single interval
  expect_equal(nrow(ip_peak_detector("d", c("CO2", "N2"))), 2L)
  # single species, multiple intervals
  expect_equal(nrow(ip_peak_detector("d", "CO2", c(0, 100), c(100, Inf))), 2L)
  # species x intervals cross product
  d <- ip_peak_detector("d", c("CO2", "N2"), c(0, 100), c(100, Inf))
  expect_equal(nrow(d), 4L)
  expect_setequal(d$species, c("CO2", "N2"))
  # type and detect/details are carried onto every row
  expect_true(all(d$type == "d"))
  expect_true(all(purrr::map_lgl(d$detect, is.function)))
})

# ip_peak_detector() validation ----

test_that("ip_peak_detector() rejects invalid intervals", {
  # mismatched start/stop lengths
  expect_error(
    ip_peak_detector("d", start.s = c(0, 1), stop.s = 10),
    "same length"
  )
  # start not before stop
  expect_error(
    ip_peak_detector("d", start.s = 10, stop.s = 5),
    "start.s < stop.s"
  )
  expect_error(
    ip_peak_detector("d", start.s = 5, stop.s = 5),
    "start.s < stop.s"
  )
  # NA bounds
  expect_error(
    ip_peak_detector("d", start.s = NA_real_, stop.s = 5),
    "must not be"
  )
  # detect must be a function, details must be a single string
  expect_error(ip_peak_detector("d", detect = "nope"), "must be a function")
  expect_error(ip_peak_detector("d", details = 5), "must be a single string")
  expect_error(
    ip_peak_detector("d", details = c("a", "b")),
    "must be a single string"
  )
  # bgrd must be an ip_bgrd_detector
  expect_error(ip_peak_detector("d", bgrd = "nope"), "ip_bgrd_detector")
})

test_that("ip_peak_detector() enforces non-overlapping intervals per species", {
  # overlapping -> error
  expect_error(
    ip_peak_detector("d", start.s = c(0, 100), stop.s = c(120, Inf)),
    "cannot overlap"
  )
  # touching (stop == next start) -> allowed
  expect_no_error(ip_peak_detector(
    "d",
    start.s = c(0, 100),
    stop.s = c(100, Inf)
  ))
  # unsorted input is still validated correctly
  expect_error(
    ip_peak_detector("d", start.s = c(100, 0), stop.s = c(Inf, 120)),
    "cannot overlap"
  )
})

# c() combining ----

test_that("c() combines, sorts, and re-checks detectors", {
  a <- ip_peak_detector("d", "CO2", 0, 100)
  b <- ip_peak_detector("d", "CO2", 100, Inf)
  combined <- c(a, b)
  expect_s3_class(combined, "ip_peak_detector")
  expect_equal(nrow(combined), 2L)
  # sorted by type, species, then start.s
  expect_identical(combined$start.s, c(0, 100))

  # different species don't conflict
  expect_equal(
    nrow(c(
      ip_peak_detector("d", "CO2", 0, 50),
      ip_peak_detector("d", "N2", 0, 50)
    )),
    2L
  )
})

test_that("c() forbids overlapping intervals per species across all types", {
  a <- ip_peak_detector("d", "CO2", 0, 120)
  b <- ip_peak_detector("d", "CO2", 100, Inf)
  # same species + overlapping time -> error (same type)
  expect_error(c(a, b), "cannot overlap")
  # overlapping time for the same species across DIFFERENT types -> also error
  expect_error(
    c(
      ip_peak_detector("d1", "CO2", 0, 120),
      ip_peak_detector("d2", "CO2", 100, Inf)
    ),
    "cannot overlap"
  )
  # different species is allowed (even with overlapping time)
  expect_no_error(c(
    ip_peak_detector("d1", "CO2", 0, 120),
    ip_peak_detector("d2", "N2", 100, Inf)
  ))
  # only detectors can be combined
  expect_error(c(a, tibble::tibble(x = 1)), "must be")
})

test_that("c() catches all-species (NA) intervals overlapping a concrete species", {
  # the all-species detector at [80, Inf] overlaps the CO2/N2 intervals in time
  expect_error(
    c(
      ip_peak_detector(
        "d",
        c("CO2", "N2"),
        start.s = c(0, 100),
        stop.s = c(100, Inf)
      ),
      ip_peak_detector("d", start.s = 80)
    ),
    "cannot overlap"
  )
  # all-species + concrete that do NOT overlap in time is allowed
  expect_no_error(c(
    ip_peak_detector("d", "CO2", start.s = 0, stop.s = 80),
    ip_peak_detector("d", start.s = 80)
  ))
  # the offending all-species interval is annotated in the error message
  msg <- tryCatch(
    c(
      ip_peak_detector("d", "CO2", start.s = 0, stop.s = 100),
      ip_peak_detector("d", start.s = 80)
    ),
    error = conditionMessage
  )
  expect_match(msg, "all species")
})

# print() ----

test_that("print() renders one line per species with intervals joined by '; '", {
  d <- ip_peak_detector("mydetector", c("CO2", "N2"), c(0, 100), c(100, Inf))
  out <- cli::cli_fmt(print(d))
  expect_true(any(grepl("mydetector", out)))
  expect_true(any(grepl("CO2", out)))
  expect_true(any(grepl("N2", out)))
  # one bullet line per species (CO2 and N2), each joining its two intervals
  bullet_lines <- out[grepl(cli::symbol$bullet, out, fixed = TRUE)]
  expect_length(bullet_lines, 2L)
  expect_true(all(grepl(";", bullet_lines)))
  # invisibly returns the detector
  expect_identical(withVisible(print(d))$visible, FALSE)
})

test_that("print() groups multiple types into separate cli rules", {
  d <- c(
    ip_peak_detector("alpha", "CO2", 0, 100),
    ip_peak_detector("beta", "N2", 0, 100)
  )
  out <- cli::cli_fmt(print(d))
  expect_true(any(grepl("alpha", out)))
  expect_true(any(grepl("beta", out)))
  # one rule line per type
  expect_equal(sum(grepl("peak detector", out)), 2L)
})

test_that("print() shows the details after 'with' when provided", {
  # no details -> no "with" suffix
  expect_false(any(grepl("with", cli::cli_fmt(print(ip_peak_detector("d"))))))

  d <- ip_peak_detector("d", details = "threshold 500 mV")
  expect_true(any(grepl("with threshold 500 mV", cli::cli_fmt(print(d)))))
})

test_that("print() labels NA species as 'all species' and Inf stop as 'end'", {
  out <- cli::cli_fmt(print(ip_peak_detector("d")))
  expect_true(any(grepl("all species", out)))
  expect_true(any(grepl("end", out)))
})
# ip_detect_peaks() ----

test_that("ip_detect_peaks() applies a detector and stores peak_traces", {
  agg <- fake_aggregated(tibble::tibble(
    species = c("CO2", "CO2", "N2"),
    time.s = c(0, 10, 0),
    intensity.mV = c(1, 5, 3)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  out <- run_detect(agg, det)$result
  expect_s3_class(out, "ir_aggregated_data")
  expect_true("peak_traces" %in% names(out))
  # only the CO2 traces were passed to the detector
  expect_equal(nrow(out$peak_traces), 2L)
  expect_true(all(out$peak_traces$species == "CO2"))
  expect_true("peak" %in% names(out$peak_traces))
})

test_that("ip_detect_peaks() requires a tp column in the traces", {
  # build the object directly so no tp is added
  agg <- structure(
    list(
      traces = tibble::tibble(
        species = "CO2",
        time.s = c(0, 10),
        intensity.mV = c(1, 5)
      )
    ),
    class = "ir_aggregated_data"
  )
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  expect_error(ip_detect_peaks(agg, det), "must have a .*tp")
})

test_that("ip_detect_peaks() preserves trace columns (e.g. tp) in peak_traces", {
  # extra trace columns such as isoreader2's `tp` (time point) must carry through
  # to peak_traces for downstream filtering/processing
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    tp = 1:3,
    time.s = c(0, 10, 20),
    intensity.mV = c(1, 5, 2)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  pt <- run_detect(agg, det)$result$peak_traces
  expect_true("tp" %in% names(pt))
  expect_equal(pt$tp, 1:3)
})

test_that("ip_detect_peaks() summarizes peaks into a $peaks dataset", {
  # debug mode keeps the time-point index (*.idx) columns in the $peaks dataset
  withr::local_options(isoprocessor2.debug = TRUE)
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:5,
    time.s = c(0, 10, 20, 30, 40),
    intensity.mV = c(1, 5, 9, 4, 2)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L)) # all rows -> peak 1
  pk <- run_detect(agg, det)$result$peaks
  expect_equal(nrow(pk), 1L)
  expect_equal(pk$peak, 1L)
  # start/end as time-point index (tp) and as time (s)
  expect_equal(pk$start.idx, 1L)
  expect_equal(pk$end.idx, 5L)
  expect_equal(pk$start.s, 0)
  expect_equal(pk$end.s, 40)
  # apex = the maximum intensity (9, at tp = 3 / time = 20)
  expect_equal(pk$apex.idx, 3L)
  expect_equal(pk$apex.s, 20)
  # amplitude column carries the intensity unit and the max intensity
  expect_true("amplitude.mV" %in% names(pk))
  expect_equal(pk[["amplitude.mV"]], 9)
})

test_that("ip_detect_peaks() keeps *.idx columns only in debug mode", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:5,
    time.s = c(0, 10, 20, 30, 40),
    intensity.mV = c(1, 5, 9, 4, 2)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  # default (debug off): only the times (*.s), no *.idx
  pk <- run_detect(agg, det)$result$peaks
  expect_true(all(c("start.s", "apex.s", "end.s") %in% names(pk)))
  expect_false(any(c("start.idx", "apex.idx", "end.idx") %in% names(pk)))
  # debug on: the *.idx columns are kept
  withr::local_options(isoprocessor2.debug = TRUE)
  pk_dbg <- run_detect(agg, det)$result$peaks
  expect_true(all(c("start.idx", "apex.idx", "end.idx") %in% names(pk_dbg)))
})

test_that("ip_detect_peaks() amplitude column matches the intensity unit", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:3,
    time.s = c(0, 10, 20),
    intensity.V = c(0.1, 0.5, 0.2)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  pk <- run_detect(agg, det)$result$peaks
  expect_true("amplitude.V" %in% names(pk))
  expect_equal(pk[["amplitude.V"]], 0.5)
})

test_that("summarize_peak_traces() uses background-subtracted height for apex/amplitude", {
  # raw intensity peaks at tp 2 (25), but background is high there, so the
  # background-subtracted height peaks at tp 3
  pt <- tibble::tibble(
    species = "CO2",
    mass = "44",
    peak = 1L,
    detection_mass = TRUE,
    tp = 1:3,
    time.s = c(0, 10, 20),
    intensity.mV = c(12, 25, 22),
    bgrd.mV = c(0, 20, 0)
  )
  pk <- summarize_peak_traces(pt)
  # height = intensity - bgrd = c(12, 5, 22) -> apex at tp 3
  expect_equal(pk$apex.idx, 3L)
  expect_equal(pk$apex.s, 20)
  expect_equal(pk[["amplitude.mV"]], 22) # 22 - 0
  # the background at the apex is reported
  expect_true("bgrd.mV" %in% names(pk))
  expect_equal(pk[["bgrd.mV"]], 0)
})

test_that("summarize_peak_traces() reports background-subtracted trapezoidal area", {
  # triangle peak (height 10 above a constant background of 2) over 4 s
  pt <- tibble::tibble(
    species = "CO2",
    mass = "44",
    peak = 1L,
    detection_mass = TRUE,
    tp = 1:5,
    time.s = 0:4,
    intensity.mV = c(2, 7, 12, 7, 2),
    bgrd.mV = 2
  )
  pk <- summarize_peak_traces(pt)
  # height above bg = 0,5,10,5,0; trapezoid area = 20 mV*s
  expect_true("area.mVs" %in% names(pk))
  expect_equal(pk[["area.mVs"]], 20)
  # the background area is also reported: bg = 2 over 4 s = 8 mV*s
  expect_true("bgrd.mVs" %in% names(pk))
  expect_equal(pk[["bgrd.mVs"]], 8)
})

test_that("summarize_peak_traces() requires a background column", {
  pt <- tibble::tibble(
    species = "CO2",
    mass = "44",
    peak = 1L,
    detection_mass = TRUE,
    tp = 1:3,
    time.s = c(0, 10, 20),
    intensity.mV = c(1, 5, 2)
  )
  expect_error(summarize_peak_traces(pt), "background column")
})

test_that("summarize_peak_traces() requires a detection_mass column", {
  pt <- tibble::tibble(
    species = "CO2",
    mass = "44",
    peak = 1L,
    tp = 1:3,
    time.s = c(0, 10, 20),
    intensity.mV = c(1, 5, 2),
    bgrd.mV = 0
  )
  expect_error(summarize_peak_traces(pt), "detection_mass")
})

test_that("ip_detect_peaks() requires the detect function to return detection_mass", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = 0:2,
    intensity.mV = c(1, 2, 1)
  ))
  # detect returns peak + background but no detection_mass
  det <- ip_peak_detector(
    "test",
    "CO2",
    detect = function(t) {
      t$peak <- 1L
      t$bgrd.mV <- 0
      t
    }
  )
  expect_error(ip_detect_peaks(agg, det), "detection_mass")
})

test_that("trapezoidal_area() and isodat_area() integrate as described", {
  # plateau over uniform spacing: the trapezoidal rule weights the endpoints at
  # half, isodat counts them in full
  expect_equal(trapezoidal_area(0:4, c(0, 20, 20, 20, 0)), 60)
  expect_equal(isodat_area(0:4, c(0, 20, 20, 20, 0)), 48)
  # the trapezoidal rule is interval-safe (non-uniform spacing)
  expect_equal(trapezoidal_area(c(0, 1, 3), c(0, 10, 0)), 15)
  # fewer than two points -> 0
  expect_equal(trapezoidal_area(0, 5), 0)
  expect_equal(isodat_area(0, 5), 0)
})

test_that("summarize_peaks_initial() flags rectangular peaks (no background)", {
  # peak 1 = square pulse (area fills the box -> rectangular), peak 2 = triangle
  # (area = half the box -> analytical). the initial summary uses the raw signal,
  # so no background column is needed.
  pt <- tibble::tibble(
    species = "CO2",
    mass = "44",
    peak = rep(1:2, each = 5),
    detection_mass = TRUE,
    tp = 1:10,
    time.s = c(0:4, 0:4),
    intensity.mV = c(rep(10, 5), c(0, 5, 10, 5, 0))
  )
  pk <- summarize_peaks_initial(pt)
  expect_equal(pk$rectangular, c(TRUE, FALSE)) # ratios 1.0 and 0.5 vs 0.55
  # the initial summary does not carry ref_peak (added later by ip_detect_peaks)
  expect_false("ref_peak" %in% names(pk))
  # the raw amplitude/area are reported; area only for the detection mass
  expect_equal(pk[["amplitude.mV"]], c(10, 10))
  expect_equal(pk[["area.mVs"]], c(40, 20))

  # the cutoff is adjustable: at 0.4 the triangle (ratio 0.5) is rectangular too
  expect_equal(
    summarize_peaks_initial(pt, rectangularity_factor = 0.4)$rectangular,
    c(TRUE, TRUE)
  )
})

test_that("summarize_peaks_initial() rectangularity comes from the detection mass", {
  # detection mass (44) is a square pulse, the other mass (45) a triangle; both
  # masses should inherit the detection mass's rectangular flag
  pt <- tibble::tibble(
    species = "CO2",
    mass = c(rep("44", 5), rep("45", 5)),
    peak = 1L,
    detection_mass = c(rep(TRUE, 5), rep(FALSE, 5)),
    tp = 1:10,
    time.s = c(0:4, 0:4),
    intensity.mV = c(rep(10, 5), c(0, 5, 10, 5, 0))
  )
  pk <- summarize_peaks_initial(pt)
  expect_equal(nrow(pk), 2L)
  expect_true(all(pk$rectangular)) # both masses flagged from the square detection mass
  # the raw area is computed only for the detection mass (NA for the other mass)
  expect_equal(pk[["area.mVs"]][pk$detection_mass], 40)
  expect_true(is.na(pk[["area.mVs"]][!pk$detection_mass]))
})

test_that("ip_detect_peaks() adds ref_peak mirroring rectangular when requested", {
  tr <- tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:10,
    time.s = c(0:4, 0:4),
    intensity.mV = c(rep(10, 5), c(0, 5, 10, 5, 0)),
    peak = rep(1:2, each = 5),
    detection_mass = TRUE
  )
  agg <- structure(list(traces = tr), class = "ir_aggregated_data")
  det <- ip_peak_detector("test", "CO2", detect = function(t) t)
  # peak 1 rectangular, peak 2 analytical; ref_peak mirrors rectangular
  pk <- run_detect(agg, det)$result$peaks
  expect_equal(pk$rectangular, c(TRUE, FALSE))
  expect_equal(pk$ref_peak, pk$rectangular)
  # ref_peak is omitted when flagging is disabled
  pk_no_ref <- run_detect(
    agg,
    det,
    flag_rectangular_peaks_as_ref_peaks = FALSE
  )$result$peaks
  expect_false("ref_peak" %in% names(pk_no_ref))
  expect_true("rectangular" %in% names(pk_no_ref))
})

test_that("ip_detect_peaks() reports a rectangular/analytical summary", {
  tr <- tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:10,
    time.s = c(0:4, 0:4),
    intensity.mV = c(rep(10, 5), c(0, 5, 10, 5, 0)),
    peak = rep(1:2, each = 5),
    detection_mass = TRUE,
    bgrd.mV = 0
  )
  agg <- structure(list(traces = tr), class = "ir_aggregated_data")
  det <- ip_peak_detector("test", "CO2", detect = function(t) t)

  out <- run_detect(agg, det)
  msg <- paste(out$msgs, collapse = " ")
  expect_match(msg, "1 analytical and 1 rectangular/ref")
  expect_match(msg, "trapezoidal method")
  expect_match(msg, "time shifts via the parabolic time shift detector")

  # without ref flagging the label drops "/ref"; the area + time shift methods
  # are reported
  out2 <- run_detect(
    agg,
    det,
    flag_rectangular_peaks_as_ref_peaks = FALSE,
    area = "isodat",
    time_shift = ip_no_time_shift_detector()
  )
  msg2 <- paste(out2$msgs, collapse = " ")
  expect_match(msg2, "rectangular,")
  expect_no_match(msg2, "rectangular/ref")
  expect_match(msg2, "isodat method")
  expect_match(msg2, "time shifts via the no time shift detector")
})

test_that("ip_detect_peaks() validates rectangularity arguments", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = 0:2,
    intensity.mV = c(1, 2, 1)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  expect_error(
    ip_detect_peaks(agg, det, rectangularity_factor = -1),
    "positive"
  )
  expect_error(
    ip_detect_peaks(agg, det, flag_rectangular_peaks_as_ref_peaks = NA),
    "TRUE or FALSE"
  )
})

test_that("ip_detect_peaks() area argument selects the integration method", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    mass = "44",
    time.s = 0:4,
    intensity.mV = c(0, 20, 20, 20, 0)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  # default is trapezoidal
  expect_equal(run_detect(agg, det)$result$peaks[["area.mVs"]], 60)
  expect_equal(
    run_detect(agg, det, area = "trapezoidal")$result$peaks[["area.mVs"]],
    60
  )
  expect_equal(
    run_detect(agg, det, area = "isodat")$result$peaks[["area.mVs"]],
    48
  )
  expect_error(ip_detect_peaks(agg, det, area = "simpson"), "must be one of")
})

test_that("ip_detect_peaks() peak summary groups by species/mass/peak (slope detector)", {
  tt <- seq(0, 600)
  yy <- 200 *
    exp(-((tt - 150)^2) / (2 * 12^2)) +
    150 * exp(-((tt - 400)^2) / (2 * 15^2))
  traces <- tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = seq_along(tt),
    time.s = tt,
    intensity.mV = yy
  )
  agg <- structure(list(traces = traces), class = "ir_aggregated_data")
  withr::local_options(isoprocessor2.debug = TRUE) # keep the *.idx columns
  pk <- run_detect(agg, ip_isodat_default_detector("CO2"))$result$peaks
  expect_equal(nrow(pk), 2L)
  expect_setequal(
    names(pk)[1:4],
    c("species", "mass", "detection_mass", "peak")
  )
  # the two peak apexes sit at the gaussian centers
  expect_equal(pk$apex.s, c(150, 400))
  expect_true(all(c("start.idx", "end.idx", "apex.idx") %in% names(pk)))
})

test_that("ip_detect_peaks() yields an empty $peaks when nothing is detected", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    mass = "44",
    tp = 1:3,
    time.s = c(0, 10, 20),
    intensity.mV = c(1, 2, 3)
  ))
  det <- ip_peak_detector("test", "CO2", detect = function(t) NULL)
  out <- run_detect(agg, det)$result
  expect_equal(nrow(out$peaks), 0L)
})

test_that("ip_detect_peaks() expands all-species (NA) detectors to every species", {
  agg <- fake_aggregated(tibble::tibble(
    species = c("CO2", "N2"),
    time.s = c(0, 0),
    intensity.mV = c(1, 2)
  ))
  det <- ip_peak_detector("test", detect = mark_peaks(1L)) # species = NA
  out <- run_detect(agg, det)$result
  expect_setequal(out$peak_traces$species, c("CO2", "N2"))
})

test_that("ip_detect_peaks() drops NA time.s/intensity rows and sorts by time", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(20, NA, 0, 10),
    intensity.mV = c(2, 9, 1, NA)
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  pt <- run_detect(agg, det)$result$peak_traces
  # the NA-time and NA-intensity rows are gone, remaining sorted ascending
  expect_equal(pt$time.s, c(0, 20))
})

test_that("ip_detect_peaks() restricts traces to [start.s, stop.s)", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(0, 10, 20, 30),
    intensity.mV = c(1, 2, 3, 4)
  ))
  # capture the time window the detect function actually receives
  seen <- NULL
  det <- ip_peak_detector(
    "test",
    "CO2",
    start.s = 10,
    stop.s = 30,
    detect = function(t) {
      seen <<- t$time.s
      out <- dplyr::mutate(t, peak = 1L)
      out$detection_mass <- TRUE
      out$bgrd.mV <- 0
      out
    }
  )
  run_detect(agg, det)
  expect_equal(seen, c(10, 20)) # >= 10 and < 30
})

test_that("ip_detect_peaks() treats NULL / empty detect results as no peaks", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(0, 10),
    intensity.mV = c(1, 2)
  ))
  det <- ip_peak_detector("test", "CO2", detect = function(t) NULL)
  out <- run_detect(agg, det)
  expect_equal(nrow(out$result$peak_traces), 0L)
  msg <- paste(out$msgs, collapse = " ")
  # nothing detected -> no per-entry detection line, and the summary reports 0
  expect_false(any(grepl("between", out$msgs)))
  expect_match(msg, "a total of 0 peaks")
})

test_that("ip_detect_peaks() errors if the detect result lacks a peak column", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(0, 10),
    intensity.mV = c(1, 2)
  ))
  det <- ip_peak_detector("test", "CO2", detect = function(t) t) # no peak column
  expect_error(ip_detect_peaks(agg, det), "must return a tibble with a")
})

test_that("ip_detect_peaks() reports detection with species, interval and details", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = c(0, 10, 20),
    intensity.mV = c(1, 2, 3)
  ))
  det <- ip_peak_detector(
    "test",
    "CO2",
    detect = mark_peaks(2L),
    details = "a test detector"
  )
  msgs <- run_detect(agg, det)$msgs
  expect_true(any(grepl("detected 2 CO2 peaks", msgs)))
  expect_true(any(grepl("with a test detector", msgs)))
})

test_that("ip_detect_peaks() validates its inputs", {
  agg <- fake_aggregated(tibble::tibble(
    species = "CO2",
    time.s = 0,
    intensity.mV = 1
  ))
  det <- ip_peak_detector("test", "CO2", detect = mark_peaks(1L))
  expect_error(
    ip_detect_peaks(tibble::tibble(x = 1), det),
    "aggregated isofiles"
  )
  expect_error(ip_detect_peaks(agg, tibble::tibble(x = 1)), "ip_peak_detector")
  # no traces dataset
  no_traces <- structure(
    list(metadata = tibble::tibble(x = 1)),
    class = "ir_aggregated_data"
  )
  expect_error(ip_detect_peaks(no_traces, det), "does not include a")
})
