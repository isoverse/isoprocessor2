# ip_peak_detector() construction ----

test_that("ip_peak_detector() builds the expected structure", {
  d <- ip_peak_detector("peaks")
  expect_s3_class(d, "ip_peak_detector")
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 1L)
  # type is the first column
  expect_identical(
    names(d),
    c("type", "species", "start.s", "stop.s", "detect", "details")
  )
  # column types
  expect_type(d$type, "character")
  expect_type(d$species, "character")
  expect_type(d$start.s, "double")
  expect_type(d$stop.s, "double")
  expect_type(d$detect, "list")
  expect_type(d$details, "character")
  expect_true(is.function(d$detect[[1]]))
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
  expect_identical(ip_peak_detector("d", start.min = 1, stop.min = 2)$start.s, 60)
  expect_identical(ip_peak_detector("d", start.min = 1, stop.min = 2)$stop.s, 120)
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
  expect_error(ip_peak_detector("d", start.s = c(0, 1), stop.s = 10), "same length")
  # start not before stop
  expect_error(ip_peak_detector("d", start.s = 10, stop.s = 5), "start.s < stop.s")
  expect_error(ip_peak_detector("d", start.s = 5, stop.s = 5), "start.s < stop.s")
  # NA bounds
  expect_error(ip_peak_detector("d", start.s = NA_real_, stop.s = 5), "must not be")
  # detect must be a function, details must be a single string
  expect_error(ip_peak_detector("d", detect = "nope"), "must be a function")
  expect_error(ip_peak_detector("d", details = 5), "must be a single string")
  expect_error(ip_peak_detector("d", details = c("a", "b")), "must be a single string")
})

test_that("ip_peak_detector() enforces non-overlapping intervals per species", {
  # overlapping -> error
  expect_error(
    ip_peak_detector("d", start.s = c(0, 100), stop.s = c(120, Inf)),
    "cannot overlap"
  )
  # touching (stop == next start) -> allowed
  expect_no_error(ip_peak_detector("d", start.s = c(0, 100), stop.s = c(100, Inf)))
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
    nrow(c(ip_peak_detector("d", "CO2", 0, 50), ip_peak_detector("d", "N2", 0, 50))),
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
      ip_peak_detector("d", c("CO2", "N2"), start.s = c(0, 100), stop.s = c(100, Inf)),
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

test_that("print() shows the plain-text details when provided", {
  # no details -> no parenthetical
  expect_false(any(grepl("\\(", cli::cli_fmt(print(ip_peak_detector("d"))))) )

  d <- ip_peak_detector("d", details = "threshold 500 mV")
  expect_true(any(grepl("threshold 500 mV", cli::cli_fmt(print(d)))))
})

test_that("print() labels NA species as 'all species' and Inf stop as 'end'", {
  out <- cli::cli_fmt(print(ip_peak_detector("d")))
  expect_true(any(grepl("all species", out)))
  expect_true(any(grepl("end", out)))
})
