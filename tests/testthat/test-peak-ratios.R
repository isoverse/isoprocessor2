# tests for per-peak ratio calculations (peak_ratios.R)

# a simple peaks table: 3 masses x 2 peaks of one species
ratio_peaks <- function() {
  tibble::tibble(
    species = "CO2",
    mass = c("44", "45", "46", "44", "45", "46"),
    detection_mass = rep(c(TRUE, FALSE, FALSE), 2),
    peak = rep(1:2, each = 3),
    area.mVs = c(200, 2.1, 0.8, 20, 0.22, 0.09)
  )
}

test_that("calculate_peak_ratios() uses the smallest mass as the default base", {
  res <- calculate_peak_ratios(ratio_peaks())
  expect_true(all(c("ratio_name", "ratio") %in% names(res)))
  # base mass (44) rows are NA in both columns
  expect_true(all(is.na(res$ratio_name[res$mass == "44"])))
  expect_true(all(is.na(res$ratio[res$mass == "44"])))
  # ratio = area(mass) / area(base mass) within each peak
  expect_equal(res$ratio[res$mass == "45" & res$peak == 1], 2.1 / 200)
  expect_equal(res$ratio[res$mass == "46" & res$peak == 1], 0.8 / 200)
  expect_equal(res$ratio[res$mass == "45" & res$peak == 2], 0.22 / 20)
  # ratio_name is "<mass>/<base>"
  expect_equal(res$ratio_name[res$mass == "45" & res$peak == 1], "45/44")
  expect_equal(res$ratio_name[res$mass == "46" & res$peak == 1], "46/44")
})

test_that("calculate_peak_ratios() respects a per-species base mass override", {
  res <- calculate_peak_ratios(ratio_peaks(), CO2 = 45)
  # now 45 is the base (NA), and 44/46 are ratioed against it
  expect_true(all(is.na(res$ratio[res$mass == "45"])))
  expect_equal(res$ratio[res$mass == "44" & res$peak == 1], 200 / 2.1)
  expect_equal(res$ratio[res$mass == "46" & res$peak == 1], 0.8 / 2.1)
  expect_equal(res$ratio_name[res$mass == "44" & res$peak == 1], "44/45")
})

test_that("calculate_peak_ratios() computes ratios per species", {
  peaks <- tibble::tibble(
    species = c("CO2", "CO2", "N2", "N2"),
    mass = c("44", "45", "28", "29"),
    peak = 1L,
    area.mVs = c(100, 1, 200, 2)
  )
  res <- calculate_peak_ratios(peaks)
  # CO2 base = 44, N2 base = 28 (each the species' smallest mass)
  expect_equal(res$ratio_name[res$mass == "45"], "45/44")
  expect_equal(res$ratio_name[res$mass == "29"], "29/28")
  expect_equal(res$ratio[res$mass == "45"], 0.01)
  expect_equal(res$ratio[res$mass == "29"], 0.01)
})

test_that("calculate_peak_ratios() groups by uidx/analysis when present", {
  peaks <- tibble::tibble(
    analysis = c("a", "a", "b", "b"),
    species = "CO2",
    mass = c("44", "45", "44", "45"),
    peak = 1L,
    area.mVs = c(100, 1, 50, 2)
  )
  res <- calculate_peak_ratios(peaks)
  # each analysis uses its own base mass area
  expect_equal(res$ratio[res$analysis == "a" & res$mass == "45"], 1 / 100)
  expect_equal(res$ratio[res$analysis == "b" & res$mass == "45"], 2 / 50)
})

test_that("calculate_peak_ratios() warns and skips an unfound base mass", {
  peaks <- tibble::tibble(
    species = "CO2",
    mass = c("44", "45"),
    peak = 1L,
    area.mVs = c(100, 1)
  )
  expect_warning(
    res <- calculate_peak_ratios(peaks, CO2 = 99),
    "was not found"
  )
  expect_true(all(is.na(res$ratio)))
  expect_true(all(is.na(res$ratio_name)))
})

test_that("calculate_peak_ratios() is idempotent (overwrites existing columns)", {
  res1 <- calculate_peak_ratios(ratio_peaks())
  res2 <- calculate_peak_ratios(res1)
  expect_equal(res2$ratio, res1$ratio)
  # the ratio columns are not duplicated
  expect_equal(sum(names(res2) == "ratio"), 1L)
  expect_equal(sum(names(res2) == "ratio_name"), 1L)
})

test_that("calculate_peak_ratios() places the columns after the area column", {
  res <- calculate_peak_ratios(ratio_peaks())
  pos_area <- which(names(res) == "area.mVs")
  expect_equal(names(res)[pos_area + 1:2], c("ratio_name", "ratio"))
})

test_that("calculate_peak_ratios() handles an empty peaks table", {
  empty <- ratio_peaks()[0, ]
  res <- calculate_peak_ratios(empty)
  expect_equal(nrow(res), 0L)
  expect_true(all(c("ratio_name", "ratio") %in% names(res)))
  expect_type(res$ratio_name, "character")
  expect_type(res$ratio, "double")
})

# ip_calculate_deltas_vs_ref_gas() ----

# a peaks table with a ref_peak flag (peaks 1-2 = reference, 3-4 = sample) wrapped
# in a minimal ir_aggregated_data
delta_agg <- function(peaks) {
  structure(list(peaks = peaks), class = "ir_aggregated_data")
}

ref_peaks <- function() {
  tibble::tibble(
    species = "CO2",
    mass = rep(c("44", "45"), 4),
    detection_mass = rep(c(TRUE, FALSE), 4),
    peak = rep(1:4, each = 2),
    ref_peak = rep(c(TRUE, TRUE, FALSE, FALSE), each = 2),
    area.mVs = c(100, 1.0, 100, 1.1, 100, 1.2, 100, 0.9)
  )
}

test_that("ip_calculate_deltas_vs_ref_gas() computes deltas vs the mean ref peak ratio", {
  out <- suppressMessages(ip_calculate_deltas_vs_ref_gas(delta_agg(ref_peaks())))
  expect_s3_class(out, "ir_aggregated_data")
  pk <- out$peaks
  # the ratio columns are also added (calculate_peak_ratios runs first)
  expect_true(all(
    c("ratio_name", "ratio", "ref_ratio", "delta_name", "delta") %in% names(pk)
  ))

  m45 <- pk[pk$mass == "45", ]
  # ref_ratio = mean ratio of the reference peaks (1 & 2): mean(0.01, 0.011)
  expect_equal(unique(m45$ref_ratio), 0.0105)
  # delta = (ratio / ref_ratio - 1) * 1000
  expect_equal(m45$delta, (m45$ratio / 0.0105 - 1) * 1000)
  expect_equal(m45$delta[m45$peak == 3], (0.012 / 0.0105 - 1) * 1000)
  expect_equal(m45$delta[m45$peak == 4], (0.009 / 0.0105 - 1) * 1000)
  # delta_name = "δ" + ratio_name
  expect_equal(unique(m45$delta_name), "δ45/44")

  # base mass rows have NA throughout
  m44 <- pk[pk$mass == "44", ]
  expect_true(all(is.na(m44$ratio)))
  expect_true(all(is.na(m44$ref_ratio)))
  expect_true(all(is.na(m44$delta)))
  expect_true(all(is.na(m44$delta_name)))
})

test_that("ip_calculate_deltas_vs_ref_gas() uses a per-analysis reference ratio", {
  peaks <- tibble::tibble(
    analysis = rep(c("a", "b"), each = 4),
    species = "CO2",
    mass = rep(c("44", "45"), 4),
    detection_mass = rep(c(TRUE, FALSE), 4),
    peak = rep(1:2, each = 2, times = 2),
    ref_peak = rep(c(TRUE, FALSE), each = 2, times = 2),
    area.mVs = c(100, 1, 100, 1.5, 100, 2, 100, 2.5)
  )
  pk <- suppressMessages(ip_calculate_deltas_vs_ref_gas(delta_agg(peaks)))$peaks
  # analysis a: ref_ratio = 0.01 -> sample peak 2 delta = (0.015/0.01 - 1)*1000 = 500
  # analysis b: ref_ratio = 0.02 -> sample peak 2 delta = (0.025/0.02 - 1)*1000 = 250
  expect_equal(
    pk$ref_ratio[pk$analysis == "a" & pk$mass == "45"],
    c(0.01, 0.01)
  )
  expect_equal(
    pk$ref_ratio[pk$analysis == "b" & pk$mass == "45"],
    c(0.02, 0.02)
  )
  expect_equal(
    pk$delta[pk$analysis == "a" & pk$peak == 2 & pk$mass == "45"],
    500
  )
  expect_equal(
    pk$delta[pk$analysis == "b" & pk$peak == 2 & pk$mass == "45"],
    250
  )
})

test_that("ip_calculate_deltas_vs_ref_gas() forwards base mass overrides", {
  pk <- suppressMessages(
    ip_calculate_deltas_vs_ref_gas(delta_agg(ref_peaks()), CO2 = 45)
  )$peaks
  # base mass 45 -> ratios/deltas now for mass 44 (the 45 rows are the base)
  expect_equal(unique(pk$delta_name[pk$mass == "44"]), "δ44/45")
  expect_true(all(is.na(pk$delta[pk$mass == "45"])))
})

test_that("ip_calculate_deltas_vs_ref_gas() warns when there are no reference peaks", {
  peaks <- dplyr::mutate(ref_peaks(), ref_peak = FALSE)
  expect_warning(
    out <- ip_calculate_deltas_vs_ref_gas(delta_agg(peaks)),
    "no peaks are flagged as"
  )
  expect_true(all(is.na(out$peaks$delta)))
})

test_that("ip_calculate_deltas_vs_ref_gas() places the columns and is idempotent", {
  out1 <- suppressMessages(ip_calculate_deltas_vs_ref_gas(delta_agg(ref_peaks())))
  # columns sit after ratio
  nm <- names(out1$peaks)
  expect_equal(
    nm[which(nm == "ratio") + 1:3],
    c("ref_ratio", "delta_name", "delta")
  )
  # re-running overwrites (no duplicate columns, same values)
  out2 <- suppressMessages(ip_calculate_deltas_vs_ref_gas(out1))
  expect_equal(sum(names(out2$peaks) == "delta"), 1L)
  expect_equal(out2$peaks$delta, out1$peaks$delta)
})

test_that("ip_calculate_deltas_vs_ref_gas() validates its inputs", {
  expect_error(
    ip_calculate_deltas_vs_ref_gas(tibble::tibble(x = 1)),
    "aggregated isofiles"
  )
  # no peaks dataset
  no_peaks <- structure(
    list(traces = tibble::tibble(x = 1)),
    class = "ir_aggregated_data"
  )
  expect_error(ip_calculate_deltas_vs_ref_gas(no_peaks), "does not include a")
  # peaks without a ref_peak flag
  no_ref <- delta_agg(dplyr::select(ref_peaks(), -"ref_peak"))
  expect_error(ip_calculate_deltas_vs_ref_gas(no_ref), "ref_peak")
})

test_that("calculate_peak_ratios() validates its inputs", {
  # ... base masses must be named by species
  expect_error(calculate_peak_ratios(ratio_peaks(), 44), "named by species")
  # ... base masses must be single numbers
  expect_error(
    calculate_peak_ratios(ratio_peaks(), CO2 = c(44, 45)),
    "single number"
  )
  # required columns
  expect_error(
    calculate_peak_ratios(dplyr::select(ratio_peaks(), -"area.mVs")),
    "area"
  )
  expect_error(
    calculate_peak_ratios(dplyr::select(ratio_peaks(), -"mass")),
    "mass"
  )
  expect_error(
    calculate_peak_ratios(dplyr::select(ratio_peaks(), -"species")),
    "species"
  )
})
