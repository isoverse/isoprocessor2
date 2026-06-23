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
