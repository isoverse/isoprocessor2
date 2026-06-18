test_that("ip_options() returns all defined options", {
  opts <- ip_options()
  expect_type(opts, "list")
  expect_true(all(c("debug", "auto_use_ansi") %in% names(opts)))
  # defaults
  expect_false(ip_get_option("debug"))
  expect_true(ip_get_option("auto_use_ansi"))
})

test_that("ip_get_options() filters by pattern", {
  expect_named(ip_get_options("ansi"), "auto_use_ansi")
  expect_length(ip_get_options("does_not_exist"), 0L)
})

test_that("ip_get_option() errors for undefined options", {
  expect_error(ip_get_option("not_a_real_option"), "not defined")
})

test_that("ip_options() sets and resets values", {
  old <- ip_options(debug = TRUE)
  on.exit(ip_options(old), add = TRUE)
  expect_true(ip_get_option("debug"))
  # type checking rejects invalid values
  expect_error(ip_options(debug = "yes"), "invalid value")
})
