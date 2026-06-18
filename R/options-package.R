#' Package options
#'
#' These options are best set via [ip_options()] and queried via [ip_get_option()].
#' However, the base functions [options()] and [getOption()] work as well but require
#' an `isoprocessor2.` prefix (the package name and a dot) for the option name. Setting
#' an option to a value of `NULL` means that the default is used. [ip_get_options()]
#' is available as an additional convenience function to retrieve a subset of options
#' with a regular expression pattern.
#'
#' @examples
#' # All default options
#' ip_get_options()
#'
#' @param ... set package options, syntax identical to [options()]
#' @return `ip_options()` and `ip_get_options()` return a named list of option
#'   values; `ip_get_option()` returns the value of the single requested option.
#' @describeIn ip_options set/get option values
#' @export
ip_options <- function(...) {
  pkg_options(pkg = "isoprocessor2", pkg_options = get_pkg_options(), ...)
}

#' @param pattern to retrieve multiple options (as a list) with a shared pattern
#' @describeIn ip_options get a subset of option values that fit a pattern
#' @export
ip_get_options <- function(pattern = NULL) {
  pkg_options <- ip_options()
  if (!is.null(pattern)) {
    pkg_options <- pkg_options[grepl(pattern, names(pkg_options))]
  }
  return(pkg_options)
}

#' @describeIn ip_options retrieve the current value of one option (option must be defined for the package)
#' @param x name of the specific option to retrieve
#' @export
ip_get_option <- function(x) {
  get_pkg_option(
    option = x,
    pkg = "isoprocessor2",
    pkg_options = get_pkg_options()
  )
}

#' @rdname ip_options
#' @format NULL
#' @usage NULL
#' @section Options for the isoprocessor2 package:
get_pkg_options <- function() {
  list(
    #' - `debug`: turn on debug mode
    debug = define_pkg_option(default = FALSE, check_fn = is_scalar_logical),
    #' - `auto_use_ansi`: whether to automatically enable correct rendering of stylized (ansi) output in HTML reports from notebooks that call `library(isoprocessor2)`. Can be turned off by calling `isoprocessor2::ip_options(auto_use_ansi = FALSE)` **before** calling `library(isoprocessor2)`.
    auto_use_ansi = define_pkg_option(
      default = TRUE,
      check_fn = is_scalar_logical
    )
  )
}
