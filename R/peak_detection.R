# peak detector class ========

#' Create a peak detector
#'
#' Create an `ip_peak_detector` object describing where (and how) peaks should be
#' detected in a set of chromatographic traces. An `ip_peak_detector` is a tibble
#' with one row per species/interval combination and the following columns:
#'
#' - `type` (character): the type of the detector. A single call uses one type,
#'   but detectors with **different** types can be row-bound together (via [c()])
#'   into a single object; the type is stored as a column (rather than an
#'   attribute) precisely so that this is possible.
#' - `species` (character): the species the detector applies to, or `NA` to apply
#'   it to **all** species.
#' - `start.s` / `stop.s` (double): the time interval (in seconds) the detector
#'   covers.
#' - `detect` (list of functions): a function that takes a tibble of traces and
#'   returns a tibble of traces (the actual detection step, run later).
#' - `details` (character): an optional plain-text description of the detector,
#'   shown by [print()]. Defaults to `NA`.
#'
#' The time interval can be provided either in seconds (`start.s` / `stop.s`) or
#' in minutes (`start.min` / `stop.min`); seconds take precedence and are what is
#' stored. `species`, and the `start`/`stop` interval bounds, may be vectors:
#' multiple species and multiple intervals are expanded into one row per
#' species/interval combination. Intervals must not overlap within a `species`
#' (across all detector types); touching intervals, where one stops exactly where
#' the next starts, are allowed. An all-species (`NA`) detector is checked against
#' every concrete species, since it applies to all of them.
#'
#' Combine multiple detectors with [c()], which row-binds them, sorts by `type`,
#' `species` and `start.s`, and re-checks that intervals do not overlap within
#' each species.
#'
#' @param type the type of the detector (required, a single non-empty string).
#'   Stored as the first column so detectors with different types can be combined.
#' @param species character vector of species the detector applies to. The
#'   default `NA_character_` applies the detector to all species. A vector
#'   produces one set of intervals per species.
#' @param start.s,stop.s the start/stop of the detection interval(s) in seconds.
#'   Must be the same length. Take precedence over `start.min` / `stop.min`.
#' @param start.min,stop.min the start/stop of the detection interval(s) in
#'   minutes (converted to seconds for `start.s` / `stop.s` if those are not
#'   provided directly).
#' @param detect a function that takes a tibble of traces and returns a tibble of
#'   traces. Defaults to an identity function (no-op).
#' @param details an optional plain-text description of the detector (a single
#'   string), shown by [print()]. Defaults to `NA_character_`.
#' @return an `ip_peak_detector` tibble (one row per species/interval).
#' @examples
#' # a single detector covering all species over the whole trace
#' ip_peak_detector("all peaks")
#'
#' # two non-overlapping intervals (in minutes) for a specific species
#' ip_peak_detector("CO2 peaks", "CO2", start.min = c(0, 5), stop.min = c(5, Inf))
#' @export
ip_peak_detector <- function(
  type,
  species = NA_character_,
  start.s = 60 * start.min,
  stop.s = 60 * stop.min,
  start.min = 0,
  stop.min = Inf,
  detect = function(traces) traces,
  details = NA_character_
) {
  # argument checks
  check_arg(
    type,
    !missing(type) && is_scalar_character(type) && !is.na(type) && nzchar(type),
    "must be a non-empty string"
  )
  check_arg(
    species,
    is.character(species) && length(species) >= 1L,
    "must be a character vector (use NA for all species)"
  )
  check_arg(start.s, is.numeric(start.s), "must be numeric")
  check_arg(stop.s, is.numeric(stop.s), "must be numeric")
  check_arg(detect, is_function(detect), "must be a function")
  check_arg(details, is_scalar_character(details), "must be a single string")

  # start/stop must be the same length
  if (length(start.s) != length(stop.s)) {
    cli_abort(c(
      "{.field start.s}/{.field start.min} and {.field stop.s}/{.field stop.min} must have the same length",
      "i" = "{.field start.s} has length {length(start.s)}, {.field stop.s} has length {length(stop.s)}"
    ))
  }

  # intervals must be valid (no NA, start strictly before stop)
  if (anyNA(start.s) || anyNA(stop.s)) {
    cli_abort("{.field start.s} and {.field stop.s} must not be {.val {NA}}")
  }
  invalid <- which(!(start.s < stop.s))
  if (length(invalid) > 0L) {
    cli_abort(c(
      "each detection interval must have {.field start.s} < {.field stop.s}",
      "x" = "{qty(invalid)}interval{?s} {invalid}: {.field start.s} = {start.s[invalid]}, {.field stop.s} = {stop.s[invalid]}"
    ))
  }

  # assemble one row per species/interval combination
  intervals <- tibble(start.s = as.double(start.s), stop.s = as.double(stop.s))
  detector <- as.character(species) |>
    purrr::map(function(sp) {
      dplyr::mutate(intervals, species = sp, .before = 1L)
    }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      type = type,
      detect = rep(list(detect), dplyr::n()),
      details = details
    ) |>
    dplyr::relocate("type") |>
    dplyr::arrange(.data$type, .data$species, .data$start.s) |>
    new_ip_peak_detector()

  # intervals must not overlap within a species (across all types)
  check_peak_detector_overlaps(detector)

  return(detector)
}

# construct an ip_peak_detector from a (validated) tibble
new_ip_peak_detector <- function(x) {
  class(x) <- unique(c("ip_peak_detector", class(x)))
  return(x)
}

# is x an ip_peak_detector?
is_peak_detector <- function(x) {
  is(x, "ip_peak_detector")
}

# find overlapping intervals within a single (already validated) set of
# intervals; returns a tibble of the overlapping interval pairs (carrying the
# per-interval `all_species` flag), or NULL if none. intervals are overlapping
# if, once sorted by start, one interval stops *after* the next one starts
# (touching, i.e. stop == next start, is allowed).
find_interval_overlaps <- function(start.s, stop.s, all_species = FALSE) {
  ord <- order(start.s)
  start.s <- start.s[ord]
  stop.s <- stop.s[ord]
  all_species <- rep_len(all_species, length(ord))[ord]
  bad <- which(utils::head(stop.s, -1L) > utils::tail(start.s, -1L))
  if (length(bad) == 0L) {
    return(NULL)
  }
  tibble(
    start1 = start.s[bad],
    stop1 = stop.s[bad],
    all1 = all_species[bad],
    start2 = start.s[bad + 1L],
    stop2 = stop.s[bad + 1L],
    all2 = all_species[bad + 1L]
  )
}

# check that detection intervals do not overlap within any species (across all
# detector types); all-species (NA) intervals are checked against every concrete
# species. throws an informative error listing the offending intervals if they do.
check_peak_detector_overlaps <- function(detector, .call = caller_env()) {
  # rows applying to all species (NA) must be checked against every concrete
  # species, since an all-species detector overlaps a specific-species one in time
  all_species <- dplyr::filter(detector, is.na(.data$species))
  concrete <- dplyr::filter(detector, !is.na(.data$species))
  species_list <- unique(concrete$species)

  # groups to check: each concrete species combined with the all-species rows;
  # if there are no concrete species, just check the all-species rows themselves
  if (length(species_list) > 0L) {
    keys <- species_list
    parts <- purrr::map(
      species_list,
      function(sp) {
        dplyr::bind_rows(
          dplyr::filter(concrete, .data$species == sp),
          all_species
        )
      }
    )
  } else {
    keys <- NA_character_
    parts <- list(all_species)
  }

  for (i in seq_along(parts)) {
    part <- parts[[i]]
    overlaps <- find_interval_overlaps(
      part$start.s,
      part$stop.s,
      is.na(part$species)
    )
    if (!is.null(overlaps)) {
      concrete_group <- !is.na(keys[[i]])
      sp <- if (concrete_group) keys[[i]] else "all species"
      # annotate an interval as "(all species)" when it comes from the
      # all-species detector and we are reporting under a concrete species
      annotate <- function(start.s, stop.s, all) {
        label <- format_peak_interval(start.s, stop.s)
        if (concrete_group && all) {
          format_inline("{label} ({.emph all species})")
        } else {
          label
        }
      }
      msgs <- overlaps |>
        purrr::pmap_chr(function(start1, stop1, all1, start2, stop2, all2) {
          format_inline(
            "{col_magenta(sp)}: {annotate(start1, stop1, all1)} overlaps {annotate(start2, stop2, all2)}"
          )
        })
      cli_abort(
        c(
          "peak detector intervals cannot overlap",
          set_names(msgs, rep("x", length(msgs)))
        ),
        call = .call
      )
    }
  }
  return(invisible(detector))
}

# combine peak detectors ========

#' @export
c.ip_peak_detector <- function(...) {
  objs <- list(...)
  if (!all(purrr::map_lgl(objs, is_peak_detector))) {
    cli_abort(
      "all objects to combine must be {.cls ip_peak_detector} (from {.fn ip_peak_detector})"
    )
  }
  combined <- dplyr::bind_rows(objs) |>
    dplyr::arrange(.data$type, .data$species, .data$start.s) |>
    new_ip_peak_detector()

  # intervals must still not overlap within a species after combining
  check_peak_detector_overlaps(combined)

  return(combined)
}

# print peak detectors ========

# format a single detection interval for printing
format_peak_interval <- function(start.s, stop.s) {
  from <- if (start.s <= 0) {
    format_inline("{.field start}")
  } else {
    format_inline("{.field {secs_to_text(start.s)}}")
  }
  to <- if (is.infinite(stop.s)) {
    format_inline("{.field end}")
  } else {
    format_inline("{.field {secs_to_text(stop.s)}}")
  }
  format_inline("{from} {symbol$arrow_right} {to}")
}

#' @export
print.ip_peak_detector <- function(x, ...) {
  if (nrow(x) == 0L) {
    cli_rule(center = "{.strong empty peak detector}")
    return(invisible(x))
  }

  # one cli rule per detector type, one bullet line per species (with all of
  # that species' intervals joined by "; ")
  for (typ in unique(x$type)) {
    rows <- x[x$type == typ, ]
    cli_rule(center = "{.strong {.emph {typ}} peak detector}")
    lines <- unique(rows$species) |>
      purrr::map_chr(function(sp) {
        sp_rows <- if (is.na(sp)) {
          rows[is.na(rows$species), ]
        } else {
          rows[!is.na(rows$species) & rows$species == sp, ]
        }
        species_label <- if (is.na(sp)) {
          format_inline("{.emph all species}")
        } else {
          format_inline("{col_magenta(sp)}")
        }
        intervals <- seq_len(nrow(sp_rows)) |>
          purrr::map_chr(function(i) {
            interval_label <- format_peak_interval(
              sp_rows$start.s[[i]],
              sp_rows$stop.s[[i]]
            )
            details <- sp_rows$details[[i]]
            if (!is.na(details) && nzchar(details)) {
              format_inline("{interval_label} ({details})")
            } else {
              interval_label
            }
          })
        paste0(
          symbol$bullet,
          " ",
          species_label,
          ": ",
          paste(intervals, collapse = "; ")
        )
      })
    lines |> cli_bullets_raw() |> cli()
  }
  return(invisible(x))
}

#' @exportS3Method knitr::knit_print
knit_print.ip_peak_detector <- function(x, ...) {
  print(x, ...)
}
