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
              paste0(interval_label, " with ", details)
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
    # cli_bullets_raw does not escape < / > (they are taken as markup); double
    # them so literal angle brackets (e.g. "max height < 50%") render correctly
    lines <- gsub("<", "<<", lines, fixed = TRUE)
    lines <- gsub(">", ">>", lines, fixed = TRUE)
    lines |> cli_bullets_raw() |> cli()
  }
  return(invisible(x))
}

#' @exportS3Method knitr::knit_print
knit_print.ip_peak_detector <- function(x, ...) {
  print(x, ...)
}

# apply a peak detector ========

#' Detect peaks in aggregated isotope data
#'
#' Apply a peak detector (see [ip_peak_detector()]) to aggregated isotope data
#' (from `isoreader2::ir_aggregate_isofiles()`). For each detector entry the
#' matching traces (for that species, with `start.s <= time.s < stop.s`) are
#' passed to the detector's `detect` function, and the combined detected peaks are
#' stored in the `peak_traces` dataset of the returned object.
#'
#' Detector entries that apply to all species (`species = NA`) are expanded to
#' every species present in the data's `traces`. Traces with a missing `time.s`
#' or intensity are dropped, and traces are sorted by `time.s` before detection.
#'
#' Peaks are found on the detector's detection mass, but the returned
#' `peak_traces` (and `peaks`) cover **all** masses over each detected peak, with
#' a logical `detection_mass` column marking the mass detection ran on.
#'
#' A per-peak summary is also stored in the `peaks` dataset, with one row per
#' `uidx`/`analysis`/`species`/`mass`/`detection_mass`/`peak` (whichever of these
#' are present) and columns `start.idx`/`end.idx` (the `tp` time-point index),
#' `start.s`/`end.s` (time), `apex.idx`/`apex.s` (index/time of the maximum
#' intensity), and `amplitude.<unit>` (the maximum intensity, in the intensity
#' column's unit).
#'
#' @param aggregated_data an `ir_aggregated_data` object (from
#'   `isoreader2::ir_aggregate_isofiles()`) that includes a `traces` dataset.
#' @param detector an [ip_peak_detector] object.
#' @return the `aggregated_data` with `peak_traces` (one row per data point) and
#'   `peaks` (one row per peak) datasets added.
#' @export
ip_detect_peaks <- function(aggregated_data, detector) {
  # input checks
  check_arg(
    aggregated_data,
    !missing(aggregated_data) && is(aggregated_data, "ir_aggregated_data"),
    "must be a set of aggregated isofiles (use isoreader2::ir_aggregate_isofiles())"
  )
  check_arg(
    detector,
    !missing(detector) && is_peak_detector(detector),
    format_inline(
      "must be an {.cls ip_peak_detector}, e.g. {.emph ip_slope_based_peak_detector()}"
    )
  )
  if (
    !"traces" %in% names(aggregated_data) ||
      !is.data.frame(aggregated_data[["traces"]])
  ) {
    cli_abort(
      "the aggregated data does not include a {.field traces} dataset to detect peaks in"
    )
  }

  # prepare traces: require tp + time.s + an intensity column, drop incomplete
  # rows, and sort by ascending time
  traces <- aggregated_data[["traces"]]
  if (!"tp" %in% names(traces)) {
    cli_abort(
      "the {.field traces} dataset must have a {.field tp} (time point) column"
    )
  }
  if (!"time.s" %in% names(traces)) {
    cli_abort("the {.field traces} dataset must have a {.field time.s} column")
  }
  intensity_cols <- grep("^intensity\\.", names(traces), value = TRUE)
  if (length(intensity_cols) == 0L) {
    cli_abort(
      "the {.field traces} dataset must have an {.field intensity.*} column"
    )
  }
  traces <- traces |>
    dplyr::filter(
      !is.na(.data$time.s),
      dplyr::if_all(dplyr::all_of(intensity_cols), ~ !is.na(.x))
    ) |>
    dplyr::arrange(.data$time.s)

  # expand all-species (NA) detector entries to every species in the data
  available_species <- unique(traces$species)
  detector <- seq_len(nrow(detector)) |>
    purrr::map(function(i) {
      row <- detector[i, ]
      if (!is.na(row$species)) {
        return(row)
      }
      if (length(available_species) == 0L) {
        return(NULL)
      }
      row[rep(1L, length(available_species)), ] |>
        dplyr::mutate(species = available_species)
    }) |>
    dplyr::bind_rows()

  # run each detector entry over its species and time window
  results <- list()
  for (i in seq_len(nrow(detector))) {
    row <- detector[i, ]
    species <- row$species
    start.s <- row$start.s
    stop.s <- row$stop.s
    details <- row$details[[1]]

    entry_traces <- traces |>
      dplyr::filter(
        .data$species == !!species,
        .data$time.s >= !!start.s,
        .data$time.s < !!stop.s
      )

    entry_start <- start_info()
    detected <- (row$detect[[1]])(entry_traces)

    # nothing detected for this entry
    if (is.null(detected) || nrow(detected) == 0L) {
      next
    }

    # a real return must carry a peak column
    if (!"peak" %in% names(detected)) {
      cli_abort(
        "the detector's {.field detect} function must return a tibble with a {.field peak} column (for {.emph {species}})"
      )
    }
    n_peaks <- suppressWarnings(max(detected$peak, na.rm = TRUE))
    if (!is.finite(n_peaks)) {
      next
    }

    # report the detection ({qty} pins the pluralization quantity to n_peaks);
    # a detector may attach a "filter_info" note about peaks it filtered out
    filter_info <- attr(detected, "filter_info")
    filter_note <- if (!is.null(filter_info)) {
      paste0(" (", filter_info, ")")
    } else {
      ""
    }
    msg <- format_inline(
      "detected {n_peaks} {col_magenta(species)} {qty(n_peaks)}peak{?s}{filter_note} between {format_peak_interval(start.s, stop.s)}"
    )
    if (!is.na(details) && nzchar(details)) {
      msg <- paste0(msg, " with ", details)
    }
    finish_info(msg, start = entry_start)

    results[[length(results) + 1L]] <- detected
  }

  # store the combined detected peaks and their per-peak summary
  aggregated_data[["peak_traces"]] <- dplyr::bind_rows(results)
  aggregated_data[["peaks"]] <- summarize_peak_traces(
    aggregated_data[["peak_traces"]]
  )
  return(aggregated_data)
}

# summarize a peak_traces tibble (one row per data point) into a peaks tibble
# (one row per mass per peak). grouping is by whichever of
# uidx/analysis/species/mass/detection_mass/peak are present. start/end are
# reported as both the time-point index (`tp`, as *.idx) and the time (`time.s`,
# as *.s). when a background column (`bgrd.<unit>`) is present, the apex and
# amplitude are taken from the background-subtracted height (intensity - bgrd)
# and the background at the apex is reported in a `bgrd.<unit>` column; otherwise
# the raw intensity is used. the amplitude column is named after the intensity
# unit (e.g. `amplitude.mV`). returns an empty tibble if there are no peak traces.
summarize_peak_traces <- function(peak_traces) {
  if (nrow(peak_traces) == 0L) {
    return(tibble())
  }

  # locate the intensity column to derive the amplitude (and its unit)
  intensity_col <- grep("^intensity\\.", names(peak_traces), value = TRUE)[1L]
  if (is.na(intensity_col)) {
    cli_abort(
      "{.field peak_traces} must have an {.field intensity.*} column to summarize peaks"
    )
  }
  amplitude_col <- sub("^intensity", "amplitude", intensity_col)

  # the time-point index column (provided by isoreader2) is required for the
  # *.idx summary columns
  if (!"tp" %in% names(peak_traces)) {
    cli_abort(
      "{.field peak_traces} must have a {.field tp} (time point) column to summarize peaks"
    )
  }

  # background-subtracted height (intensity - bgrd) when a background column is
  # present, otherwise the raw intensity
  bgrd_col <- grep("^bgrd\\.", names(peak_traces), value = TRUE)[1L]
  has_bgrd <- !is.na(bgrd_col)
  peak_traces$.height <- if (has_bgrd) {
    peak_traces[[intensity_col]] - peak_traces[[bgrd_col]]
  } else {
    peak_traces[[intensity_col]]
  }

  group_cols <- intersect(
    c("uidx", "analysis", "species", "mass", "detection_mass", "peak"),
    names(peak_traces)
  )

  # apex/amplitude from the height; the bgrd-at-apex column is added when present
  summary_exprs <- rlang::exprs(
    start.idx = .data$tp[which.min(.data$time.s)],
    end.idx = .data$tp[which.max(.data$time.s)],
    start.s = min(.data$time.s),
    end.s = max(.data$time.s),
    apex.idx = .data$tp[which.max(.data[[".height"]])],
    apex.s = .data$time.s[which.max(.data[[".height"]])]
  )
  summary_exprs[[amplitude_col]] <- rlang::expr(max(.data[[".height"]]))
  if (has_bgrd) {
    summary_exprs[[bgrd_col]] <- rlang::expr(
      .data[[!!bgrd_col]][which.max(.data[[".height"]])]
    )
  }

  peak_traces |>
    dplyr::summarize(.by = dplyr::all_of(group_cols), !!!summary_exprs)
}
