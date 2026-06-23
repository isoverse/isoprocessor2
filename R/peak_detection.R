# peak detector class ========

#' Peak detectors
#'
#' Peak detectors describe where (and how) peaks should be detected in a set of
#' chromatographic traces; apply them with [ip_detect_peaks()]. Use the generic
#' `ip_peak_detector()` to build a custom detector, the
#' `ip_slope_based_peak_detector()` for slope-based detection, or the
#' `ip_isodat_default_detector()` for the common Isodat preset.
#'
#' `ip_peak_detector()` creates an `ip_peak_detector` object: a tibble with one
#' row per species/interval combination and the following columns:
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
#' are present): `start.idx`/`end.idx` (the `tp` time-point index) and
#' `start.s`/`end.s` (time) of the peak bounds, `apex.idx`/`apex.s` (index/time
#' of the maximum background-subtracted height), `amplitude.<unit>` (that maximum
#' background-subtracted height), `bgrd.<unit>` (the background at the apex),
#' `area.<unit>s` (the integrated background-subtracted peak area) and
#' `bgrd.<unit>s` (the integrated background area).
#'
#' Each peak is also flagged `rectangular`: a peak is considered rectangular when
#' its (background-subtracted) area fills more than `rectangularity_factor` of the
#' bounding box, i.e. `area / ((end.s - start.s) * amplitude) >
#' rectangularity_factor`. Rectangularity is determined from the detection mass
#' and applied to all masses of the peak. Reference (e.g. on/off) peaks tend to be
#' rectangular, so when `flag_rectangular_peaks_as_ref_peaks = TRUE` (the default)
#' an additional `ref_peak` column is added that mirrors `rectangular`.
#'
#' Once the peaks (and their rectangularity) are known, the `time_shift` detector
#' (see [ip_time_shift_detector()]) is applied to add a `time_shift.s` column: the
#' per-mass apex time offset relative to the detection mass. The default
#' [ip_parabolic_time_shift_detector()] measures it from a parabolic fit of each
#' apex; rectangular peaks always get `time_shift.s = 0`.
#'
#' The `area` argument selects how the areas are integrated:
#' - `"trapezoidal"` (the default) applies the proper trapezoidal rule and is
#'   interval-safe: it uses the actual time spacing between points, so it is
#'   correct even when the points are not evenly spaced.
#' - `"isodat"` reproduces the area calculation used by Isodat (and Qtegra): it
#'   assumes all time intervals are identical and effectively counts the first
#'   and last data point in full, whereas the trapezoidal rule counts them at
#'   half weight. This usually matters little because peaks tend towards zero at
#'   their ends, but it differs slightly for peaks that do not return to baseline.
#'   Use it to reproduce Isodat values exactly.
#'
#' @param aggregated_data an `ir_aggregated_data` object (from
#'   `isoreader2::ir_aggregate_isofiles()`) that includes a `traces` dataset.
#' @param detector an [ip_peak_detector] object.
#' @param area the peak-area integration method, `"trapezoidal"` (default) or
#'   `"isodat"`; see Details.
#' @param rectangularity_factor the rectangularity cutoff (default `0.55`): a peak
#'   is flagged `rectangular` when `area / ((end.s - start.s) * amplitude)`
#'   exceeds this value (see Details).
#' @param flag_rectangular_peaks_as_ref_peaks whether to add a `ref_peak` column
#'   to the `peaks` dataset that mirrors the `rectangular` flag (default `TRUE`).
#' @param time_shift an [ip_time_shift_detector] used to add the `time_shift.s`
#'   column to the `peaks` dataset (default [ip_parabolic_time_shift_detector()]).
#' @return the `aggregated_data` with `peak_traces` (one row per data point) and
#'   `peaks` (one row per peak) datasets added.
#' @export
ip_detect_peaks <- function(
  aggregated_data,
  detector,
  area = c("trapezoidal", "isodat"),
  rectangularity_factor = 0.55,
  flag_rectangular_peaks_as_ref_peaks = TRUE,
  time_shift = ip_parabolic_time_shift_detector()
) {
  # input checks
  check_arg(
    aggregated_data,
    !missing(aggregated_data) && is(aggregated_data, "ir_aggregated_data"),
    "must be a set of aggregated isofiles (use isoreader2::ir_aggregate_isofiles())"
  )
  area <- arg_match(area)
  area_fn <- switch(area, trapezoidal = trapezoidal_area, isodat = isodat_area)
  check_arg(
    rectangularity_factor,
    is.numeric(rectangularity_factor) &&
      length(rectangularity_factor) == 1L &&
      !is.na(rectangularity_factor) &&
      rectangularity_factor > 0,
    "must be a single positive number"
  )
  check_arg(
    flag_rectangular_peaks_as_ref_peaks,
    is_scalar_logical(flag_rectangular_peaks_as_ref_peaks) &&
      !is.na(flag_rectangular_peaks_as_ref_peaks),
    "must be TRUE or FALSE"
  )
  check_arg(
    time_shift,
    is_time_shift_detector(time_shift),
    format_inline(
      "must be an {.cls ip_time_shift_detector}, e.g. {.emph ip_parabolic_time_shift_detector()}"
    )
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
  overall_start <- start_info()
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

    # a real return must carry peak and detection_mass columns
    if (!"peak" %in% names(detected)) {
      cli_abort(
        "the detector's {.field detect} function must return a tibble with a {.field peak} column (for {.emph {species}})"
      )
    }
    if (!"detection_mass" %in% names(detected)) {
      cli_abort(
        "the detector's {.field detect} function must return a tibble with a {.field detection_mass} column (for {.emph {species}})"
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
      "detected {.pkg {n_peaks} {col_magenta(species)} {qty(n_peaks)}peak{?s}}{filter_note} between {format_peak_interval(start.s, stop.s)}"
    )
    if (!is.na(details) && nzchar(details)) {
      msg <- paste0(msg, " with ", details)
    }
    finish_info(msg, start = entry_start)

    results[[length(results) + 1L]] <- detected
  }

  # store the combined detected peaks and their per-peak summary
  aggregated_data[["peak_traces"]] <- dplyr::bind_rows(results)
  peaks <- summarize_peak_traces(
    aggregated_data[["peak_traces"]],
    area_fn = area_fn,
    rectangularity_factor = rectangularity_factor,
    flag_ref = flag_rectangular_peaks_as_ref_peaks
  )
  # apply the time shift detector now that the peaks (and their rectangularity)
  # are known, adding a time_shift.s column to the per-peak summary
  if (nrow(peaks) > 0L) {
    peaks <- time_shift$detect(traces, peaks)
  }
  aggregated_data[["peaks"]] <- peaks

  # final summary: how many rectangular (= reference) vs analytical peaks, the
  # area integration method and the time shift detector, counting each peak once
  if (nrow(peaks) > 0L && "rectangular" %in% names(peaks)) {
    peak_keys <- intersect(
      c("uidx", "analysis", "species", "peak"),
      names(peaks)
    )
    peak_flags <- dplyr::distinct(
      peaks,
      dplyr::across(dplyr::all_of(c(peak_keys, "rectangular")))
    )
    n_rect <- sum(peak_flags$rectangular)
    n_analytical <- sum(!peak_flags$rectangular)
  } else {
    n_rect <- 0L
    n_analytical <- 0L
  }
  n_total <- n_rect + n_analytical
  rect_label <- if (flag_rectangular_peaks_as_ref_peaks) {
    "rectangular/ref"
  } else {
    "rectangular"
  }
  finish_info(
    format_inline(
      "detected a total of {.strong {n_total}} {qty(n_total)}peak{?s}: ",
      "{.strong {col_green(n_analytical)}} analytical and ",
      "{.strong {col_yellow(n_rect)}} {rect_label}, ",
      "area integrated with the {.emph {area}} method, ",
      "with {format_time_shift_detector(time_shift)}"
    ),
    start = overall_start
  )
  return(aggregated_data)
}

# summarize a peak_traces tibble (one row per data point) into a peaks tibble
# (one row per mass per peak). grouping is by whichever of
# uidx/analysis/species/mass/detection_mass/peak are present. start/end are
# reported as both the time-point index (`tp`, as *.idx) and the time (`time.s`,
# as *.s). the apex/amplitude/area are taken from the background-subtracted height
# (intensity - bgrd), which requires a background column (`bgrd.<unit>`). reported
# per peak: `amplitude.<unit>` (the max background-subtracted height), `bgrd.<unit>`
# (the background at the apex), `area.<unit>s` (the background-subtracted peak
# area) and `bgrd.<unit>s` (the background area), both integrated with `area_fn`
# (trapezoidal_area or isodat_area). also flags each peak `rectangular` (from its
# detection-mass area/amplitude vs `rectangularity_factor`, applied to all masses)
# and, when `flag_ref`, copies it to `ref_peak`. empty tibble if no peak traces.
summarize_peak_traces <- function(
  peak_traces,
  area_fn = trapezoidal_area,
  rectangularity_factor = 0.55,
  flag_ref = TRUE
) {
  if (nrow(peak_traces) == 0L) {
    return(tibble())
  }

  # locate the intensity column to derive the amplitude/area (and their unit)
  intensity_col <- grep("^intensity\\.", names(peak_traces), value = TRUE)[1L]
  if (is.na(intensity_col)) {
    cli_abort(
      "{.field peak_traces} must have an {.field intensity.*} column to summarize peaks"
    )
  }
  unit <- sub("^intensity\\.", "", intensity_col)
  amplitude_col <- paste0("amplitude.", unit)
  area_col <- paste0("area.", unit, "s") # e.g. intensity.mV -> area.mVs

  # the background column is required (added by the peak detector's background
  # detector); the apex value goes to bgrd.<unit> and its area to bgrd.<unit>s
  bgrd_col <- sub("^intensity\\.", "bgrd.", intensity_col) # e.g. bgrd.mV
  bgrd_area_col <- paste0("bgrd.", unit, "s") # e.g. bgrd.mVs
  if (!bgrd_col %in% names(peak_traces)) {
    cli_abort(c(
      "{.field peak_traces} must have a {.field {bgrd_col}} background column to summarize peaks",
      "i" = "the peak detector did not introduce the required background column (did its {.field bgrd_detector} run?)"
    ))
  }

  # the time-point index column (provided by isoreader2) is required for the
  # *.idx summary columns
  if (!"tp" %in% names(peak_traces)) {
    cli_abort(
      "{.field peak_traces} must have a {.field tp} (time point) column to summarize peaks"
    )
  }

  # the detection mass flag is required (the rectangularity is taken from it)
  if (!"detection_mass" %in% names(peak_traces)) {
    cli_abort(c(
      "{.field peak_traces} must have a {.field detection_mass} column to summarize peaks",
      "i" = "the peak detection function must return a {.field detection_mass} column"
    ))
  }

  # temp columns: background-subtracted height and a copy of the background (the
  # copy avoids the apex bgrd.<unit> summary overwriting the source it is read from)
  peak_traces$.height <- peak_traces[[intensity_col]] - peak_traces[[bgrd_col]]
  peak_traces$.bgrd <- peak_traces[[bgrd_col]]

  group_cols <- intersect(
    c("uidx", "analysis", "species", "mass", "detection_mass", "peak"),
    names(peak_traces)
  )

  summary_exprs <- rlang::exprs(
    start.idx = .data$tp[which.min(.data$time.s)],
    end.idx = .data$tp[which.max(.data$time.s)],
    start.s = min(.data$time.s),
    end.s = max(.data$time.s),
    apex.idx = .data$tp[which.max(.data[[".height"]])],
    apex.s = .data$time.s[which.max(.data[[".height"]])]
  )
  # peak height and background at the apex
  summary_exprs[[amplitude_col]] <- rlang::expr(max(.data[[".height"]]))
  summary_exprs[[bgrd_col]] <- rlang::expr(
    .data[[".bgrd"]][which.max(.data[[".height"]])]
  )
  # integrated areas (method via area_fn): background-subtracted peak area and
  # background area
  summary_exprs[[area_col]] <- rlang::expr(
    (!!area_fn)(.data$time.s, .data[[".height"]])
  )
  summary_exprs[[bgrd_area_col]] <- rlang::expr(
    (!!area_fn)(.data$time.s, .data[[".bgrd"]])
  )

  peaks <- peak_traces |>
    dplyr::summarize(.by = dplyr::all_of(group_cols), !!!summary_exprs)

  # rectangularity flag: a peak is "rectangular" if its (background-subtracted)
  # area fills more than `rectangularity_factor` of the start-to-end x amplitude
  # box, i.e. area / ((end.s - start.s) * amplitude) > rectangularity_factor. it
  # is determined from each peak's detection-mass row and applied to all masses.
  rect_rows <- dplyr::filter(peaks, .data$detection_mass)
  peak_keys <- intersect(c("uidx", "analysis", "species", "peak"), names(peaks))
  rect_flags <- rect_rows |>
    dplyr::mutate(
      rectangular = .data[[area_col]] >
        rectangularity_factor *
          (.data$end.s - .data$start.s) *
          .data[[amplitude_col]]
    ) |>
    dplyr::select(dplyr::all_of(peak_keys), "rectangular")
  peaks <- dplyr::left_join(peaks, rect_flags, by = peak_keys)
  if (flag_ref) {
    # reference peaks are (by default) the rectangular peaks
    peaks$ref_peak <- peaks$rectangular
  }
  peaks
}

# Isodat-style area integral (the difference vs the trapezoidal rule is explained
# in the ip_detect_peaks() documentation): assumes uniform time intervals and
# effectively counts the first/last point in full rather than at half weight.
isodat_area <- function(time.s, intensity) {
  if (length(time.s) < 2L) {
    return(0)
  }
  sum(intensity) * diff(range(time.s)) / length(time.s)
}

# trapezoidal area integral of `intensity` over `time.s`; interval-safe (uses the
# actual time spacing between points) and applies the proper trapezoidal rule
trapezoidal_area <- function(time.s, intensity) {
  if (length(time.s) < 2L) {
    return(0)
  }
  sum(
    diff(time.s) *
      (utils::head(intensity, -1L) + utils::tail(intensity, -1L)) /
      2
  )
}
