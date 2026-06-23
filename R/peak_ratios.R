# peak ratio calculations ==========

# Calculate per-peak intensity (area) ratios for a `peaks` table.
#
# Mirrors isoreader2's ir_calculate_ratios(), but for the peaks summary: instead
# of an intensity column at each time point it uses the integrated peak area
# (`area.<unit>s`) and computes, for every peak, the ratio of each mass's area to
# the area of that peak's base mass. Two columns are added to (and returned with)
# `peaks`:
#   - `ratio_name` (e.g. "45/44"): the mass over the base mass,
#   - `ratio`: area(mass) / area(base mass) within the same
#     uidx/analysis/species/peak group.
# The base mass row of each peak (and any species whose requested base mass is not
# present) gets `NA` in both columns; all rows are kept. Re-running overwrites the
# two columns (idempotent).
#
# The base mass for a species is, by default, the numerically lowest mass measured
# for that species. Override it per species via `...` (named by species, e.g.
# `SO2 = 64, N2 = 28`), exactly like ir_calculate_ratios().
#
# @param peaks a peaks table (from ip_detect_peaks()) with `mass`, `species` and an
#   `area.<unit>s` column.
# @param ... named base masses for individual species (e.g. `CO2 = 44`); species
#   not listed use their numerically lowest measured mass.
# @return `peaks` with the `ratio_name` and `ratio` columns added after the area
#   column.
calculate_peak_ratios <- function(peaks, ...) {
  check_arg(
    peaks,
    !missing(peaks) && is.data.frame(peaks),
    "must be a peaks table (a data frame)"
  )

  # parse base mass overrides from ... (named by species), like ir_calculate_ratios
  base_masses <- rlang::list2(...)
  if (length(base_masses) > 0L) {
    if (!rlang::is_named(base_masses) || any(!nzchar(names(base_masses)))) {
      cli_abort(c(
        "base masses provided in {.arg ...} must be named by species",
        "i" = "e.g. {.code calculate_peak_ratios(peaks, CO2 = 44, N2 = 28)}"
      ))
    }
    ok <- purrr::map_lgl(base_masses, \(x) is.numeric(x) && length(x) == 1L)
    if (any(!ok)) {
      cli_abort(
        "each base mass in {.arg ...} must be a single number (problem with {.field {names(base_masses)[!ok]}})"
      )
    }
    base_masses <- purrr::map_dbl(base_masses, as.numeric)
  }

  # always (re)introduce the ratio columns, even for an empty table
  peaks <- dplyr::select(peaks, -dplyr::any_of(c("ratio_name", "ratio")))
  if (nrow(peaks) == 0L) {
    peaks$ratio_name <- character(0)
    peaks$ratio <- numeric(0)
    return(peaks)
  }

  # the area column supplies the ratio numerator/denominator
  area_col <- grep("^area\\.", names(peaks), value = TRUE)[1L]
  if (is.na(area_col)) {
    cli_abort(
      "the {.field peaks} table must have an {.field area.*} column to calculate ratios"
    )
  }
  if (!"mass" %in% names(peaks)) {
    cli_abort(
      "the {.field peaks} table must have a {.field mass} column to calculate ratios"
    )
  }
  if (!"species" %in% names(peaks)) {
    cli_abort(
      "the {.field peaks} table must have a {.field species} column to calculate ratios"
    )
  }

  # one base area per peak: group by whichever of uidx/analysis/species/peak exist
  group_cols <- intersect(
    c("uidx", "analysis", "species", "peak"),
    names(peaks)
  )

  # the base mass (as it appears in the `mass` column) for each species: an
  # override from `...` if given, else the numerically lowest measured mass
  species_list <- unique(peaks$species)
  base <- purrr::map_chr(species_list, function(sp) {
    sp_masses <- unique(peaks$mass[peaks$species == sp])
    if (sp %in% names(base_masses)) {
      hit <- sp_masses[as.numeric(sp_masses) == base_masses[[sp]]]
      if (length(hit) == 0L) {
        cli_warn(c(
          "the requested base mass {base_masses[[sp]]} for species {.field {sp}} was not found in the {.field peaks} - no ratios calculated for it",
          "i" = "available {qty(length(sp_masses))}mass{?es} for {.field {sp}}: {.field {sort(as.numeric(sp_masses))}}"
        ))
        return(NA_character_)
      }
      hit[1L]
    } else {
      sp_masses[which.min(as.numeric(sp_masses))]
    }
  })
  names(base) <- species_list

  # the base mass area for each peak group
  peaks <- dplyr::mutate(peaks, .base_mass = base[.data$species])
  base_area <- peaks |>
    dplyr::filter(.data$mass == .data[[".base_mass"]]) |>
    dplyr::select(dplyr::all_of(c(group_cols, area_col)))
  names(base_area)[names(base_area) == area_col] <- ".base_area"
  base_area <- dplyr::distinct(base_area)

  # ratio = area(mass) / area(base mass); NA for the base mass row (and where no
  # base mass could be determined). all rows are kept.
  peaks <- peaks |>
    dplyr::left_join(base_area, by = group_cols) |>
    dplyr::mutate(
      .no_ratio = is.na(.data[[".base_mass"]]) |
        .data$mass == .data[[".base_mass"]],
      ratio_name = dplyr::if_else(
        .data$.no_ratio,
        NA_character_,
        sprintf("%s/%s", .data$mass, .data[[".base_mass"]])
      ),
      ratio = dplyr::if_else(
        .data$.no_ratio,
        NA_real_,
        .data[[area_col]] / .data$.base_area
      )
    ) |>
    dplyr::select(-".base_mass", -".base_area", -".no_ratio")

  dplyr::relocate(
    peaks,
    "ratio_name",
    "ratio",
    .after = dplyr::all_of(area_col)
  )
}

# the reference-gas ratio for a group: the mean ratio over the rows flagged as a
# reference peak (ignoring NA ratios, e.g. base mass rows). NA if the group has no
# usable reference peak.
ref_gas_ratio <- function(ratio, ref_peak) {
  vals <- ratio[which(ref_peak)]
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) {
    return(NA_real_)
  }
  mean(vals)
}

# delta calculations ==========

#' Calculate peak delta values versus the reference gas
#'
#' Calculate the per-peak isotope delta values relative to the reference gas for a
#' `peaks` dataset (from [ip_detect_peaks()]). This first calculates the peak area
#' ratios (via the same dynamics as the per-peak ratio calculation: each mass's
#' `area.<unit>s` over its peak's base mass, with the base mass overridable per
#' species through `...`, defaulting to the numerically lowest mass), adding the
#' `ratio_name` and `ratio` columns.
#'
#' It then determines, for each `ratio_name` within each analysis (`uidx` /
#' `analysis`), the reference ratio `ref_ratio` as the mean `ratio` of all peaks
#' flagged `ref_peak` (see [ip_detect_peaks()]'s
#' `flag_rectangular_peaks_as_ref_peaks`). The delta of every peak is then
#' `delta = (ratio / ref_ratio - 1) * 1000`, with a `delta_name` of
#' `paste0("δ", ratio_name)` (e.g. `"δ45/44"`). Base mass rows (and peaks
#' in an analysis without any reference peak) get `NA`.
#'
#' @param aggregated_data an `ir_aggregated_data` with a `peaks` dataset (from
#'   [ip_detect_peaks()]) that includes a `ref_peak` flag.
#' @param ... named base masses for individual species (e.g. `CO2 = 44`), passed
#'   on to the ratio calculation. Species not listed use their numerically lowest
#'   measured mass.
#' @return the `aggregated_data` with `ratio_name`/`ratio`, `ref_ratio` and
#'   `delta_name`/`delta` columns added to the `peaks` dataset.
#' @export
ip_calculate_deltas_vs_ref_gas <- function(aggregated_data, ...) {
  check_arg(
    aggregated_data,
    !missing(aggregated_data) && is(aggregated_data, "ir_aggregated_data"),
    "must be a set of aggregated isofiles (use isoreader2::ir_aggregate_isofiles())"
  )
  if (
    !"peaks" %in% names(aggregated_data) ||
      !is.data.frame(aggregated_data[["peaks"]])
  ) {
    cli_abort(
      "the aggregated data does not include a {.field peaks} dataset (run {.fn ip_detect_peaks} first)"
    )
  }

  overall_start <- start_info()

  # the peak area ratios first (adds ratio_name / ratio)
  peaks <- calculate_peak_ratios(aggregated_data[["peaks"]], ...)

  # the reference peak flag is required to find the reference gas
  if (!"ref_peak" %in% names(peaks)) {
    cli_abort(c(
      "the {.field peaks} dataset must have a {.field ref_peak} flag to calculate deltas vs the reference gas",
      "i" = "run {.fn ip_detect_peaks} with {.code flag_rectangular_peaks_as_ref_peaks = TRUE}"
    ))
  }

  # drop pre-existing delta columns so re-running overwrites cleanly
  peaks <- dplyr::select(
    peaks,
    -dplyr::any_of(c("ref_ratio", "delta_name", "delta"))
  )

  if (nrow(peaks) == 0L) {
    peaks$ref_ratio <- numeric(0)
    peaks$delta_name <- character(0)
    peaks$delta <- numeric(0)
    aggregated_data[["peaks"]] <- peaks
    finish_info(
      format_inline(
        "calculated no {.field δ} values (no peaks)"
      ),
      start = overall_start
    )
    return(aggregated_data)
  }

  if (!any(peaks$ref_peak %in% TRUE)) {
    cli_warn(
      "no peaks are flagged as {.field ref_peak} - all {.field delta} values will be {.val {NA}}"
    )
  }

  # the reference ratio per (uidx/analysis, ratio_name) = mean ratio of the
  # reference peaks; then delta = (ratio / ref_ratio - 1) * 1000
  group_cols <- c(intersect(c("uidx", "analysis"), names(peaks)), "ratio_name")
  peaks <- peaks |>
    dplyr::mutate(
      .by = dplyr::all_of(group_cols),
      ref_ratio = ref_gas_ratio(.data$ratio, .data$ref_peak)
    ) |>
    dplyr::mutate(
      delta_name = dplyr::if_else(
        is.na(.data$ratio_name),
        NA_character_,
        paste0("δ", .data$ratio_name) # δ = the lowercase delta sign
      ),
      delta = (.data$ratio / .data$ref_ratio - 1) * 1000
    ) |>
    dplyr::relocate("ref_ratio", "delta_name", "delta", .after = "ratio")

  aggregated_data[["peaks"]] <- peaks

  n_deltas <- sum(!is.na(peaks$delta))
  finish_info(
    format_inline(
      "calculated {numbers_to_text(n_deltas)} {qty(n_deltas)}{.field δ} value{?s} vs the reference gas ",
      "(the mean of the {.field ref_peak} ratios per analysis), ",
      "added {.field delta_name}/{.field delta} to {.field peaks}"
    ),
    start = overall_start
  )
  return(aggregated_data)
}
