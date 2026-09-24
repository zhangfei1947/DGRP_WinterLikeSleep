#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(behavr)
  library(damr)
  library(data.table)
  library(sleepr)
})

data.table::setDTthreads(2)
args <- commandArgs(trailingOnly = TRUE)
input_root <- if (length(args) >= 1) args[[1]] else "../data/raw_dam"
output_root <- if (length(args) >= 2) args[[2]] else file.path(getwd(), "minute_rebuild")
batch_regex <- if (length(args) >= 3) args[[3]] else "^Batch[[:alnum:]]+_DGRP$"
force <- if (length(args) >= 4) tolower(args[[4]]) %in% c("true", "t", "1", "yes") else FALSE

pipeline_version <- "2026-08-26-publication-v3"
analysis_days <- 6L
lookback_minutes <- 10L
coverage_threshold <- 0.95
transition_min_denominator <- 30L
sleep_min_denominator <- 30L

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
batch_output_root <- file.path(output_root, "by_batch")
dir.create(batch_output_root, recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  paste0("pipeline_version=", pipeline_version),
  "phenotype_time_origin=lights_on_ZT0",
  "isolation_occurred_approximately_17h_before_ZT0=true",
  "pre_ZT0_acclimation_excluded_from_phenotype=true",
  paste0("sleep_annotation_lookback_minutes=", lookback_minutes),
  "sleep_definition=all bins in an inactivity run of at least 5 consecutive one-minute bins",
  "light_dark_hours=8:16",
  "temperature_C=21.5",
  "all_monitors_same_chamber=true",
  "batches_run_at_different_calendar_times=true",
  "death_rule=sleepr::curate_dead_animals defaults, followed by lifespan >= 6 days",
  paste0("window_coverage_threshold=", coverage_threshold),
  "DGRP809_GWAS_status=source_typo_corrected_to_DGRP805"
), file.path(output_root, "analysis_config.txt"))

batch_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
batch_dirs <- batch_dirs[grepl(batch_regex, basename(batch_dirs))]
batch_dirs <- batch_dirs[grepl("^Batch[[:alnum:]]+_DGRP$", basename(batch_dirs))]
if (!length(batch_dirs)) stop("No matching Batch*_DGRP directories")

pseudo_datetime_shift <- function(x, seconds) {
  y <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (anyNA(y)) stop("Could not parse start_datetime")
  format(y + seconds, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
}

read_removed <- function(path) {
  if (!file.exists(path)) return(character())
  x <- readLines(path, warn = FALSE)
  if (length(x) <= 1L) return(character())
  trimws(x[-1L])
}

safe_divide <- function(numerator, denominator) {
  fifelse(is.finite(denominator) & denominator > 0, numerator / denominator, NA_real_)
}

logit_corrected <- function(numerator, denominator, min_denominator = transition_min_denominator) {
  fifelse(
    is.finite(denominator) & denominator >= min_denominator,
    log((numerator + 0.5) / (denominator - numerator + 0.5)),
    NA_real_
  )
}

window_definitions <- data.table(
  window = c(
    "ZT00_04", "ZT04_08", "ZT08_12", "ZT12_16", "ZT16_20", "ZT20_24",
    "ZT00_08", "ZT08_16", "ZT16_24", "ZT08_24", "ZT00_24"
  ),
  start_minute = c(0, 240, 480, 720, 960, 1200, 0, 480, 960, 480, 0),
  end_minute = c(240, 480, 720, 960, 1200, 1440, 480, 960, 1440, 1440, 1440)
)
window_definitions[, expected_bins := end_minute - start_minute]

summarize_window <- function(dt, retained_ids, definition) {
  window_name <- definition$window[[1]]
  window_start <- definition$start_minute[[1]]
  window_end <- definition$end_minute[[1]]
  expected_bins <- definition$expected_bins[[1]]

  z <- dt[
    id %chin% retained_ids & zt_minute >= window_start & zt_minute < window_end,
    .(id, day, minute, activity, moving, asleep)
  ]
  setorder(z, id, day, minute)

  z[, `:=`(
    previous_minute = shift(minute),
    previous_moving = shift(moving),
    previous_asleep = shift(asleep)
  ), by = .(id, day)]
  z[, consecutive := !is.na(previous_minute) & minute - previous_minute == 1L]
  z[, sleep_run_id := cumsum(
    is.na(previous_minute) |
      minute - previous_minute != 1L |
      asleep != shift(asleep)
  ), by = .(id, day)]

  base <- z[, .(
    observed_bins = uniqueN(minute),
    duplicate_bins = .N - uniqueN(minute),
    sleep_minutes_observed = sum(asleep, na.rm = TRUE),
    sleep_fraction_observed = mean(asleep, na.rm = TRUE),
    total_activity = sum(activity, na.rm = TRUE),
    awake_minutes = sum(!asleep, na.rm = TRUE),
    waking_activity = safe_divide(sum(activity[!asleep], na.rm = TRUE), sum(!asleep, na.rm = TRUE)),
    wake_transitions = sum(consecutive & !previous_moving & moving, na.rm = TRUE),
    inactive_opportunities = sum(consecutive & !previous_moving, na.rm = TRUE),
    doze_transitions = sum(consecutive & previous_moving & !moving, na.rm = TRUE),
    active_opportunities = sum(consecutive & previous_moving, na.rm = TRUE),
    sleep_to_wake_transitions = sum(consecutive & previous_asleep & !asleep, na.rm = TRUE),
    asleep_opportunities = sum(consecutive & previous_asleep, na.rm = TRUE)
  ), by = .(id, day)]

  bouts <- z[asleep == TRUE, .(bout_minutes = .N), by = .(id, day, sleep_run_id)]
  bout_summary <- bouts[, .(
    longest_sleep_bout_minutes = max(bout_minutes),
    n_sleep_bouts = .N
  ), by = .(id, day)]
  base <- merge(base, bout_summary, by = c("id", "day"), all.x = TRUE)

  grid <- CJ(id = retained_ids, day = seq_len(analysis_days), unique = TRUE)
  out <- merge(grid, base, by = c("id", "day"), all.x = TRUE)
  out[, `:=`(
    window = window_name,
    expected_bins = expected_bins,
    observed_bins = fcoalesce(observed_bins, 0L),
    duplicate_bins = fcoalesce(duplicate_bins, 0L),
    coverage = fcoalesce(observed_bins, 0L) / expected_bins
  )]
  out[, metric_valid := coverage >= coverage_threshold & duplicate_bins == 0L]
  out[, sleep_minutes_scaled := fifelse(metric_valid, sleep_fraction_observed * expected_bins, NA_real_)]
  out[, p_wake := fifelse(metric_valid, safe_divide(wake_transitions, inactive_opportunities), NA_real_)]
  out[, p_doze := fifelse(metric_valid, safe_divide(doze_transitions, active_opportunities), NA_real_)]
  out[, p_sleep_to_wake := fifelse(metric_valid, safe_divide(sleep_to_wake_transitions, asleep_opportunities), NA_real_)]
  out[, p_wake_logit := fifelse(metric_valid, logit_corrected(wake_transitions, inactive_opportunities), NA_real_)]
  out[, p_doze_logit := fifelse(metric_valid, logit_corrected(doze_transitions, active_opportunities), NA_real_)]
  out[, p_sleep_to_wake_logit := fifelse(metric_valid, logit_corrected(sleep_to_wake_transitions, asleep_opportunities), NA_real_)]
  out[, longest_bout_log1p := fifelse(metric_valid, log1p(longest_sleep_bout_minutes), NA_real_)]
  out[, waking_activity_log1p := fifelse(
    metric_valid & awake_minutes >= transition_min_denominator,
    log1p(waking_activity),
    NA_real_
  )]
  out[, bouts_per_sleep_hour := fifelse(
    metric_valid & sleep_minutes_observed >= sleep_min_denominator,
    n_sleep_bouts / (sleep_minutes_observed / 60),
    NA_real_
  )]
  out[, fragmentation_log1p := log1p(bouts_per_sleep_hour)]
  setcolorder(out, c("id", "day", "window", "expected_bins", "observed_bins", "coverage", "metric_valid"))
  out
}

make_contrasts <- function(window_metrics, window_name, value_field, prefix) {
  z <- window_metrics[window == window_name, .(id, day, value = get(value_field))]
  wide <- dcast(z, id ~ day, value.var = "value")
  for (day in seq_len(analysis_days)) {
    day_name <- as.character(day)
    if (!day_name %in% names(wide)) wide[, (day_name) := NA_real_]
  }
  setnames(wide, as.character(seq_len(analysis_days)), paste0("D", seq_len(analysis_days)))
  day_matrix <- as.matrix(wide[, paste0("D", seq_len(analysis_days)), with = FALSE])
  finite_rows <- function(cols) apply(day_matrix[, cols, drop = FALSE], 1, function(x) all(is.finite(x)))
  d6d1 <- rep(NA_real_, nrow(wide))
  valid <- finite_rows(c(1, 6))
  d6d1[valid] <- day_matrix[valid, 6] - day_matrix[valid, 1]

  late2 <- rep(NA_real_, nrow(wide))
  valid <- finite_rows(c(1, 2, 5, 6))
  late2[valid] <- rowMeans(day_matrix[valid, 5:6, drop = FALSE]) - rowMeans(day_matrix[valid, 1:2, drop = FALSE])

  late3 <- rep(NA_real_, nrow(wide))
  valid <- finite_rows(1:6)
  late3[valid] <- rowMeans(day_matrix[valid, 4:6, drop = FALSE]) - rowMeans(day_matrix[valid, 1:3, drop = FALSE])

  slope5 <- rep(NA_real_, nrow(wide))
  valid <- finite_rows(1:6)
  ols_weights <- (1:6 - mean(1:6)) / sum((1:6 - mean(1:6))^2) * 5
  slope5[valid] <- as.numeric(day_matrix[valid, , drop = FALSE] %*% ols_weights)

  out <- data.table(id = wide$id)
  out[, (paste0(prefix, "_D6_D1")) := d6d1]
  out[, (paste0(prefix, "_Late2_Early2")) := late2]
  out[, (paste0(prefix, "_Late3_Early3")) := late3]
  out[, (paste0(prefix, "_Slope5")) := slope5]
  out
}

process_batch <- function(batch_dir) {
  batch <- basename(batch_dir)
  message(format(Sys.time()), " processing ", batch)
  out_dir <- file.path(batch_output_root, batch)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  version_path <- file.path(out_dir, "pipeline_version.txt")
  required_outputs <- file.path(out_dir, c(
    "fly_qc.csv", "daily_qc.csv", "window_metrics.csv", "individual_phenotypes.csv",
    "existing_summary_comparison.csv"
  ))
  if (!force && file.exists(version_path) && identical(readLines(version_path, warn = FALSE)[1], pipeline_version) && all(file.exists(required_outputs))) {
    message("  using validated cache for ", batch)
    return(invisible(TRUE))
  }

  loading_path <- file.path(batch_dir, paste0("loadinginfo_", batch, ".csv"))
  summary_path <- file.path(batch_dir, paste0("summary_", batch, ".csv"))
  removed1_path <- file.path(batch_dir, paste0("removed_list1_", batch, ".csv"))
  removed2_path <- file.path(batch_dir, paste0("removed_list2_", batch, ".csv"))
  if (!file.exists(loading_path) || !file.exists(summary_path)) stop("Missing canonical input for ", batch)

  info <- fread(loading_path, na.strings = character(), showProgress = FALSE)
  info[, Light_original := as.character(Light)]
  info[, Light := "08:16"]
  info[, monitor := as.character(monitor)]

  linked_original <- damr::link_dam_metadata(info, result_dir = batch_dir)
  canonical_meta <- as.data.table(linked_original)[, .(
    id = as.character(id), monitor = as.character(monitor), region_id = as.integer(region_id),
    Batch = as.character(Batch), Genotype = as.character(Genotype), Sex = as.character(Sex),
    Temperature = as.character(Temperature), status = as.character(status),
    start_datetime = as.character(start_datetime), stop_datetime = as.character(stop_datetime)
  )]
  if (canonical_meta[, anyDuplicated(paste(monitor, region_id))]) stop("Non-unique monitor-region mapping in ", batch)

  info_lookback <- copy(info)
  info_lookback[, start_datetime := pseudo_datetime_shift(start_datetime, -lookback_minutes * 60)]
  linked_lookback <- damr::link_dam_metadata(info_lookback, result_dir = batch_dir)
  loaded <- damr::load_dam(linked_lookback)
  loaded_meta <- as.data.table(behavr::meta(loaded))[, .(
    loaded_id = as.character(id), monitor = as.character(monitor), region_id = as.integer(region_id)
  )]
  id_map <- merge(loaded_meta, canonical_meta[, .(id, monitor, region_id)], by = c("monitor", "region_id"), all.x = TRUE)
  if (id_map[, anyNA(id)] || id_map[, anyDuplicated(loaded_id)]) stop("Could not map lookback IDs to canonical IDs in ", batch)
  id_lookup <- setNames(id_map$id, id_map$loaded_id)

  dt <- as.data.table(loaded)[, .(
    loaded_id = as.character(id), t = as.numeric(t), activity = as.integer(activity)
  )]
  dt[, id := unname(id_lookup[loaded_id])]
  dt[, analysis_t := t - lookback_minutes * 60]
  dt[, minute_error_seconds := abs(analysis_t - round(analysis_t / 60) * 60)]
  dt[, minute := as.integer(round(analysis_t / 60))]
  dt <- dt[, .(
    activity = sum(activity, na.rm = TRUE),
    minute_error_seconds = max(minute_error_seconds, na.rm = TRUE)
  ), by = .(loaded_id, id, minute)]
  setorder(dt, id, minute)
  dt[, moving := activity > 0L]
  dt[, previous_minute_for_sleep := shift(minute), by = id]
  dt[, sleep_state_run := cumsum(
    is.na(previous_minute_for_sleep) |
      minute - previous_minute_for_sleep != 1L |
      moving != shift(moving)
  ), by = id]
  dt[, asleep := !moving & .N >= 5L, by = .(id, sleep_state_run)]
  dt[, analysis_t := minute * 60]
  dt <- dt[analysis_t >= 0 & analysis_t < analysis_days * 86400]
  dt[, `:=`(day = minute %/% 1440L + 1L, zt_minute = minute %% 1440L)]
  setorder(dt, id, minute)

  duplicate_by_id <- dt[, .N, by = .(id, minute)][N > 1L, .(duplicate_rows = sum(N - 1L)), by = id]
  if (nrow(duplicate_by_id)) dt <- unique(dt, by = c("id", "minute"))

  summary_existing <- fread(summary_path, showProgress = FALSE)
  summary_existing[, id := as.character(id)]
  retained_ids <- intersect(canonical_meta$id, summary_existing$id)
  removed1 <- read_removed(removed1_path)
  removed2 <- read_removed(removed2_path)

  daily_observed <- dt[, .(
    observed_bins = uniqueN(minute),
    duplicate_bins = .N - uniqueN(minute),
    max_minute_error_seconds = max(minute_error_seconds, na.rm = TRUE),
    light_bins = uniqueN(minute[zt_minute < 480L]),
    dark_bins = uniqueN(minute[zt_minute >= 480L]),
    first_minute = min(minute),
    last_minute = max(minute)
  ), by = .(id, day)]
  daily_grid <- CJ(id = canonical_meta$id, day = seq_len(analysis_days), unique = TRUE)
  daily_qc <- merge(daily_grid, daily_observed, by = c("id", "day"), all.x = TRUE)
  daily_qc[, `:=`(
    observed_bins = fcoalesce(observed_bins, 0L),
    duplicate_bins = fcoalesce(duplicate_bins, 0L),
    light_bins = fcoalesce(light_bins, 0L),
    dark_bins = fcoalesce(dark_bins, 0L)
  )]
  daily_qc[, `:=`(
    total_coverage = observed_bins / 1440,
    light_coverage = light_bins / 480,
    dark_coverage = dark_bins / 960
  )]
  daily_qc <- merge(daily_qc, canonical_meta[, .(id, Batch, Genotype, monitor, region_id)], by = "id", all.x = TRUE)

  fly_coverage <- daily_qc[, .(
    min_total_coverage = min(total_coverage),
    min_light_coverage = min(light_coverage),
    min_dark_coverage = min(dark_coverage),
    D1_light_coverage = light_coverage[day == 1L],
    D1_dark_coverage = dark_coverage[day == 1L],
    D6_light_coverage = light_coverage[day == 6L],
    D6_dark_coverage = dark_coverage[day == 6L],
    days_with_any_missing_bins = sum(total_coverage < 1),
    max_minute_error_seconds = suppressWarnings(max(max_minute_error_seconds, na.rm = TRUE))
  ), by = id]
  fly_qc <- merge(canonical_meta, fly_coverage, by = "id", all.x = TRUE)
  fly_qc <- merge(fly_qc, duplicate_by_id, by = "id", all.x = TRUE)
  fly_qc[, duplicate_rows := fcoalesce(duplicate_rows, 0L)]
  fly_qc[, `:=`(
    retained_by_existing_gigem = id %chin% retained_ids,
    removed_by_curate_dead_animals = id %chin% removed1,
    removed_by_lifespan_rule = id %chin% removed2
  )]
  fly_qc[, exclusion_reason := fcase(
    removed_by_curate_dead_animals, "curate_dead_animals_24h_moving_le_1pct",
    removed_by_lifespan_rule, "lifespan_lt_6_days_after_curation",
    !retained_by_existing_gigem, "not_in_existing_summary_other",
    default = "retained"
  )]
  fly_qc[, primary_coverage_eligible :=
           retained_by_existing_gigem & D1_dark_coverage >= coverage_threshold &
           D6_dark_coverage >= coverage_threshold & duplicate_rows == 0L]
  fwrite(fly_qc, file.path(out_dir, "fly_qc.csv"))
  fwrite(daily_qc, file.path(out_dir, "daily_qc.csv"))

  window_metrics <- rbindlist(lapply(seq_len(nrow(window_definitions)), function(index) {
    summarize_window(dt, retained_ids, window_definitions[index])
  }), use.names = TRUE, fill = TRUE)
  window_metrics <- merge(
    window_metrics,
    canonical_meta[, .(id, Batch, Genotype, Sex, Temperature, monitor, region_id)],
    by = "id", all.x = TRUE
  )
  setcolorder(window_metrics, c("id", "Batch", "Genotype", "monitor", "region_id", "day", "window"))
  fwrite(window_metrics, file.path(out_dir, "window_metrics.csv"))

  phenotype_specs <- list(
    list("ZT08_24", "sleep_minutes_scaled", "NightSleep"),
    list("ZT00_08", "sleep_minutes_scaled", "DaySleep"),
    list("ZT08_12", "sleep_minutes_scaled", "EarlyNightSleep_ZT08_12"),
    list("ZT12_16", "sleep_minutes_scaled", "MidNightSleep_ZT12_16"),
    list("ZT16_20", "sleep_minutes_scaled", "LateMidNightSleep_ZT16_20"),
    list("ZT20_24", "sleep_minutes_scaled", "LateNightSleep_ZT20_24"),
    list("ZT08_16", "sleep_minutes_scaled", "FirstHalfNightSleep_ZT08_16"),
    list("ZT16_24", "sleep_minutes_scaled", "SecondHalfNightSleep_ZT16_24"),
    list("ZT08_24", "longest_bout_log1p", "NightLongestBout_log1p"),
    list("ZT08_24", "p_wake_logit", "NightPWake_logit"),
    list("ZT08_24", "p_doze_logit", "NightPDoze_logit"),
    list("ZT08_24", "waking_activity_log1p", "NightWakingActivity_log1p"),
    list("ZT08_24", "fragmentation_log1p", "NightFragmentation_log1p")
  )
  contrast_tables <- lapply(phenotype_specs, function(spec) {
    make_contrasts(window_metrics, spec[[1]], spec[[2]], spec[[3]])
  })
  individual_phenotypes <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), contrast_tables)
  individual_phenotypes <- merge(
    canonical_meta[, .(id, Batch, Genotype, Sex, Temperature, monitor, region_id)],
    individual_phenotypes,
    by = "id", all.y = TRUE
  )
  fwrite(individual_phenotypes, file.path(out_dir, "individual_phenotypes.csv"))

  night_rebuilt <- dcast(
    window_metrics[window == "ZT08_24", .(id, day, rebuilt_night_sleep = sleep_minutes_scaled)],
    id ~ day,
    value.var = "rebuilt_night_sleep"
  )
  setnames(night_rebuilt, as.character(seq_len(analysis_days)), paste0("rebuilt_D", seq_len(analysis_days)))
  comparison_columns <- c("id", paste0("Day", seq_len(analysis_days), "_Sleep_Time_D"))
  comparison <- merge(summary_existing[, ..comparison_columns], night_rebuilt, by = "id", all.x = TRUE)
  for (day in seq_len(analysis_days)) {
    comparison[, (paste0("difference_D", day)) := get(paste0("rebuilt_D", day)) - get(paste0("Day", day, "_Sleep_Time_D"))]
  }
  comparison[, Batch := batch]
  fwrite(comparison, file.path(out_dir, "existing_summary_comparison.csv"))

  writeLines(pipeline_version, version_path)
  rm(loaded, dt, window_metrics)
  invisible(gc())
  message(format(Sys.time()), " completed ", batch)
  invisible(TRUE)
}

for (batch_dir in batch_dirs) process_batch(batch_dir)

collect_file <- function(filename) {
  paths <- file.path(batch_output_root, basename(batch_dirs), filename)
  paths <- paths[file.exists(paths)]
  rbindlist(lapply(paths, fread, showProgress = FALSE), use.names = TRUE, fill = TRUE)
}

all_fly_qc <- collect_file("fly_qc.csv")
all_daily_qc <- collect_file("daily_qc.csv")
all_window_metrics <- collect_file("window_metrics.csv")
all_individual_phenotypes <- collect_file("individual_phenotypes.csv")
all_existing_comparison <- collect_file("existing_summary_comparison.csv")

fwrite(all_fly_qc, file.path(output_root, "fly_qc_all_batches.csv"))
fwrite(all_daily_qc, file.path(output_root, "daily_qc_all_batches.csv"))
fwrite(all_window_metrics, file.path(output_root, "window_metrics_all_batches.csv"))
fwrite(all_individual_phenotypes, file.path(output_root, "individual_phenotypes_all_batches.csv"))
fwrite(all_existing_comparison, file.path(output_root, "existing_summary_comparison_all_batches.csv"))

phenotype_id_columns <- c("id", "Batch", "Genotype", "Sex", "Temperature", "monitor", "region_id")
phenotype_columns <- setdiff(names(all_individual_phenotypes), phenotype_id_columns)
all_individual_phenotypes[, (phenotype_columns) := lapply(.SD, as.numeric), .SDcols = phenotype_columns]
individual_long <- melt(
  all_individual_phenotypes,
  id.vars = phenotype_id_columns,
  measure.vars = phenotype_columns,
  variable.name = "phenotype",
  value.name = "value",
  variable.factor = FALSE
)
individual_long <- individual_long[is.finite(value)]

line_batch_long <- individual_long[, .(
  value = median(value),
  n_flies = .N,
  mad = mad(value, constant = 1),
  q25 = quantile(value, 0.25),
  q75 = quantile(value, 0.75)
), by = .(Batch, Genotype, phenotype)]
line_batch_long[, robust_se := mad / sqrt(n_flies)]
fwrite(line_batch_long, file.path(output_root, "line_batch_phenotypes_long.csv"))

line_raw_long <- line_batch_long[Genotype != "CS", .(
  value = mean(value),
  n_batches = .N,
  total_flies = sum(n_flies),
  min_batch_flies = min(n_flies)
), by = .(Genotype, phenotype)]
line_raw_long[, DGRP_line := sub("^DGRP", "", Genotype)]
line_raw_long[, gwas_resource_eligible := TRUE]
setcolorder(line_raw_long, c("DGRP_line", "Genotype", "phenotype", "value"))
fwrite(line_raw_long, file.path(output_root, "line_phenotypes_raw_long.csv"))

retention_by_genotype <- all_fly_qc[, .(
  planned_flies = .N,
  retained_flies = sum(retained_by_existing_gigem),
  primary_eligible_flies = sum(primary_coverage_eligible),
  removed_dead_rule = sum(removed_by_curate_dead_animals),
  removed_lifespan_rule = sum(removed_by_lifespan_rule)
), by = .(Batch, Genotype)]
retention_by_genotype[, `:=`(
  retention_rate = retained_flies / planned_flies,
  primary_eligibility_rate = primary_eligible_flies / planned_flies
)]
fwrite(retention_by_genotype, file.path(output_root, "retention_by_batch_genotype.csv"))

comparison_fields <- grep("^difference_D", names(all_existing_comparison), value = TRUE)
comparison_summary <- rbindlist(lapply(comparison_fields, function(field) {
  data.table(
    day = as.integer(sub("^difference_D", "", field)),
    n = sum(is.finite(all_existing_comparison[[field]])),
    median_difference_min = median(all_existing_comparison[[field]], na.rm = TRUE),
    median_abs_difference_min = median(abs(all_existing_comparison[[field]]), na.rm = TRUE),
    max_abs_difference_min = max(abs(all_existing_comparison[[field]]), na.rm = TRUE)
  )
}))
fwrite(comparison_summary, file.path(output_root, "existing_summary_comparison_summary.csv"))

message("Wrote minute-level reconstruction outputs to ", normalizePath(output_root))
