suppressPackageStartupMessages({
  library(stringr)
  library(mice)
  library(ggmice)
  library(ggplot2)
})

source("src/ext/common.r")

ext_source("src/q3/timetoevent.r")

q3_impute_m <- 30
q3_impute_minpuc <- 0
q3_impute_mincor <- .1
q3_impute_method <- "cart"

q3_impute_model_name <- str_c(
  coalesce(q3_impute_method, "default"),
  str_c("m", q3_impute_m),
  str_c("mincor", q3_impute_mincor),
  str_c("minpuc", q3_impute_minpuc),
  sep = "-"
)

q3_impute_root_dir <- file.path(q3_output_root_dir, "impute2", q3_impute_model_name)

q3_output_imputed_data_path <- file.path(q3_impute_root_dir, "global_imp.rds")
q3_output_imputed_data_recent_path <- file.path(q3_impute_root_dir, "recent_imp.rds")

q3_impute_data <- function(
    data, method = q3_impute_method, m = q3_impute_m,
    mincor = q3_impute_mincor, minpuc = q3_impute_minpuc,
    exclude.in = NULL, exclude.out = NULL, ...) {
  pmat <- quickpred(data, mincor = mincor, minpuc = minpuc, exclude = exclude.in)
  mids <- mice(data, method = method, predictorMatrix = pmat, m = m, ...)
  if (!is.null(exclude.out)) {
    mids.l <- complete(mids, action = "long", include = TRUE)
    original <- filter(mids.l, .imp == 0)
    mids.l[, exclude.out] <- original[, exclude.out]
    mids <- as.mids(mids.l)
  }
  mids
}

q3_save_impute_diagnostics <- function(mids, prefix, exclude = NULL) {
  pmat <- mids$predictorMatrix
  mask <- !(colnames(pmat) %in% exclude)
  plot_pred(pmat[mask, mask], label = FALSE) + theme(
    axis.text.x = element_text(angle = 90, hjust = 0, size = 3),
    axis.text.y = element_text(size = 3),
  )
  ggsave(str_c(prefix, "predictor-matrix.png"), bg = "white")

  ggmice(mids, aes(age_at_onset, group = .imp)) +
    geom_density() +
    labs(title = "Age at onset: observed vs imputed")
  ggsave(str_c(prefix, "dist-age-at-onset.png"))

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    labs(title = "Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "dist-diagnostic_delay.png"))

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    labs(title = "DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "dist-deltafs.png"))

  ggmice(mids, aes(baseline_vc_rel, group = .imp)) +
    geom_density() +
    labs(title = "Baseline %VC: observed vs imputed")
  ggsave(str_c(prefix, "dist-baseline_vc_rel.png"))

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    facet_wrap(~bulbar_onset, scales = "free") +
    labs(title = "Bulbar onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "bulbar_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    facet_wrap(~spinal_onset, scales = "free") +
    labs(title = "Spinal onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "spinal_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    facet_wrap(~cognitive_onset, scales = "free") +
    labs(title = "Cognitive onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "cognitive_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    facet_wrap(~respiratory_onset, scales = "free") +
    labs(title = "Respiratory onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "respiratory_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    facet_wrap(~generalized_onset, scales = "free") +
    labs(title = "Generalized onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "generalized_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    facet_wrap(~clinical_phenotype, scales = "free") +
    labs(title = "Clinical phenotype vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "clinical_phenotype-vs-diagnostic_delay.png"))

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    facet_wrap(~bulbar_onset, scales = "free") +
    labs(title = "Bulbar onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "bulbar_onset-vs-delta_fs.png"))

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    facet_wrap(~spinal_onset, scales = "free") +
    labs(title = "Spinal onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "spinal_onset-vs-delta_fs.png"))

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    facet_wrap(~cognitive_onset, scales = "free") +
    labs(title = "Cognitive onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "cognitive_onset-vs-delta_fs.png"))

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    facet_wrap(~respiratory_onset, scales = "free") +
    labs(title = "Respiratory onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "respiratory_onset-vs-delta_fs.png"))

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    facet_wrap(~generalized_onset, scales = "free") +
    labs(title = "Generalized onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "generalized_onset-vs-delta_fs.png"))

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    facet_wrap(~clinical_phenotype, scales = "free") +
    labs(title = "Clinical phenotype vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "clinical_phenotype-vs-delta_fs.png"))
}

time_cols <- colnames(select(q3_data_w, starts_with("time_")))
status_cols <- colnames(select(q3_data_w, starts_with("status_")))

onset_to_milestones_after_diagnosis_cols <- colnames(q3_data_w %>% select(
  matches("_onset_") & !matches("_onset_diagnosis")
))

cumhaz_onset_to_milestones_after_diagnosis_cols <- colnames(
  q3_data_w %>% select(starts_with("cumhaz_onset_") & !matches("cumhaz_onset_diagnosis"))
)

exclude.in <- c(
  "id", time_cols,
  onset_to_milestones_after_diagnosis_cols,
  "status_birth_onset", "cumhaz_birth_onset",
  "status_diagnosis_kings_2", "cumhaz_diagnosis_kings_2",
  "status_diagnosis_mitos_2", "cumhaz_diagnosis_mitos_2",
  "status_diagnosis_kings_5", "cumhaz_diagnosis_kings_5",
  "status_diagnosis_mitos_5", "cumhaz_diagnosis_mitos_5"
)

exclude.out <- c(
  cumhaz_onset_to_milestones_after_diagnosis_cols
)

if (!exists("q3_data_w.mids") || !exists("q3_data_recent_w.mids")) {
  if (file.exists(q3_output_imputed_data_path) & file.exists(q3_output_imputed_data_recent_path)) {
    q3_show_progress("Loading cached imputed data", {
      q3_data_w.mids <- readRDS(q3_output_imputed_data_path)
      q3_data_recent_w.mids <- readRDS(q3_output_imputed_data_recent_path)
    })
  } else {
    set.seed(1234)
    dir.create(q3_impute_root_dir, showWarnings = FALSE, recursive = TRUE)

    q3_show_progress("Imputing data for the entire cohort", {
      q3_data_w.mids <- q3_impute_data(q3_data_w, exclude.in = exclude.in, exclude.out = exclude.out)
      saveRDS(q3_data_w.mids, q3_output_imputed_data_path)
    })

    q3_show_progress("Imputing data for the recent cohort", {
      q3_data_recent_w.mids <- q3_impute_data(q3_data_recent_w, exclude.in = exclude.in, exclude.out = exclude.out)
      saveRDS(q3_data_recent_w.mids, q3_output_imputed_data_recent_path)
    })
  }
}

ext_interactive({
  q3_save_impute_diagnostics(q3_data_w.mids, file.path(q3_impute_root_dir, "global-"), exclude = exclude.out)
  q3_save_impute_diagnostics(q3_data_recent_w.mids, file.path(q3_impute_root_dir, "recent-"), exclude = exclude.out)
})
