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
q3_impute_global_data_dir <- file.path(q3_impute_root_dir, "global", "")
q3_impute_recent_data_dir <- file.path(q3_impute_root_dir, "recent", "")

q3_output_imputed_global_data_path <- file.path(q3_impute_global_data_dir, "data_imp.rds")
q3_output_imputed_recent_data_path <- file.path(q3_impute_recent_data_dir, "data_imp.rds")

q3_impute_data <- function(
    data, method = q3_impute_method, m = q3_impute_m,
    mincor = q3_impute_mincor, minpuc = q3_impute_minpuc,
    exclude = NULL, ...) {
  pmat <- quickpred(data, mincor = mincor, minpuc = minpuc, exclude = exclude.in)
  mice(data, method = method, predictorMatrix = pmat, m = m, ...)
}

q3_save_impute_diagnostics <- function(mids, prefix, exclude = NULL, exclude.out = NULL, pmat.label = FALSE, pmat.square = TRUE, pmat.rotate.x = 90, pmat.size = NULL, ...) {
  pmat <- mids$predictorMatrix
  mask <- !colnames(pmat) %in% exclude
  plot_pred(pmat[mask, mask], method = mids$method[mask], label = pmat.label, square = pmat.square) +
    theme(
      axis.text.x = element_text(angle = pmat.rotate.x, hjust = 0, size = pmat.size),
      axis.text.y = element_text(size = pmat.size),
    )
  ggsave(str_c(prefix, "predictor-matrix.png"), bg = "white", ...)

  numeric_cols <- complete(mids) %>%
    select(where(is.numeric) & !all_of(c(exclude, exclude.out))) %>%
    colnames()

  for (col in numeric_cols) {
    ggmice(mids, aes(!!sym(col), group = .imp)) +
      geom_density() +
      labs(title = col)
    ggsave(str_c(prefix, col, ".png"))
  }

  ggmice(mids, aes(age_at_onset, group = .imp)) +
    geom_density() +
    labs(title = "Age at onset: observed vs imputed")
  ggsave(str_c(prefix, "dist-age-at-onset.png"), ...)

  ggmice(mids, aes(diagnostic_delay, group = .imp)) +
    geom_density() +
    labs(title = "Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "dist-diagnostic_delay.png"), ...)

  ggmice(mids, aes(delta_fs, group = .imp)) +
    geom_density() +
    labs(title = "DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "dist-deltafs.png"), ...)

  ggmice(mids, aes(baseline_vc_rel, group = .imp)) +
    geom_density() +
    labs(title = "Baseline %VC: observed vs imputed")
  ggsave(str_c(prefix, "dist-baseline_vc_rel.png"), ...)
}

time_cols <- colnames(select(q3_data_w, starts_with("time_")))
status_cols <- colnames(select(q3_data_w, starts_with("status_")))
cumhaz_cols <- colnames(select(q3_data_w, starts_with("cumhaz_")))

onset_to_milestones_after_diagnosis_cols <- colnames(q3_data_w %>% select(
  matches("_onset_") & !matches("_onset_diagnosis")
))

cumhaz_onset_to_milestones_after_diagnosis_cols <- colnames(
  q3_data_w %>% select(starts_with("cumhaz_onset_") & !matches("cumhaz_onset_diagnosis"))
)

exclude.in <- c(
  "id", # uninformative
  time_cols, # redundant with cumhaz
  onset_to_milestones_after_diagnosis_cols, # collinear with onset-to-diagnosis + diagnosis-to-event
  "status_birth_onset", "cumhaz_birth_onset", # uninformative/collinear with age at onset
  "status_diagnosis_kings_2", "cumhaz_diagnosis_kings_2", # collinear with time to kings 1
  "status_diagnosis_mitos_2", "cumhaz_diagnosis_mitos_2", # collinear with time to mitos 1
  "status_diagnosis_kings_5", "cumhaz_diagnosis_kings_5", # uninformative/collinear with time to death
  "status_diagnosis_mitos_5", "cumhaz_diagnosis_mitos_5" # uninformative/collinear with time to death
)

exclude.out <- c(
  "year_of_diagnosis",
  onset_to_milestones_after_diagnosis_cols
)

if (!exists("q3_data_w.mids") || !exists("q3_data_recent_w.mids")) {
  if (file.exists(q3_output_imputed_global_data_path) & file.exists(q3_output_imputed_recent_data_path)) {
    q3_show_progress("Loading cached imputed data", {
      q3_data_w.mids <- readRDS(q3_output_imputed_global_data_path)
      q3_data_recent_w.mids <- readRDS(q3_output_imputed_recent_data_path)
    })
  } else {
    set.seed(1234)
    dir.create(q3_impute_root_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(q3_impute_global_data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(q3_impute_recent_data_dir, recursive = TRUE, showWarnings = FALSE)

    q3_show_progress("Imputing data for the entire cohort", {
      q3_data_w.mids <- q3_impute_data(q3_data_w, exclude = exclude.in)
      saveRDS(q3_data_w.mids, q3_output_imputed_global_data_path)
    })

    q3_show_progress("Imputing data for the recent cohort", {
      q3_data_recent_w.mids <- q3_impute_data(q3_data_recent_w, exclude = exclude.in)
      saveRDS(q3_data_recent_w.mids, q3_output_imputed_recent_data_path)
    })
  }
}

ext_interactive({
  # q3_save_impute_diagnostics(
  #  q3_data_w.mids, file.path(q3_impute_global_data_dir, ""),
  #  exclude = exclude.out, exclude.out = cumhaz_cols,
  #  pmat.label = FALSE, pmat.rotate.x = 90, pmat.size = 3
  # )

  dir.create(q3_impute_recent_data_dir, recursive = TRUE, showWarnings = FALSE)
  q3_save_impute_diagnostics(
    q3_data_recent_w.mids, file.path(q3_impute_recent_data_dir, ""),
    exclude = exclude.out, exclude.out = cumhaz_cols, pmat.label = FALSE,
    pmat.rotate.x = 90, pmat.size = 3
  )
})
