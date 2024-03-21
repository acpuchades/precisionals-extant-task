suppressPackageStartupMessages({
  library(stringr)
  library(mice)
  library(ggmice)
  library(ggplot2)
})

source("src/ext/common.r")

ext_source("src/q3/timetoevent.r")

q3_impute_m <- 20
q3_impute_minpuc <- 0
q3_impute_mincor <- .1
q3_impute_method <- "cart"

q3_impute_dir_name <- str_glue("{q3_impute_method}-m{q3_impute_m}-mincor{q3_impute_mincor}-minpuc{q3_impute_minpuc}")
q3_impute_root_dir <- file.path(q3_output_root_dir, "impute2", q3_impute_dir_name)

q3_output_imputed_data_path <- file.path(q3_impute_root_dir, "global_imp.rds")
q3_output_imputed_data_recent_path <- file.path(q3_impute_root_dir, "recent_imp.rds")

q3_fix_imputed_fields <- function(df) {
  df %>% mutate(across(c(
    riluzole_use,
    bulbar_onset,
    spinal_onset,
    cognitive_onset,
    respiratory_onset,
    generalized_onset
  ), as.logical))
}

q3_impute_data <- function(data, method = q3_impute_method, m = q3_impute_m, mincor = q3_impute_mincor, minpuc = q3_impute_minpuc, ...) {
  input <- data %>% select(-c(
    starts_with("cumhaz_onset_") & !matches("cumhaz_onset_diagnosis"),
    starts_with("date_"), starts_with(c("time", "status"))
  ))

  predmat <- quickpred(input, mincor = mincor)
  futuremice(input, predictorMatrix = predmat, method = method, m = m, mincor = mincor, minpuc = minpuc, ...) %>%
    complete(action = "long", include = TRUE) %>%
    left_join(
      data %>%
        select(id, site, starts_with("date_"), starts_with(c("time", "status"))),
      by = c("id", "site")
    ) %>%
    q3_fix_imputed_fields() %>%
    as.mids()
}

q3_save_impute_diagnostics <- function(mids, prefix) {
  ggmice(mids, aes(age_at_onset)) +
    geom_density() +
    labs(title = "Age at onset: observed vs imputed")
  ggsave(str_c(prefix, "dist-age-at-onset.png"))

  ggmice(mids, aes(diagnostic_delay)) +
    geom_density() +
    labs(title = "Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "dist-diagnostic_delay.png"))

  ggmice(mids, aes(delta_fs)) +
    geom_density() +
    labs(title = "DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "dist-deltafs.png"))

  ggmice(mids, aes(baseline_vc_rel)) +
    geom_density() +
    labs(title = "Baseline %VC: observed vs imputed")
  ggsave(str_c(prefix, "dist-baseline_vc_rel.png"))

  ggmice(mids, aes(bulbar_onset, diagnostic_delay)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Bulbar onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "bulbar_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(spinal_onset, diagnostic_delay)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Spinal onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "spinal_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(cognitive_onset, diagnostic_delay)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Cognitive onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "cognitive_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(respiratory_onset, diagnostic_delay)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Respiratory onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "respiratory_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(generalized_onset, diagnostic_delay)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Generalized onset vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "generalized_onset-vs-diagnostic_delay.png"))

  ggmice(mids, aes(clinical_phenotype, diagnostic_delay)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Clinical phenotype vs Diagnostic delay: observed vs imputed")
  ggsave(str_c(prefix, "clinical_phenotype-vs-diagnostic_delay.png"))

  ggmice(mids, aes(bulbar_onset, delta_fs)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Bulbar onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "bulbar_onset-vs-delta_fs.png"))

  ggmice(mids, aes(spinal_onset, delta_fs)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Spinal onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "spinal_onset-vs-delta_fs.png"))

  ggmice(mids, aes(cognitive_onset, delta_fs)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Cognitive onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "cognitive_onset-vs-delta_fs.png"))

  ggmice(mids, aes(respiratory_onset, delta_fs)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Respiratory onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "respiratory_onset-vs-delta_fs.png"))

  ggmice(mids, aes(generalized_onset, delta_fs)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Generalized onset vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "generalized_onset-vs-delta_fs.png"))

  ggmice(mids, aes(clinical_phenotype, delta_fs)) +
    geom_violin() +
    coord_flip() +
    labs(title = "Clinical phenotype vs DeltaFS: observed vs imputed")
  ggsave(str_c(prefix, "clinical_phenotype-vs-delta_fs.png"))

  plot_pred(mids$predictorMatrix, label = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0))
  ggsave(str_c(prefix, "predictor-matrix.png"), bg = "white")
}

if (file.exists(q3_output_imputed_data_path) & file.exists(q3_output_imputed_data_recent_path)) {
  q3_show_progress("Loading cached imputed data", {
    q3_data_w.mids <- readRDS(q3_output_imputed_data_path)
    q3_data_recent_w.mids <- readRDS(q3_output_imputed_data_recent_path)
  })
} else {
  dir.create(q3_impute_root_dir, showWarnings = FALSE, recursive = TRUE)

  q3_show_progress("Imputing data for the entire cohort", {
    q3_data_w.mids <- q3_impute_data(q3_data_w)
    saveRDS(q3_data_w.mids, q3_output_imputed_data_path)
  })

  q3_show_progress("Imputing data for the recent cohort", {
    q3_data_recent_w.mids <- q3_impute_data(q3_data_recent_w)
    saveRDS(q3_data_recent_w.mids, q3_output_imputed_data_recent_path)
  })
}

ext_interactive({
  q3_save_impute_diagnostics(q3_data_w.mids, file.path(q3_impute_root_dir, "global-"))
  q3_save_impute_diagnostics(q3_data_recent_w.mids, file.path(q3_impute_root_dir, "recent-"))
})