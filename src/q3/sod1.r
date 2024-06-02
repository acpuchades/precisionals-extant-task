library(coxme)
library(magrittr)
library(stringr)
library(survival)
library(writexl)

source("src/ext/common.r")

ext_source("src/q3/helpers.r")
ext_source("src/q3/timetoevent.r")

q3_analyze_t2e_by_sod1_variant <- function(data, prefix) {
  sink(str_c(prefix, "sod1-analysis.txt"))

  cat("# SUMMARY TABLE\n\n")
  sod1_variant_freq <- data.frame(with(data, table(sod1_variant_p)))
  sod1_variant_freq <- subset(sod1_variant_freq, sod1_variant_p != "None")
  print(sod1_variant_freq)
  cat("\n\n")

  sod1_time_to_onset <- with(
    data %>%
      q3_select_event("birth", "onset", censor_after_epochs = 100 * 12) %>%
      mutate(across(c(site, site_of_onset), fct_drop)),
    coxme(Surv(time, status) ~ sex + site_of_onset + sod1_variant_p + (1 | site))
  )

  cat("# TIME TO ONSET\n\n")
  print(tidy(sod1_time_to_onset))
  cat("\n\n")

  sod1_time_to_death <- with(
    data %>%
      q3_select_event("onset", "death", censor_after_epochs = 10 * 12) %>%
      mutate(across(c(site, site_of_onset), fct_drop)),
    coxme(Surv(time, status) ~ sex + site_of_onset + age_at_onset + sod1_variant_p + (1 | site))
  )

  cat("# TIME TO DEATH FROM ONSET\n\n")
  print(tidy(sod1_time_to_death))
  sink()
}

# Note: for some reason, we don't have data on SOD1 for patients from Trinity, so we can't
#       simply include only those marked as negative without excluding the entire Ireland
#       cohort (we have no comparator for the only variant carrier reported there).

q3_sod1_dataset_global <- q3_data_w %>%
  filter(
    !is.na(sod1_status) | !is.na(sod1_variant_p),
    site != "Cohort 7", site_of_onset != "Generalized",
  ) %>%
  mutate(sod1_variant_p = fct_relevel(case_when(
    !is.na(sod1_variant_p) ~ sod1_variant_p,
    sod1_status == "Negative" ~ "None",
  ), "None"))

q3_sod1_dataset_recent <- q3_data_recent_w %>%
  filter(
    !is.na(sod1_status) | !is.na(sod1_variant_p),
    site != "Cohort 7", site_of_onset != "Generalized",
  ) %>%
  mutate(sod1_variant_p = fct_relevel(case_when(
    !is.na(sod1_variant_p) ~ sod1_variant_p,
    sod1_status == "Negative" ~ "None",
  ), "None"))

ext_interactive({
  q3_sod1_output_dir <- file.path(q3_output_root_dir, "sod1")
  dir.create(q3_sod1_output_dir, recursive = TRUE, showWarnings = FALSE)

  q3_show_progress("Analyzing SOD1 variant data on the global cohort", {
    q3_analyze_t2e_by_sod1_variant(q3_sod1_dataset_global, file.path(q3_sod1_output_dir, "global-"))
  })

  q3_show_progress("Analyzing SOD1 variant data on the recent cohort", {
    q3_analyze_t2e_by_sod1_variant(q3_sod1_dataset_recent, file.path(q3_sod1_output_dir, "recent-"))
  })
})
