library(coxme)
library(magrittr)
library(survival)
library(writexl)

source("src/ext/common.r")

ext_source("src/q3/helpers.r")
ext_source("src/q3/timetoevent.r")

q3_sod1_dataset <- q3_data_w %>%
  filter(sod1_status == "Positive" | !is.na(sod1_variant_p) | altered_genes == 0) %>%
  mutate(sod1_variant_p = case_when(
    !is.na(sod1_variant_p) ~ sod1_variant_p,
    # sod1_status == "Negative" ~ "none",
    altered_genes == 0 ~ "none",
  ))

q3_sod1_dataset_recent <- q3_data_recent_w %>%
  filter(sod1_status == "Positive" | !is.na(sod1_variant_p) | altered_genes == 0) %>%
  mutate(sod1_variant_p = case_when(
    !is.na(sod1_variant_p) ~ sod1_variant_p,
    # sod1_status == "Negative" ~ "none",
    altered_genes == 0 ~ "none",
  ))

q3_sod1_output_dir <- file.path(q3_output_root_dir, "sod1")
dir.create(q3_sod1_output_dir, recursive = TRUE, showWarnings = FALSE)

sink(file.path(q3_sod1_output_dir, "entire-cohort.txt"))

cat("# SUMMARY TABLE\n\n")
sod1_variant_freq <- data.frame(with(q3_sod1_dataset, table(sod1_variant_p)))
sod1_variant_freq <- subset(sod1_variant_freq, sod1_variant_p != "none")
print(sod1_variant_freq)
cat("\n\n")

q3_sod1_onset_global <- with(
  q3_select_event(q3_sod1_dataset, "birth", "onset", censor_after_epochs = 100 * 12),
  coxme(Surv(time, status) ~ sex + site_of_onset + sod1_variant_p + (1 | site))
)

cat("# TIME TO ONSET (Entire cohort)\n\n")
q3_sod1_onset_global %>%
  summary() %>%
  print()
cat("\n\n")

q3_sod1_survival_global <- with(
  q3_select_event(q3_sod1_dataset, "onset", "death", censor_after_epochs = 10 * 12),
  coxme(Surv(time, status) ~ sex + site_of_onset + age_at_onset + sod1_variant_p + (1 | site))
)

cat("# TIME TO DEATH FROM ONSET (Entire cohort)\n\n")
q3_sod1_survival_global %>%
  summary() %>%
  print()
sink()

sink(file.path(q3_sod1_output_dir, "recent-cohort.txt"))

cat("# SUMMARY TABLE\n\n")
sod1_variant_freq_recent <- data.frame(with(q3_sod1_dataset_recent, table(sod1_variant_p)))
sod1_variant_freq_recent <- subset(sod1_variant_freq, sod1_variant_p != "none")
print(sod1_variant_freq_recent)
cat("\n\n")

q3_sod1_onset_recent <- with(
  q3_select_event(q3_sod1_dataset_recent, "birth", "onset", censor_after_epochs = 100 * 12),
  coxme(Surv(time, status) ~ sex + site_of_onset + sod1_variant_p + (1 | site))
)

cat("# TIME TO ONSET (Recent cohort)\n\n")
q3_sod1_onset_recent %>%
  summary() %>%
  print()
cat("\n\n")

png("output/q3/sod1/coxzph-onset-recent.png")
plot(cox.zph(q3_sod1_onset_recent), col = "red")
dev.off()

q3_sod1_survival_recent <- with(
  q3_select_event(q3_sod1_dataset_recent, "onset", "death", censor_after_epochs = 10 * 12),
  coxme(Surv(time, status) ~ sex + site_of_onset + age_at_onset + sod1_variant_p + (1 | site))
)

cat("# TIME TO DEATH FROM ONSET (Recent cohort)\n\n")
q3_sod1_survival_recent %>%
  summary() %>%
  print()
sink()
