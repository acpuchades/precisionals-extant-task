library(coxme)
library(magrittr)
library(survival)
library(writexl)

source("src/ext/common.r")

ext_source("src/q3/impute2.r")
ext_source("src/q3/helpers.r")

q3_mcox_output_dir <- file.path(q3_output_root_dir, "mcox2")
dir.create(q3_mcox_output_dir, showWarnings = FALSE, recursive = TRUE)

q3_birth_onset.mcox <- with(
  q3_select_event(q3_data_w.mids, "birth", "onset", censor_after_epochs = 100 * 12, event_required = TRUE),
  coxme(Surv(time, status) ~
    sex + site_of_onset + c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site))
)

png(file.path(q3_mcox_output_dir, "birth-to-onset.png"), width = 1800, height = 1800)
q3_birth_onset.mcox.zph <- q3_plot_coxph(q3_birth_onset.mcox$analyses[[1]])
dev.off()

write_xlsx(
  summary(pool(q3_birth_onset.mcox), exponentiate = TRUE, conf.int = TRUE),
  file.path(q3_mcox_output_dir, "birth-to-onset.xlsx")
)

q3_onset_diagnosis.mcox <- with(
  q3_select_event(q3_data_w.mids, "onset", "diagnosis", censor_after_epochs = 10 * 12, event_required = TRUE),
  coxme(Surv(time, status) ~
    site_of_onset + sex + age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) +
    c9orf72_status + sod1_status + fus_status + tardbp_status +
    (1 | site))
)

png(file.path(q3_mcox_output_dir, "onset-to-diagnosis.png"), width = 1800, height = 1800)
q3_onset_diagnosis.zph <- q3_plot_coxph(q3_onset_diagnosis.mcox$analyses[[1]])
dev.off()

write_xlsx(
  summary(pool(q3_onset_diagnosis.mcox), exponentiate = TRUE, conf.int = TRUE),
  file.path(q3_mcox_output_dir, "onset-to-diagnosis.xlsx")
)

q3_diagnosis_respiratory_onset.mcox <- with(
  q3_select_event(q3_data_w.mids, "diagnosis", "respiratory_onset", censor_after_epochs = 10 * 12),
  coxme(Surv(time, status) ~
    site_of_onset + sex + age_at_onset + baseline_vc_rel +
    I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
    c9orf72_status + sod1_status + fus_status + tardbp_status +
    (1 | site))
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-respiratory_onset.png"), width = 1800, height = 1800)
q3_diagnosis_respiratory_onset.zph <- q3_plot_coxph(q3_diagnosis_respiratory_onset.mcox$analyses[[1]])
dev.off()

write_xlsx(
  summary(pool(q3_diagnosis_respiratory_onset.mcox), exponentiate = TRUE, conf.int = TRUE),
  file.path(q3_mcox_output_dir, "diagnosis-to-respiratory_onset.xlsx")
)

q3_diagnosis_walking_support.mcox <- with(
  q3_select_event(q3_data_w.mids, "diagnosis", "walking_support", censor_after_epochs = 10 * 12),
  coxme(Surv(time, status) ~
    site_of_onset + sex + age_at_onset + baseline_vc_rel +
    I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
    c9orf72_status + sod1_status + fus_status + tardbp_status +
    (1 | site))
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-walking_support.png"), width = 1800, height = 1800)
q3_diagnosis_walking_support.zph <- q3_plot_coxph(q3_diagnosis_walking_support.mcox$analyses[[1]])
dev.off()

write_xlsx(
  summary(pool(q3_diagnosis_walking_support.mcox), exponentiate = TRUE, conf.int = TRUE),
  file.path(q3_mcox_output_dir, "diagnosis-to-walking_support.xlsx")
)

q3_diagnosis_niv.mcox <- with(
  q3_select_event(q3_data_w.mids, "diagnosis", "niv", censor_after_epochs = 10 * 12),
  coxme(Surv(time, status) ~
    site_of_onset + sex + age_at_onset + baseline_vc_rel +
    I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
    c9orf72_status + sod1_status + fus_status + tardbp_status +
    (1 | site))
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-niv.png"), width = 1800, height = 1800)
q3_diagnosis_niv.zph <- q3_plot_coxph(q3_diagnosis_niv.mcox$analyses[[1]])
dev.off()

write_xlsx(
  summary(pool(q3_diagnosis_niv.mcox), exponentiate = TRUE, conf.int = TRUE),
  file.path(q3_mcox_output_dir, "diagnosis-to-niv.xlsx")
)

q3_diagnosis_gastrostomy.mcox <- with(
  q3_select_event(q3_data_w.mids, "diagnosis", "gastrostomy", censor_after_epochs = 10 * 12),
  coxme(Surv(time, status) ~
    site_of_onset + sex + age_at_onset + baseline_vc_rel +
    I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
    c9orf72_status + sod1_status + fus_status + tardbp_status +
    (1 | site))
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-gastrostomy.png"), width = 1800, height = 1800)
q3_diagnosis_gastrostomy.zph <- q3_plot_coxph(q3_diagnosis_gastrostomy.mcox$analyses[[1]])
dev.off()

write_xlsx(
  summary(pool(q3_diagnosis_gastrostomy.mcox), exponentiate = TRUE, conf.int = TRUE),
  file.path(q3_mcox_output_dir, "diagnosis-to-gastrostomy.xlsx")
)

q3_diagnosis_death.mcox <- with(
  q3_select_event(q3_data_w.mids, "diagnosis", "death", censor_after_epochs = 10 * 12),
  coxme(Surv(time, status) ~
    site_of_onset + sex + age_at_onset + baseline_vc_rel +
    I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
    c9orf72_status + sod1_status + fus_status + tardbp_status +
    (1 | site))
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-death.png"), width = 1800, height = 1800)
q3_diagnosis_death.zph <- q3_plot_coxph(q3_diagnosis_death.mcox$analyses[[1]])
dev.off()

write_xlsx(
  summary(pool(q3_diagnosis_death.mcox), exponentiate = TRUE, conf.int = TRUE),
  file.path(q3_mcox_output_dir, "diagnosis-to-death.xlsx")
)
