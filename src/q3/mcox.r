library(coxme)
library(magrittr)
library(survival)
library(writexl)

source("src/ext/common.r")

ext_source("src/q3/impute.r")

q3_survival_data <- q3_data.imputed %>%
    filter(year_of_diagnosis >= 2010) %>%
    mutate(site = fct_drop(site)) %>%
    q3_as_survival_data()

q3_mcox_output_dir <- file.path(q3_output_root_dir, "mcox")
dir.create(q3_mcox_output_dir, showWarnings = FALSE, recursive = TRUE)

birth_to_onset.mcox <- coxme(
    Surv(time, status) ~
        sex + c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site),
    data = q3_survival_data %>% q3_select_event("birth", "onset", censor_after_epochs = 100)
)

png(file.path(q3_mcox_output_dir, "birth-to-onset.png"), width = 1800, height = 1800)
birth_to_onset.zph <- q3_plot_coxph(birth_to_onset.mcox)
dev.off()

q3_print_object(list(
    `Multilevel Cox Regression (Diagnosis -> Walking support)` = birth_to_onset.mcox,
    `Cox Proportional Hazard (Diagnosis -> Walking support)` = birth_to_onset.zph
), file.path(q3_mcox_output_dir, "birth-to-onset.txt"))

onset_to_diagnosis.mcox <- coxme(
    Surv(time, status) ~
        site_of_onset + sex + age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = q3_survival_data %>% q3_select_event("onset", "diagnosis", censor_after_epochs = 10)
)

png(file.path(q3_mcox_output_dir, "onset-to-diagnosis.png"), width = 1800, height = 1800)
onset_to_diagnosis.zph <- q3_plot_coxph(onset_to_diagnosis.mcox)
dev.off()

q3_print_object(list(
    `Multilevel Cox Regression (Onset -> Diagnosis)` = onset_to_diagnosis.mcox,
    `Cox Proportional Hazard (Onset -> Diagnosis)` = onset_to_diagnosis.zph
), file.path(q3_mcox_output_dir, "onset-to-diagnosis.txt"))

diagnosis_to_respiratory_onset.mcox <- coxme(
    Surv(time, status) ~
        site_of_onset + sex + age_at_onset + baseline_vc_rel +
        I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "respiratory_onset", censor_after_epochs = 10)
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-respiratory_onset.png"), width = 1800, height = 1800)
diagnosis_to_respiratory_onset.zph <- q3_plot_coxph(diagnosis_to_respiratory_onset.mcox)
dev.off()

q3_print_object(list(
    `Multilevel Cox Regression (Diagnosis -> Respiratory onset)` = diagnosis_to_respiratory_onset.mcox,
    `Cox Proportional Hazard (Diagnosis -> Respiratory onset)` = diagnosis_to_respiratory_onset.zph
), file.path(q3_mcox_output_dir, "diagnosis-to-respiratory_onset.txt"))

diagnosis_to_walking_support.mcox <- coxme(
    Surv(time, status) ~
        site_of_onset + sex + age_at_onset + baseline_vc_rel +
        I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "walking_support", censor_after_epochs = 10)
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-walking_support.png"), width = 1800, height = 1800)
diagnosis_to_walking_support.zph <- q3_plot_coxph(diagnosis_to_walking_support.mcox)
dev.off()

q3_print_object(list(
    `Multilevel Cox Regression (Diagnosis -> Walking support)` = diagnosis_to_walking_support.mcox,
    `Cox Proportional Hazard (Diagnosis -> Walking support)` = diagnosis_to_walking_support.zph
), file.path(q3_mcox_output_dir, "diagnosis-to-walking_support.txt"))

diagnosis_to_niv.mcox <- coxme(
    Surv(time, status) ~
        site_of_onset + sex + age_at_onset + baseline_vc_rel +
        I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "niv", censor_after_epochs = 10)
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-niv.png"), width = 1800, height = 1800)
diagnosis_to_niv.zph <- q3_plot_coxph(diagnosis_to_niv.mcox)
dev.off()

q3_print_object(list(
    `Multilevel Cox Regression (Diagnosis -> NIV initiation)` = diagnosis_to_niv.mcox,
    `Cox Proportional Hazard (Diagnosis -> NIV initiation)` = diagnosis_to_niv.zph
), file.path(q3_mcox_output_dir, "diagnosis-to-niv.txt"))


diagnosis_to_gastrostomy.mcox <- coxme(
    Surv(time, status) ~
        site_of_onset + sex + age_at_onset + baseline_vc_rel +
        I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "gastrostomy", censor_after_epochs = 10)
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-gastrostomy.png"), width = 1800, height = 1800)
diagnosis_to_gastrostomy.zph <- q3_plot_coxph(diagnosis_to_gastrostomy.mcox)
dev.off()

q3_print_object(list(
    `Multilevel Cox Regression (Diagnosis -> Gastrostomy placement)` = diagnosis_to_gastrostomy.mcox,
    `Cox Proportional Hazard (Diagnosis -> Gastrostomy placement)` = diagnosis_to_gastrostomy.zph
), file.path(q3_mcox_output_dir, "diagnosis-to-gastrostomy.txt"))

diagnosis_to_death.mcox <- coxme(
    Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        site_of_onset + sex + c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "death", censor_after_epochs = 10)
)

png(file.path(q3_mcox_output_dir, "diagnosis-to-death.png"), width = 1800, height = 1800)
diagnosis_to_death.zph <- q3_plot_coxph(diagnosis_to_death.mcox)
dev.off()

q3_print_object(list(
    `Multilevel Cox Regression (Diagnosis -> Death)` = diagnosis_to_death.mcox,
    `Cox Proportional Hazard (Diagnosis -> Death)` = diagnosis_to_death.zph
), file.path(q3_mcox_output_dir, "diagnosis-to-death.txt"))
