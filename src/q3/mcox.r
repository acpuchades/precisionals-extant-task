library(coxme)
library(magrittr)
library(survival)
library(writexl)

source("src/q3/common.r")
source("src/q3/timetoevent.r")
source("src/q3/impute.r")

q3_as_survival_data <- function(data, unit = "years", censor_after = NULL) {
    data %<>%
        mutate(
            status = as.integer(status == "event"),
            time = duration / lubridate::duration(1, unit)
        ) %>%
        filter(time > 0)

    if (!is.null(censor_after)) {
        data %<>% mutate(
            status = if_else(status == 1 & duration <= censor_after, 1, 0),
            time = min(time, censor_after / lubridate::duration(1, unit))
        )
    }

    data
}

q3_survival_data <- q3_data.imputed %>%
    filter(
        bulbar_onset == TRUE | spinal_onset == TRUE | respiratory_onset == TRUE | cognitive_onset == TRUE,
        year_of_diagnosis >= 2010
    ) %>%
    q3_as_survival_data()

onset_to_diagnosis.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset + cognitive_onset +
        sex + delta_fs + age_at_onset + baseline_vc_rel + # clinical_phenotype +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = q3_survival_data %>% q3_select_event("onset", "diagnosis", censor_after_epochs = 10)
)

diagnosis_to_respiratory_onset.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset + cognitive_onset +
        sex + delta_fs + age_at_onset + baseline_vc_rel + diagnostic_delay + # clinical_phenotype +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "respiratory_onset", censor_after_epochs = 10)
)

diagnosis_to_walking_support.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset + cognitive_onset +
        sex + delta_fs + age_at_onset + baseline_vc_rel + diagnostic_delay + # clinical_phenotype +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "walking_support", censor_after_epochs = 10)
)

diagnosis_to_death.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset + cognitive_onset +
        sex + delta_fs + age_at_onset + baseline_vc_rel + diagnostic_delay +
        c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "death", censor_after_epochs = 10)
)

diagnosis_to_niv.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset + cognitive_onset +
        sex + delta_fs + age_at_onset + baseline_vc_rel + diagnostic_delay +
        c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "niv", censor_after_epochs = 10)
)

diagnosis_to_gastrostomy.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset + cognitive_onset +
        sex + delta_fs + age_at_onset + baseline_vc_rel + diagnostic_delay +
        c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "gastrostomy", censor_after_epochs = 10)
)

diagnosis_to_tracheostomy.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset + cognitive_onset +
        sex + delta_fs + age_at_onset + baseline_vc_rel + diagnostic_delay +
        c9orf72_status + sod1_status + fus_status + tardbp_status + (1 | site),
    data = q3_survival_data %>% q3_select_event("diagnosis", "tracheostomy", censor_after_epochs = 10)
)

q3_mcox_output_dir <- file.path(q3_output_root_dir, "mcox")
dir.create(q3_mcox_output_dir, showWarnings = FALSE, recursive = TRUE)
q3_output_model_summary(onset_to_diagnosis.mcox, file.path(q3_mcox_output_dir, "onset-to-diagnosis.txt"))
q3_output_model_summary(diagnosis_to_respiratory_onset.mcox, file.path(q3_mcox_output_dir, "diagnosis-to-respiratory_onset.txt"))
q3_output_model_summary(diagnosis_to_gastrostomy.mcox, file.path(q3_mcox_output_dir, "diagnosis-to-gastrostomy.txt"))
q3_output_model_summary(diagnosis_to_niv.mcox, file.path(q3_mcox_output_dir, "diagnosis-to-niv.txt"))
q3_output_model_summary(diagnosis_to_tracheostomy.mcox, file.path(q3_mcox_output_dir, "diagnosis-to-tracheostomy.txt"))
q3_output_model_summary(diagnosis_to_walking_support.mcox, file.path(q3_mcox_output_dir, "diagnosis-to-walking_support.txt"))
q3_output_model_summary(diagnosis_to_death.mcox, file.path(q3_mcox_output_dir, "diagnosis-to-death.txt"))
