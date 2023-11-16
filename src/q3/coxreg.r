library(coxme)
library(dplyr)
library(magrittr)
library(survival)

source("src/q3/timetoevent.r")

q3_as_survival_data <- function(data, unit = "months", censor_after = NULL) {
    data %<>%
        mutate(
            status = if_else(status == "event", 1, 0),
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

survival_data <- q3_data %>%
    filter(
        clinical_phenotype != "FTD",
        bulbar_onset | spinal_onset | respiratory_onset
    ) %>%
    select(
        -causal_gene, -age_category, -site_of_onset, -clinical_phenotype
    ) %>%
    q3_as_survival_data()

survival_data.censored <- q3_data %>%
    filter(
        clinical_phenotype != "FTD",
        bulbar_onset | spinal_onset | respiratory_onset
    ) %>%
    select(
        -causal_gene, -age_category, -site_of_onset, -clinical_phenotype
    ) %>%
    q3_as_survival_data(censor_after = dyears(5))

onset_to_death <- survival_data %>%
    q3_select_event("onset", "death")

onset_to_death.censored <- survival_data %>%
    q3_select_event("onset", "death")

diagnosis_to_death <- survival_data %>%
    q3_select_event("diagnosis", "death")

diagnosis_to_death.censored <- survival_data %>%
    q3_select_event("diagnosis", "death")

onset_to_death.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + riluzole_use + baseline_vc_rel +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = onset_to_death
)

onset_to_death.censored.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + riluzole_use + baseline_vc_rel +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = onset_to_death.censored
)

diagnosis_to_death.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + riluzole_use + baseline_vc_rel +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = diagnosis_to_death
)

diagnosis_to_death.censored.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + riluzole_use + baseline_vc_rel +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = diagnosis_to_death.censored
)

onset_to_death.mcox.zph <- cox.zph(onset_to_death.mcox)
diagnosis_to_death.mcox.zph <- cox.zph(diagnosis_to_death.mcox)
onset_to_death.censored.mcox.zph <- cox.zph(onset_to_death.censored.mcox)
diagnosis_to_death.censored.mcox.zph <- cox.zph(diagnosis_to_death.censored.mcox)

sink("output/q3/survival-from-onset.mcox.txt")
cat("# MULTILEVEL COX MULTIPLE REGRESSION ANALYSIS\n\n")
print(onset_to_death.mcox)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(onset_to_death.mcox.zph)
sink()

sink("output/q3/survival-from-diagnosis.mcox.txt")
cat("# MULTILEVEL COX MULTIPLE REGRESSION ANALYSIS\n\n")
print(diagnosis_to_death.mcox)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(diagnosis_to_death.mcox.zph)
sink()

sink("output/q3/survival-from-onset-censored.mcox.txt")
cat("# MULTILEVEL COX MULTIPLE REGRESSION ANALYSIS\n\n")
print(onset_to_death.censored.mcox)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(onset_to_death.censored.mcox.zph)
sink()

sink("output/q3/survival-from-diagnosis-censored.mcox.txt")
cat("# MULTILEVEL COX MULTIPLE REGRESSION ANALYSIS\n\n")
print(diagnosis_to_death.censored.mcox)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(diagnosis_to_death.censored.mcox.zph)
sink()

q3_save_plot(onset_to_death.mcox.zph[7], "output/q3/survival-from-onset-dx_delay.mcox.zph.png") # diagnostic_delay
q3_save_plot(onset_to_death.mcox.zph[9], "output/q3/survival-from-onset-baseline_vc_rel.mcox.zph.png") # baseline_vc_rel
q3_save_plot(onset_to_death.mcox.zph[10], "output/q3/survival-from-onset-c9_status.mcox.zph.png") # c9_status
q3_save_plot(onset_to_death.mcox.zph[12], "output/q3/survival-from-onset-fus_status.mcox.zph.png") # fus_status

q3_save_plot(diagnosis_to_death.mcox.zph[8], "output/q3/survival-from-diagnosis-riluzole_use.mcox.zph.png") # riluzole_use
q3_save_plot(diagnosis_to_death.mcox.zph[9], "output/q3/survival-from-diagnosis-baseline_vc_rel.mcox.zph.png") # baseline_vc_rel
q3_save_plot(diagnosis_to_death.mcox.zph[10], "output/q3/survival-from-diagnosis-c9_status.mcox.zph.png") # c9_status
