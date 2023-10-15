library(coxme)
library(dplyr)
library(magrittr)
library(survival)

source("src/q3/timetoevent.r")

q3_as_survival_data <- function(data, unit = "months") {
    data %>%
        mutate(
            status = if_else(status == "event", 1, 0),
            time = duration / lubridate::duration(1, unit),
            .keep = "unused"
        ) %>%
        filter(time > 0)
}

survival_data <- q3_data %>%
    filter(
        clinical_phenotype %in% c("ALS", "PLS", "PMA"),
        bulbar_onset | spinal_onset | cognitive_onset | respiratory_onset
    ) %>%
    select(
        -causal_gene, -age_category, -site_of_onset, -clinical_phenotype
    ) %>%
    q3_as_survival_data()

onset_to_death <- survival_data %>%
    q3_select_event("onset", "death")

coxm <- coxme(
    Surv(time, status) ~ sex + bulbar_onset + spinal_onset + respiratory_onset +
        age_at_onset + diagnostic_delay + baseline_vc_rel + delta_fs + c9orf72_status +
        sod1_status + fus_status + tardbp_status + riluzole_use +
        strata(progression_category) + (1 | site),
    data = onset_to_death
)

coxm.zph <- cox.zph(coxm)
