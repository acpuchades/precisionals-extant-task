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

onset_to_death <- survival_data %>%
    q3_select_event("onset", "death")

diagnosis_to_death <- survival_data %>%
    q3_select_event("diagnosis", "death")

onset_to_niv <- survival_data %>%
    q3_select_event("onset", "niv")

diagnosis_to_niv <- survival_data %>%
    q3_select_event("diagnosis", "niv")

onset_to_death.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + riluzole_use + baseline_vc_rel +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = onset_to_death
)

diagnosis_to_death.mcox <- coxme(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + riluzole_use + baseline_vc_rel +
        c9orf72_status + sod1_status + fus_status + tardbp_status +
        (1 | site),
    data = diagnosis_to_death
)

diagnosis_to_niv.mcox.1 <- coxph(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + baseline_vc_rel +
        c9orf72_status + site,
    data = diagnosis_to_niv %>% filter(!is.na(diagnosis_period))
)

diagnosis_to_niv.mcox.2 <- coxph(
    Surv(time, status) ~ bulbar_onset + spinal_onset + respiratory_onset +
        sex + delta_fs + age_at_onset + diagnostic_delay + baseline_vc_rel +
        c9orf72_status + site + strata(diagnosis_period),
    data = diagnosis_to_niv %>% filter(!is.na(diagnosis_period))
)

library(broom)
library(ggplot2)
library(ggforestplot)

tidy(diagnosis_to_niv.mcox.2) %>%
    forestplot(name = term, logodds = TRUE, estimate = estimate, se = std.error, pvalue = p.value) +
    labs(title = "NIV from diagnosis", x = "Odds ratio (95% CI)")
ggsave("output/q3/mcox/diagnosis-to-niv.mcox.png")

onset_to_death.mcox.zph <- cox.zph(onset_to_death.mcox)
diagnosis_to_death.mcox.zph <- cox.zph(diagnosis_to_death.mcox)
diagnosis_to_niv.mcox.1.zph <- cox.zph(diagnosis_to_niv.mcox.1)
diagnosis_to_niv.mcox.2.zph <- cox.zph(diagnosis_to_niv.mcox.2)

dir.create("output/q3/mcox/zph", showWarnings = FALSE, recursive = TRUE)
sink("output/q3/mcox/onset-to-death.mcox.txt")
cat("# MULTILEVEL COX MULTIPLE REGRESSION ANALYSIS\n\n")
print(onset_to_death.mcox)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(onset_to_death.mcox.zph)
sink()

sink("output/q3/mcox/diagnosis-to-death.mcox.txt")
cat("# MULTILEVEL COX MULTIPLE REGRESSION ANALYSIS\n\n")
print(diagnosis_to_death.mcox)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(diagnosis_to_death.mcox.zph)
sink()

sink("output/q3/mcox/diagnosis-to-niv.mcox.txt")
cat("# MULTIPLE COX REGRESSION ANALYSIS (global)\n\n")
print(diagnosis_to_niv.mcox.1)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(diagnosis_to_niv.mcox.1.zph)
cat("\n")
cat("# MULTIPLE COX REGRESSION ANALYSIS (strata=diagnosis_period)\n\n")
print(diagnosis_to_niv.mcox.2)
cat("\n")
cat("## PROPORTIONALITY TESTS\n\n")
print(diagnosis_to_niv.mcox.2.zph)
sink()

q3_save_plot(onset_to_death.mcox.zph[7], "output/q3/mcox/zph/onset-to-death-dx_delay.mcox.zph.png") # diagnostic_delay
q3_save_plot(onset_to_death.mcox.zph[9], "output/q3/mcox/zph/onset-to-death-baseline_vc_rel.mcox.zph.png") # baseline_vc_rel
q3_save_plot(onset_to_death.mcox.zph[10], "output/q3/mcox/zph/onset-to-death-c9_status.mcox.zph.png") # c9_status
q3_save_plot(onset_to_death.mcox.zph[12], "output/q3/mcox/zph/onset-to-death-fus_status.mcox.zph.png") # fus_status

q3_save_plot(diagnosis_to_death.mcox.zph[8], "output/q3/mcox/zph/diagnosis-to-death-riluzole_use.mcox.zph.png") # riluzole_use
q3_save_plot(diagnosis_to_death.mcox.zph[9], "output/q3/mcox/zph/diagnosis-to-death-baseline_vc_rel.mcox.zph.png") # baseline_vc_rel
q3_save_plot(diagnosis_to_death.mcox.zph[10], "output/q3/mcox/zph/diagnosis-to-death-c9_status.mcox.zph.png") # c9_status
