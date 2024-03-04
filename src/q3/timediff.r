library(ggplot2)
library(magrittr)
library(survival)

source("src/q3/impute.r")
source("src/q3/vc_niv.r")

q3_survival_data <- q3_data.imputed %>%
    mutate(site = fct_drop(site)) %>%
    q3_as_survival_data()

q3_timediff_output_dir <- file.path(q3_output_root_dir, "timediff")
dir.create(q3_timediff_output_dir, showWarnings = FALSE, recursive = TRUE)

birth_to_onset.timediff <- coxph(
    Surv(time, status) ~ year_of_diagnosis + strata(site, sex, causal_gene),
    data = q3_survival_data %>%
        q3_select_event("birth", "onset", censor_after_epochs = 100) %>%
        mutate(diagnosis_period = factor(diagnosis_period, ordered = FALSE))
)

q3_print_object(birth_to_onset.timediff, file.path(q3_timediff_output_dir, "birth-to-onset.cox.txt"))
png(file.path(q3_timediff_output_dir, "birth-to-onset.cox.png"), width = 1000, height = 1000)
q3_plot_coxph(birth_to_onset.timediff)
dev.off()

age_at_onset_vs_year_of_diagnosis.lm <- ext_main %$% lm(age_at_onset ~ year_of_diagnosis)
summary(age_at_onset_vs_year_of_diagnosis.lm) %>%
    q3_print_object(file.path(q3_timediff_output_dir, "year_of_diagnosis-vs-age_at_onset.lm.txt"))

png(file.path(q3_timediff_output_dir, "year_of_diagnosis-vs-age_at_onset.lm.png"), width = 1000, height = 1000)
par(mfrow = c(2, 2))
plot(age_at_onset_vs_year_of_diagnosis.lm)
par(mfrow = c(1, 1))
dev.off()

ext_main %>%
    drop_na(year_of_diagnosis, age_at_onset) %>%
    ggplot(aes(year_of_diagnosis, age_at_onset)) +
    geom_jitter(aes(color = site), size = 0.5, width = 0.5, height = 0.5) +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(
        title = "Year of diagnosis vs Age at onset",
        x = "Year of diagnosis", y = "Age at onset", color = "Site"
    )
ggsave(file.path(q3_timediff_output_dir, "age_at_onset-vs-year_of_diagnosis.png"))

ext_main %>%
    filter(year_of_diagnosis >= 1980) %>%
    drop_na(age_at_diagnosis, diagnosis_period) %>%
    ggplot(aes(age_at_onset, fill = diagnosis_period)) +
    geom_density(alpha = 0.3) +
    labs(
        title = "Age at onset across time",
        x = "Age at diagnosis", fill = "Time of diagnosis"
    ) +
    theme_bw()
ggsave(file.path(q3_timediff_output_dir, "age_at_onset-across-periods.histplot.png"))

ext_main %>%
    filter(year_of_diagnosis >= 1980) %>%
    drop_na(age_at_diagnosis, diagnosis_period) %>%
    ggplot(aes(x = diagnosis_period, y = age_at_onset, fill = diagnosis_period)) +
    geom_boxplot(alpha = 0.3) +
    labs(
        title = "Age at onset across time",
        y = "Age at diagnosis", x = "Time of diagnosis"
    ) +
    theme_bw() +
    theme(legend.position = "none")
ggsave(file.path(q3_timediff_output_dir, "age_at_onset-across-periods.boxplot.png"))

sink(file.path(q3_timediff_output_dir, "niv-per-period.txt"))
cat("# NIV STATUS PER PERIOD\n\n")
patients_info %>%
    filter(vital_status == "Deceased") %$%
    q3_summary_table(diagnosis_period, niv, useNA = "ifany") %>%
    print()
sink()

sink(file.path(q3_timediff_output_dir, "vc-niv-xtab-per-period.txt"))
cat("# VC AT NIV START PER PERIOD\n\n")
patients_info %>%
    filter(niv == TRUE) %$%
    q3_summary_table(diagnosis_period, vc_at_niv_interval, useNA = "ifany") %>%
    print()
cat("\n\n")

patients_info %>%
    filter(niv == TRUE) %>%
    aov(vc_at_niv ~ diagnosis_period, data = .) %>%
    summary() %>%
    print()
sink()

patients_info.vc_at_niv %>%
    drop_na(diagnosis_period) %>%
    mutate(site = fct_drop(site)) %>%
    ggplot(aes(vc_at_niv, fill = diagnosis_period)) +
    geom_density(alpha = 0.3) +
    labs(title = "Vital capacity at NIV", x = "Vital capacity (%)", y = "Density", fill = "Diagnosis period") +
    theme_bw()
ggsave(file.path(q3_timediff_output_dir, "vc-at-niv.density.png"))

patients_info %>%
    filter(niv == TRUE) %>%
    drop_na(diagnosis_period, vc_at_niv) %>%
    mutate(diagnosis_period = fct_drop(diagnosis_period)) %>%
    ggplot(aes(x = diagnosis_period, y = vc_at_niv, fill = diagnosis_period)) +
    geom_boxplot() +
    labs(title = "Vital Capacity at NIV", x = "Diagnosis period", y = "VC (%)") +
    theme_bw() +
    theme(legend.position = "none")
ggsave(file.path(q3_timediff_output_dir, "vc-at-niv.boxplot.png"))
