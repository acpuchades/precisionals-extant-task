library(ggplot2)
library(forcats)
library(magrittr)
library(stringr)
library(survival)
library(writexl)

source("src/ext/common.r")

ext_source("src/ext/alsfrs.r")
ext_source("src/ext/resp.r")
ext_source("src/q3/impute.r")
ext_source("src/q3/vc_niv.r")

q3_survival_data <- q3_data.imputed %>%
    filter(year_of_diagnosis >= 2010) %>%
    mutate(site = fct_drop(site)) %>%
    q3_as_survival_data()

onset_to_diagnosis.sitediff <- coxph(
    Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) +
        strata(site_of_onset, sex, causal_gene) +
        site,
    data = q3_survival_data %>%
        q3_select_event("onset", "diagnosis", censor_after_epochs = 10)
)

diagnosis_to_niv.sitediff <- coxph(
    Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        strata(site_of_onset) + strata(sex) + strata(causal_gene) +
        site,
    data = q3_survival_data %>% q3_select_event("diagnosis", "niv", censor_after_epochs = 10)
)

diagnosis_to_gastrostomy.sitediff <- coxph(
    Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        strata(site_of_onset, sex, causal_gene) + site,
    data = q3_survival_data %>% q3_select_event("diagnosis", "gastrostomy", censor_after_epochs = 10)
)

diagnosis_to_death.sitediff <- coxph(
    Surv(time, status) ~
        age_at_onset + baseline_vc_rel + diagnostic_delay +
        strata(site_of_onset, sex, causal_gene, progression_category) + site,
    data = q3_survival_data %>% q3_select_event("diagnosis", "death", censor_after_epochs = 10)
)

patients.current <- filter(patients_info, year_of_diagnosis >= 2010)
patients.current.niv <- filter(patients.current, niv == TRUE)

ext_interactive({
    q3_sitediff_output_dir <- file.path(q3_output_root_dir, "sitediff")
    dir.create(q3_sitediff_output_dir, showWarnings = FALSE, recursive = TRUE)

    q3_print_object(onset_to_diagnosis.sitediff, file.path(q3_sitediff_output_dir, "onset-to-diagnosis.txt"))
    png(file.path(q3_sitediff_output_dir, "onset-to-diagnosis.png"), width = 1000, height = 1000)
    q3_plot_coxph(onset_to_diagnosis.sitediff)
    dev.off()

    q3_print_object(diagnosis_to_niv.sitediff, file.path(q3_sitediff_output_dir, "diagnosis-to-niv.txt"))
    png(file.path(q3_sitediff_output_dir, "diagnosis-to-niv.png"), width = 1000, height = 1000)
    q3_plot_coxph(diagnosis_to_niv.sitediff)
    dev.off()

    q3_print_object(diagnosis_to_gastrostomy.sitediff, file.path(q3_sitediff_output_dir, "diagnosis-to-gastrostomy.txt"))
    png(file.path(q3_sitediff_output_dir, "diagnosis-to-gastrostomy.png"), width = 1000, height = 1000)
    q3_plot_coxph(diagnosis_to_gastrostomy.sitediff)
    dev.off()

    q3_print_object(diagnosis_to_death.sitediff, file.path(q3_sitediff_output_dir, "diagnosis-to-death.txt"))
    png(file.path(q3_sitediff_output_dir, "diagnosis-to-death.png"), width = 1000, height = 1000)
    q3_plot_coxph(diagnosis_to_death.sitediff)
    dev.off()

    sink(file.path(q3_sitediff_output_dir, "niv-per-site.txt"))
    cat("# NIV STATUS PER SITE (2010-2022)\n\n")
    patients.current %>%
        filter(vital_status == "Deceased") %$%
        q3_summary_table(site, niv, useNA = "ifany") %>%
        print()
    sink()

    sink(file.path(q3_sitediff_output_dir, "vc-niv-xtab-per-site.txt"))
    cat("# VC AT NIV START PER SITE (2010-2022)\n\n")
    patients.current.niv %$%
        q3_summary_table(site, vc_at_niv_interval, useNA = "ifany") %>%
        print()
    cat("\n\n")

    cat("# ANOVA: VC AT NIV START PER SITE (2010-2022)\n\n")
    aov(vc_at_niv ~ site, data = patients.current.niv) %>%
        summary() %>%
        print()
    sink()

    patients.current.niv %>%
        drop_na(vc_at_niv) %>%
        mutate(site = fct_drop(site)) %>%
        ggplot(aes(sample = vc_at_niv)) +
        geom_qq_line() +
        geom_qq() +
        facet_wrap(~site) +
        theme_bw()
    ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv.qqplot.png"))

    ext_main %>%
        drop_na(age_at_onset, year_of_diagnosis) %>%
        ggplot(aes(year_of_diagnosis, age_at_onset)) +
        geom_jitter(alpha = .3, size = .5, width = .5, height = .5) +
        geom_smooth(aes(color = site), method = "lm") +
        theme_bw() +
        facet_wrap(~site, scales = "free") +
        labs(
            title = "Year of diagnosis vs Age at onset among sites",
            x = "Year of diagnosis", y = "Age at onset"
        ) +
        theme(legend.position = "none")
    ggsave(file.path(q3_sitediff_output_dir, "age_at_onset-vs-year_of_diagnosis-per-site.png"))

    patients_info.niv %>%
        drop_na(vc_at_niv) %>%
        ggplot(aes(vc_at_niv, fill = site)) +
        geom_density(alpha = .3) +
        scale_fill_custom(drop = FALSE, breaks = str_c("Site ", c(1, 3, 4, 5, 7, 8, 9))) +
        labs(title = "Vital Capacity at NIV (Entire Cohort)", x = "Vital capacity (%)", y = "Density", fill = NULL) +
        theme_bw()
    ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-overall.density.png"))

    patients_info.niv %>%
        drop_na(vc_at_niv) %>%
        ggplot(aes(x = site, y = vc_at_niv, fill = site)) +
        geom_boxplot() +
        scale_fill_custom(drop = FALSE) +
        labs(title = "Vital Capacity at NIV (Entire Cohort)", x = NULL, y = "VC (%)") +
        theme_bw() +
        theme(legend.position = "none")
    ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-overall.boxplot.png"))

    patients.current.niv %>%
        drop_na(vc_at_niv) %>%
        ggplot(aes(vc_at_niv, fill = site)) +
        geom_density(alpha = .3) +
        scale_fill_custom(drop = FALSE, breaks = str_c("Site ", c(1, 5, 7:9))) +
        labs(title = "Vital Capacity at NIV (2010 – 2022)", x = "Vital capacity (%)", y = "Density", fill = NULL) +
        theme_bw()
    ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-after-2010.density.png"))

    patients.current.niv %>%
        drop_na(vc_at_niv) %>%
        ggplot(aes(x = site, y = vc_at_niv, fill = site)) +
        geom_boxplot() +
        scale_fill_custom(drop = FALSE) +
        labs(title = "Vital Capacity at NIV (2010 – 2022)", x = NULL, y = "VC (%)") +
        theme_bw() +
        theme(legend.position = "none")
    ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-after-2010.boxplot.png"))
})
