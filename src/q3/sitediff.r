library(ggplot2)
library(forcats)
library(stringr)
library(survival)
library(writexl)

source("src/ext/resp.r")
source("src/ext/alsfrs.r")
source("src/q3/mcox.r")

q3_sitediff_output_dir <- file.path(q3_output_root_dir, "sitediff")

birth_to_onset.coxph <- coxph(
    Surv(time, status) ~ year_of_diagnosis + strata(site),
    data = q3_survival_data %>% q3_select_event("birth", "onset", censor_after_epochs = 100)
)

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
ggsave(file.path(q3_sitediff_output_dir, "age-at-onset-vs-year-of-diagnosis.png"))

ext_main %>%
    drop_na(age_at_onset, year_of_diagnosis) %>%
    ggplot(aes(year_of_diagnosis, age_at_onset)) +
    geom_jitter(alpha = 0.3, size = 0.5, width = 0.5, height = 0.5) +
    geom_smooth(aes(color = site), method = "lm") +
    theme_bw() +
    facet_wrap(~site, scales = "free") +
    labs(
        title = "Year of diagnosis vs Age at onset among sites",
        x = "Year of diagnosis", y = "Age at onset"
    ) +
    theme(legend.position = "none")
ggsave(file.path(q3_sitediff_output_dir, "age-at-onset-vs-year-of-diagnosis-per-site.png"))

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
ggsave(file.path(q3_sitediff_output_dir, "age-at-onset-across-periods.histplot.png"))

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
ggsave(file.path(q3_sitediff_output_dir, "age-at-onset-across-periods.boxplot.png"))

onset_to_kings_3.coxph <- coxph(
    Surv(time, status) ~
        strata(site, site_of_onset, sex, causal_gene, progression_category) +
        age_at_onset + baseline_vc_rel + diagnostic_delay,
    data = q3_survival_data %>% q3_select_event("onset", "kings_3", censor_after_epochs = 10)
)

onset_to_diagnosis.coxph <- coxph(
    Surv(time, status) ~ strata(site_of_onset, sex, causal_gene, progression_category) +
        age_at_onset + baseline_vc_rel + site,
    data = q3_survival_data %>% q3_select_event("onset", "diagnosis", censor_after_epochs = 10)
)

diagnosis_to_niv.coxph <- coxph(
    Surv(time, status) ~ strata(site_of_onset, sex, causal_gene, progression_category) +
        age_at_onset + baseline_vc_rel + diagnostic_delay + site,
    data = q3_survival_data %>% q3_select_event("diagnosis", "niv", censor_after_epochs = 10)
)

diagnosis_to_gastrostomy.coxph <- coxph(
    Surv(time, status) ~
        strata(site_of_onset, sex, progression_category, causal_gene) +
        age_at_onset + baseline_vc_rel + diagnostic_delay + site,
    data = q3_survival_data %>% q3_select_event("diagnosis", "gastrostomy", censor_after_epochs = 10)
)

diagnosis_to_death.coxph <- coxph(
    Surv(time, status) ~
        strata(site_of_onset, sex, causal_gene, progression_category) +
        age_at_onset + baseline_vc_rel + diagnostic_delay + site,
    data = q3_survival_data %>% q3_select_event("diagnosis", "death", censor_after_epochs = 10)
)

dir.create(q3_sitediff_output_dir, showWarnings = FALSE, recursive = TRUE)
q3_output_model_summary(birth_to_onset.coxph, file.path(q3_sitediff_output_dir, "birth-to-onset.txt"))
q3_output_model_summary(onset_to_kings_3.coxph, file.path(q3_sitediff_output_dir, "onset-to-kings_3.txt"))
q3_output_model_summary(onset_to_diagnosis.coxph, file.path(q3_sitediff_output_dir, "onset-to-diagnosis.txt"))
q3_output_model_summary(diagnosis_to_niv.coxph, file.path(q3_sitediff_output_dir, "diagnosis-to-niv.txt"))
q3_output_model_summary(diagnosis_to_gastrostomy.coxph, file.path(q3_sitediff_output_dir, "diagnosis-to-gastrostomy.txt"))
q3_output_model_summary(diagnosis_to_death.coxph, file.path(q3_sitediff_output_dir, "diagnosis-to-death.txt"))

time_of_niv_as_reported <- ext_main.anon %>%
    select(id, niv, date_of_niv, age_at_niv)

niv_status_as_reported <- ext_main.anon %>%
    select(id, date_of_assessment = "date_of_niv", age_at_assessment = "age_at_niv", niv) %>%
    left_join(ext_baseline %>% select(id, date_of_baseline, age_at_baseline), by = "id") %>%
    mutate(.after = id, time_from_baseline = coalesce(
        as.duration(date_of_assessment - date_of_baseline),
        dyears(age_at_assessment - age_at_baseline)
    ))

niv_status_by_alsfrs <- ext_alsfrs_followups %>%
    mutate(niv = q12_respiratory_insufficiency %>% between(1, 3)) %>%
    select(id, time_from_baseline, date_of_assessment, age_at_assessment, niv)

time_of_niv_by_alsfrs <- niv_status_by_alsfrs %>%
    summarize(niv = any(niv, na.rm = TRUE), .by = id) %>%
    left_join(
        niv_status_by_alsfrs %>%
            filter(niv == TRUE) %>%
            slice_min(time_from_baseline, by = id, n = 1, with_ties = FALSE) %>%
            select(id, date_of_niv = "date_of_assessment", age_at_niv = "age_at_assessment"),
        by = "id"
    )

niv_assessments <- niv_status_as_reported %>%
    bind_rows(niv_status_by_alsfrs)

time_of_niv <- time_of_niv_as_reported %>%
    bind_rows(time_of_niv_by_alsfrs) %>%
    summarize(
        .by = id,
        niv = any(niv, na.rm = TRUE),
        age_at_niv = if_else(sum(!is.na(age_at_niv)) != 0, suppressWarnings(min(age_at_niv, na.rm = TRUE)), NA),
        date_of_niv = if_else(sum(!is.na(date_of_niv)) != 0, suppressWarnings(min(date_of_niv, na.rm = TRUE)), NA)
    )

vc_assessments <- ext_resp %>%
    select(id, date_of_assessment, age_at_assessment, vc_rel) %>%
    drop_na(vc_rel) %>%
    left_join(
        ext_baseline %>% select(id, date_of_baseline, age_at_baseline),
        by = "id"
    ) %>%
    mutate(
        .after = id, time_from_baseline = coalesce(
            as.duration(date_of_assessment - date_of_baseline),
            dyears(age_at_assessment - age_at_baseline)
        )
    ) %>%
    filter(time_from_baseline >= ddays(0)) %>%
    select(-date_of_baseline, -age_at_baseline)

vc_and_niv_assessments <- vc_assessments %>%
    bind_rows(niv_assessments %>% semi_join(vc_assessments, by = "id")) %>%
    group_by(id) %>%
    arrange(time_from_baseline, .by_group = TRUE) %>%
    fill(niv, vc_rel) %>%
    ungroup()

vc_at_niv_start <- vc_and_niv_assessments %>%
    filter(niv == TRUE) %>%
    slice_min(time_from_baseline, by = id, n = 1, with_ties = FALSE) %>%
    filter(time_from_baseline > ddays(0)) %>%
    select(id, vc_at_niv = "vc_rel") %>%
    mutate(vc_at_niv_interval = factor(case_when(
        vc_at_niv < 50 ~ "<50",
        vc_at_niv %>% between(50, 60) ~ "51-60",
        vc_at_niv %>% between(60, 70) ~ "61-70",
        vc_at_niv > 70 ~ ">70"
    ), ordered = TRUE, levels = c("<50", "51-60", "61-70", ">70")))

vc_summary <- vc_assessments %>% summarize(
    vc_min = min(vc_rel, na.rm = TRUE),
    vc_max = max(vc_rel, na.rm = TRUE),
    .by = id
)

patients_info <- q3_base %>%
    q3_add_derived_variables() %>%
    transmute(id, site, diagnosis_period) %>%
    left_join(ext_main %>% select(id, vital_status), by = "id") %>%
    left_join(time_of_niv, by = "id") %>%
    left_join(vc_summary, by = "id") %>%
    left_join(vc_at_niv_start, by = "id")

write_xlsx(patients_info, file.path(q3_sitediff_output_dir, "patients-info.xlsx"))

q3_summary_table <- function(...) {
    t <- table(...)
    mt <- margin.table(t, margin = 1)
    pt <- round(prop.table(t, margin = 1) * 100, 2)
    colnames(pt) <- str_c(colnames(pt), " (%)")
    cbind(n = mt, t, pt)
}

sink(file.path(q3_sitediff_output_dir, "niv-per-site.txt"))
cat("# NIV STATUS PER SITE\n\n")
patients_info %>%
    filter(vital_status == "Deceased") %$%
    q3_summary_table(site, niv, useNA = "ifany") %>%
    print()
sink()

sink(file.path(q3_sitediff_output_dir, "niv-per-period.txt"))
cat("# NIV STATUS PER PERIOD\n\n")
patients_info %>%
    filter(vital_status == "Deceased") %$%
    q3_summary_table(diagnosis_period, niv, useNA = "ifany") %>%
    print()
sink()

sink(file.path(q3_sitediff_output_dir, "vc-niv-xtab-per-site.txt"))
cat("# VC AT NIV START PER SITE\n\n")
patients_info %>%
    filter(niv == TRUE) %$%
    q3_summary_table(site, vc_at_niv_interval, useNA = "ifany") %>%
    print()
cat("\n\n")

patients_info %>%
    filter(niv == TRUE) %>%
    aov(vc_at_niv ~ site, data = .) %>%
    summary() %>%
    print()
sink()

sink(file.path(q3_sitediff_output_dir, "vc-niv-xtab-per-period.txt"))
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

patients_info %>%
    filter(niv == TRUE) %>%
    ggplot(aes(sample = vc_at_niv)) +
    geom_qq_line() +
    geom_qq() +
    facet_wrap(~site) +
    theme_bw()
ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv.qqplot.png"))

patients_info %>%
    filter(niv == TRUE, site != "Site 3") %>%
    mutate(site = fct_drop(site)) %>%
    ggplot(aes(vc_at_niv, fill = site)) +
    geom_density(alpha = 0.3) +
    labs(title = "Vital Capacity at NIV", x = "Vital capacity (%)", y = "Density", fill = NULL) +
    theme_bw()
ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-per-site.density.png"))

patients_info %>%
    filter(niv == TRUE, site != "Site 3") %>%
    drop_na(diagnosis_period) %>%
    mutate(site = fct_drop(site)) %>%
    ggplot(aes(vc_at_niv, fill = diagnosis_period)) +
    geom_density(alpha = 0.3) +
    labs(title = "Vital capacity at NIV", x = "Vital capacity (%)", y = "Density", fill = "Diagnosis period") +
    theme_bw()
ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-per-period.density.png"))

patients_info %>%
    filter(niv == TRUE, site != "Site 3") %>%
    drop_na(vc_at_niv) %>%
    mutate(site = fct_drop(site)) %>%
    ggplot(aes(x = site, y = vc_at_niv, fill = site)) +
    geom_boxplot() +
    labs(title = "Vital Capacity at NIV", x = NULL, y = "VC (%)") +
    theme_bw() +
    theme(legend.position = "none")
ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-per-site.boxplot.png"))

patients_info %>%
    filter(niv == TRUE) %>%
    drop_na(diagnosis_period, vc_at_niv) %>%
    mutate(diagnosis_period = fct_drop(diagnosis_period)) %>%
    ggplot(aes(x = diagnosis_period, y = vc_at_niv, fill = diagnosis_period)) +
    geom_boxplot() +
    labs(title = "Vital Capacity at NIV", x = "Diagnosis period", y = "VC (%)") +
    theme_bw() +
    theme(legend.position = "none")
ggsave(file.path(q3_sitediff_output_dir, "vc-at-niv-per-period.boxplot.png"))
