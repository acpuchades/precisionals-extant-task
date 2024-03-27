library(broom)
library(magrittr)
library(stringr)
library(survival)
library(writexl)

source("src/ext/common.r")

ext_source("src/q3/impute2.r")

q3_run_stcox_analysis <- function(data, prefix) {
    q3_birth_onset.stcox <- with(
        q3_select_event(data, "birth", "onset", censor_after_epochs = 100 * 12, event_required = TRUE),
        coxph(Surv(time, status) ~ year_of_diagnosis +
            strata(site, sex, c9orf72_status, sod1_status, tardbp_status, fus_status))
    )

    png(str_c(prefix, "birth-to-onset.png"), width = 1800, height = 1800)
    q3_birth_onset.stcox.zph <- q3_plot_coxph(q3_birth_onset.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_birth_onset.stcox), exponentiate = TRUE, conf.int = TRUE),
        str_c(prefix, "birth-to-onset.xlsx")
    )

    q3_onset_diagnosis.stcox <- with(
        q3_select_event(data, "onset", "diagnosis", censor_after_epochs = 10 * 12, event_required = TRUE),
        coxph(Surv(time, status) ~
            age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) +
            strata(site_of_onset, sex, c9orf72_status, sod1_status, tardbp_status, fus_status) +
            site)
    )

    png(str_c(prefix, "onset-to-diagnosis.png"), width = 1800, height = 1800)
    q3_onset_diagnosis.stcox.zph <- q3_plot_coxph(q3_onset_diagnosis.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_onset_diagnosis.stcox), exponentiate = TRUE, conf.int = TRUE),
        str_c(prefix, "onset-to-diagnosis.xlsx")
    )

    q3_diagnosis_gastrostomy.stcox <- with(
        q3_select_event(data, "diagnosis", "gastrostomy", censor_after_epochs = 10 * 12),
        coxph(Surv(time, status) ~
            age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
            strata(site_of_onset, sex, c9orf72_status, sod1_status, tardbp_status, fus_status) +
            site)
    )

    png(str_c(prefix, "diagnosis-to-gastrostomy.png"), width = 1800, height = 1800)
    q3_diagnosis_gastrostomy.stcox.zph <- q3_plot_coxph(q3_diagnosis_gastrostomy.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_diagnosis_gastrostomy.stcox), exponentiate = TRUE, conf.int = TRUE),
        str_c(prefix, "diagnosis-to-gastrostomy.xlsx")
    )

    q3_diagnosis_niv.stcox <- with(
        q3_select_event(data, "diagnosis", "niv", censor_after_epochs = 10 * 12),
        coxph(Surv(time, status) ~
            age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
            strata(site_of_onset, sex, c9orf72_status, sod1_status, tardbp_status, fus_status) +
            site)
    )

    png(str_c(prefix, "diagnosis-to-niv.png"), width = 1800, height = 1800)
    q3_diagnosis_niv.stcox.zph <- q3_plot_coxph(q3_diagnosis_niv.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_diagnosis_niv.stcox), exponentiate = TRUE, conf.int = TRUE),
        str_c(prefix, "diagnosis-to-niv.xlsx")
    )

    q3_diagnosis_death.stcox <- with(
        q3_select_event(data, "diagnosis", "death", censor_after_epochs = 10 * 12),
        coxph(Surv(time, status) ~
            age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
            strata(site_of_onset, sex, c9orf72_status, sod1_status, tardbp_status, fus_status) +
            site)
    )

    png(str_c(prefix, "diagnosis-to-death.png"), width = 1800, height = 1800)
    q3_diagnosis_death.stcox.zph <- q3_plot_coxph(q3_diagnosis_death.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_diagnosis_death.stcox), exponentiate = TRUE, conf.int = TRUE),
        str_c(prefix, "diagnosis-to-death.xlsx")
    )
}

ext_interactive({
    q3_show_progress("Performing stcox analysis on recent cohort", {
        q3_stcox_recent_output_dir <- file.path(q3_output_root_dir, "stcox.recent")
        dir.create(q3_stcox_recent_output_dir, recursive = TRUE, showWarnings = FALSE)
        q3_run_stcox_analysis(q3_data_recent_w.mids, file.path(q3_stcox_recent_output_dir, ""))
    })
})
