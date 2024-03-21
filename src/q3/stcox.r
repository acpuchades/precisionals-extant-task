library(broom)
library(magrittr)
library(stringr)
library(survival)
library(writexl)

source("src/ext/common.r")

ext_source("src/q3/impute2.r")

q3_birth_onset.stcox <- with(
    q3_select_event(q3_data_w.mids, "birth", "onset", censor_after_epochs = 100 * 12, event_required = TRUE),
    coxph(Surv(time, status) ~ year_of_diagnosis + strata(site, sex, causal_gene))
)

q3_onset_diagnosis.stcox <- with(
    q3_select_event(q3_data_recent_w.mids, "onset", "diagnosis", censor_after_epochs = 10 * 12, event_required = TRUE),
    coxph(Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) +
        strata(site_of_onset, sex, causal_gene) +
        site)
)

q3_diagnosis_gastrostomy.stcox <- with(
    q3_select_event(q3_data_recent_w.mids, "diagnosis", "gastrostomy", censor_after_epochs = 10 * 12),
    coxph(Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        strata(site_of_onset) + strata(sex) + strata(causal_gene) +
        site)
)

q3_diagnosis_niv.stcox <- with(
    q3_select_event(q3_data_recent_w.mids, "diagnosis", "niv", censor_after_epochs = 10 * 12),
    coxph(Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        strata(site_of_onset) + strata(sex) + strata(causal_gene) +
        site)
)

q3_diagnosis_death.stcox <- with(
    q3_select_event(q3_data_recent_w.mids, "diagnosis", "death", censor_after_epochs = 10 * 12),
    coxph(Surv(time, status) ~
        age_at_onset + baseline_vc_rel + I(delta_fs^(1 / 3)) + I(diagnostic_delay^(1 / 3)) +
        strata(site_of_onset) + strata(sex) + strata(causal_gene) +
        site)
)

ext_interactive({
    q3_stcox_output_dir <- file.path(q3_output_root_dir, "stcox")
    dir.create(q3_stcox_output_dir, showWarnings = FALSE, recursive = TRUE)

    write_xlsx(
        summary(pool(q3_birth_onset.stcox), exponentiate = TRUE, conf.int = TRUE),
        file.path(q3_stcox_output_dir, "birth-to-onset.xlsx")
    )

    png(file.path(q3_stcox_output_dir, "birth-to-onset.png"), width = 1800, height = 1800)
    q3_birth_onset.stcox.zph <- q3_plot_coxph(q3_birth_onset.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_onset_diagnosis.stcox), exponentiate = TRUE, conf.int = TRUE),
        file.path(q3_stcox_output_dir, "onset-to-diagnosis.xlsx")
    )

    png(file.path(q3_stcox_output_dir, "onset-to-diagnosis.png"), width = 1800, height = 1800)
    q3_onset_diagnosis.stcox.zph <- q3_plot_coxph(q3_onset_diagnosis.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_diagnosis_gastrostomy.stcox), exponentiate = TRUE, conf.int = TRUE),
        file.path(q3_stcox_output_dir, "diagnosis-to-gastrostomy.xlsx")
    )

    png(file.path(q3_stcox_output_dir, "diagnosis-to-gastrostomy.png"), width = 1800, height = 1800)
    q3_diagnosis_gastrostomy.stcox.zph <- q3_plot_coxph(q3_diagnosis_gastrostomy.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_diagnosis_niv.stcox), exponentiate = TRUE, conf.int = TRUE),
        file.path(q3_stcox_output_dir, "diagnosis-to-niv.xlsx")
    )

    png(file.path(q3_stcox_output_dir, "diagnosis-to-niv.png"), width = 1800, height = 1800)
    q3_diagnosis_niv.stcox.zph <- q3_plot_coxph(q3_diagnosis_niv.stcox$analyses[[1]])
    dev.off()

    write_xlsx(
        summary(pool(q3_diagnosis_death.stcox), exponentiate = TRUE, conf.int = TRUE),
        file.path(q3_stcox_output_dir, "diagnosis-to-death.xlsx")
    )

    png(file.path(q3_stcox_output_dir, "diagnosis-to-death.png"), width = 1800, height = 1800)
    q3_diagnosis_death.stcox.zph <- q3_plot_coxph(q3_diagnosis_death.stcox$analyses[[1]])
    dev.off()
})
