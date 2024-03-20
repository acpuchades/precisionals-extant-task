suppressPackageStartupMessages({
    library(mice)
    library(xfun)
    library(ggmice)
    library(ggplot2)
    library(visdat)
})

source("src/ext/common.r")

ext_source("src/q3/timetoevent.r")

q3_impute_root_dir <- file.path(q3_output_root_dir, "impute")
q3_output_imputed_data_path <- file.path(q3_output_root_dir, "timetoevent_imp.rds")
q3_output_base_mids_data_path <- file.path(q3_output_root_dir, "timetoevent_basemids.rds")
q3_output_base_imputed_data_path <- file.path(q3_output_root_dir, "timetoevent_baseimp.rds")

q3_fix_imputed_types <- function(df) {
    df %>% mutate(across(c(
        riluzole_use,
        bulbar_onset,
        spinal_onset,
        cognitive_onset,
        respiratory_onset,
        generalized_onset
    ), as.logical))
}

q3_header_cols <- c("id", "site")
q3_exclude_from_imputation <- c("date_of_diagnosis", "year_of_diagnosis")

if (file.exists(q3_output_imputed_data_path)) {
    q3_show_progress("Loading cached imputed data", {
        q3_base.mids <- readRDS(q3_output_base_mids_data_path)
        q3_data.imputed <- readRDS(q3_output_imputed_data_path)
        q3_base.imputed <- readRDS(q3_output_base_imputed_data_path)
    })
} else {
    q3_show_progress("Imputing data", {
        input <- q3_base %>% select(-all_of(q3_exclude_from_imputation))
        q3_base.mids <- mice(input, predictorMatrix = quickpred(input))
    })

    q3_base.original <- q3_base %>%
        select(all_of(c(q3_header_cols, q3_exclude_from_imputation)))

    q3_base.imputed <- q3_base.original %>%
        right_join(
            complete(q3_base.mids) %>% q3_fix_imputed_types(),
            by = q3_header_cols
        )

    q3_data.imputed <- q3_base.original %>%
        right_join(
            complete(q3_base.mids, action = "long") %>% q3_fix_imputed_types(),
            by = q3_header_cols
        ) %>%
        q3_add_derived_variables() %>%
        left_join(q3_time_to_events, by = "id", relationship = "many-to-many")

    q3_show_progress("Exporting results", {
        q3_base.mids %>% saveRDS(q3_output_base_mids_data_path)
        q3_data.imputed %>% saveRDS(q3_output_imputed_data_path)
        q3_base.imputed %>% saveRDS(q3_output_base_imputed_data_path)
    })
}

ext_interactive({
    dir.create(q3_impute_root_dir, showWarnings = FALSE, recursive = TRUE)

    vis_dat(q3_base, facet = site)
    ggsave(file.path(q3_impute_root_dir, "missing-per-site.png"), bg = "white", width = 12, height = 7, dpi = 300)

    vis_miss(q3_base) + theme(plot.margin = margin(0, 2, 0, 0, "cm"))
    ggsave(file.path(q3_impute_root_dir, "missing-overall.png"), bg = "white", width = 10, height = 7, dpi = 300)

    ggmice(q3_base.mids, aes(age_at_onset)) +
        geom_density() +
        labs(title = "Age at onset: observed vs imputed")
    ggsave("output/q3/impute/mice-dist-age_at_onset.png")

    ggmice(q3_base.mids, aes(diagnostic_delay)) +
        geom_density() +
        labs(title = "Diagnostic delay: observed vs imputed")
    ggsave("output/q3/impute/mice-dist-diagnostic_delay.png")

    ggmice(q3_base.mids, aes(delta_fs)) +
        geom_density() +
        labs(title = "DeltaFS: observed vs imputed")
    ggsave("output/q3/impute/mice-dist-deltafs.png")

    ggmice(q3_base.mids, aes(baseline_vc_rel)) +
        geom_density() +
        labs(title = "Baseline %VC: observed vs imputed")
    ggsave("output/q3/impute/mice-dist-baseline_vc_rel.png")

    ggmice(q3_base.mids, aes(as.logical(bulbar_onset), diagnostic_delay)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Bulbar onset vs Diagnostic delay: observed vs imputed")
    ggsave("output/q3/impute/mice-bulbar_onset-vs-diagnostic_delay.png")

    ggmice(q3_base.mids, aes(as.logical(spinal_onset), diagnostic_delay)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Spinal onset vs Diagnostic delay: observed vs imputed")
    ggsave("output/q3/impute/mice-spinal_onset-vs-diagnostic_delay.png")

    ggmice(q3_base.mids, aes(as.logical(cognitive_onset), diagnostic_delay)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Cognitive onset vs Diagnostic delay: observed vs imputed")
    ggsave("output/q3/impute/mice-cognitive_onset-vs-diagnostic_delay.png")

    ggmice(q3_base.mids, aes(as.logical(respiratory_onset), diagnostic_delay)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Respiratory onset vs Diagnostic delay: observed vs imputed")
    ggsave("output/q3/impute/mice-respiratory_onset-vs-diagnostic_delay.png")

    ggmice(q3_base.mids, aes(as.logical(generalized_onset), diagnostic_delay)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Generalized onset vs Diagnostic delay: observed vs imputed")
    ggsave("output/q3/impute/mice-generalized_onset-vs-diagnostic_delay.png")

    ggmice(q3_base.mids, aes(clinical_phenotype, diagnostic_delay)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Clinical phenotype vs Diagnostic delay: observed vs imputed")
    ggsave("output/q3/impute/mice-clinical_phenotype-vs-diagnostic_delay.png")

    ggmice(q3_base.mids, aes(as.logical(bulbar_onset), delta_fs)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Bulbar onset vs DeltaFS: observed vs imputed")
    ggsave("output/q3/impute/mice-bulbar_onset-vs-delta_fs.png")

    ggmice(q3_base.mids, aes(as.logical(spinal_onset), delta_fs)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Spinal onset vs DeltaFS: observed vs imputed")
    ggsave("output/q3/impute/mice-spinal_onset-vs-delta_fs.png")

    ggmice(q3_base.mids, aes(as.logical(cognitive_onset), delta_fs)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Cognitive onset vs DeltaFS: observed vs imputed")
    ggsave("output/q3/impute/mice-cognitive_onset-vs-delta_fs.png")

    ggmice(q3_base.mids, aes(as.logical(respiratory_onset), delta_fs)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Respiratory onset vs DeltaFS: observed vs imputed")
    ggsave("output/q3/impute/mice-respiratory_onset-vs-delta_fs.png")

    ggmice(q3_base.mids, aes(as.logical(generalized_onset), delta_fs)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Generalized onset vs DeltaFS: observed vs imputed")
    ggsave("output/q3/impute/mice-generalized_onset-vs-delta_fs.png")

    ggmice(q3_base.mids, aes(clinical_phenotype, delta_fs)) +
        geom_violin() +
        coord_flip() +
        labs(title = "Clinical phenotype vs DeltaFS: observed vs imputed")
    ggsave("output/q3/impute/mice-clinical_phenotype-vs-delta_fs.png")

    plot_pred(q3_base.mids$predictorMatrix) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0))
    ggsave("output/q3/impute/mice-predictor-matrix.png")
})
