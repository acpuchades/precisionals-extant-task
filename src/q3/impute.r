suppressPackageStartupMessages({
    library(mice)
    library(xfun)
})

q3_output_imputed_data_path <- file.path(q3_output_root_dir, "timetoevent_imp.rds")
q3_output_base_imputed_data_path <- file.path(q3_output_root_dir, "timetoevent_baseimp.rds")

q3_fix_imputed_types <- function(df) {
    df %>% mutate(across(c(
        riluzole_use,
        bulbar_onset,
        spinal_onset,
        cognitive_onset,
        respiratory_onset
    ), as.logical))
}

q3_exclude_from_imputation <- c(
    "id", "site", "date_of_diagnosis", "year_of_diagnosis"
)

source("src/q3/timetoevent.r")

if (file.exists(q3_output_imputed_data_path)) {
    q3_show_progress("Loading cached time to event data", {
        q3_base.imputed <- readRDS(q3_output_base_imputed_data_path)
        q3_data.imputed <- readRDS(q3_output_imputed_data_path)
    })
} else {
    q3_show_progress("Imputing data", {
        predmat <- quickpred(q3_base)
        imputed <- mice(q3_base, predictorMatrix = predmat)
        q3_base.imputed <- bind_cols(
            q3_base %>% select(all_of(q3_exclude_from_imputation)),
            complete(imputed) %>%
                select(-all_of(q3_exclude_from_imputation)) %>%
                q3_fix_imputed_types()
        )
    })

    q3_data.imputed <- q3_base.imputed %>%
        q3_add_derived_variables() %>%
        left_join(q3_time_to_events, by = "id")

    q3_show_progress("Exporting results", {
        dir.create(q3_output_root_dir, recursive = TRUE, showWarnings = FALSE)
        q3_base.imputed %>% saveRDS(q3_output_base_imputed_data_path)
        q3_data.imputed %>% saveRDS(q3_output_imputed_data_path)
    })
}
