suppressPackageStartupMessages({
    library(dplyr)
    library(forcats)
    library(lazyeval)
    library(lubridate)
    library(mice)
    library(rlang)
    library(stringr)
    library(tidyr)
})

q3_output_root_dir <- "output/q3"
q3_output_data_path <- file.path(q3_output_root_dir, "timetoevent.rds")
q3_output_imputed_data_path <- file.path(q3_output_root_dir, "timetoevent_imp.rds")

source("src/q3/common.r")

q3_fix_imputed_types <- function(df) {
    df %>% mutate(
        across(c(
            riluzole_use,
            bulbar_onset,
            spinal_onset,
            cognitive_onset,
            respiratory_onset
        ), as.logical)
    )
}

q3_add_derived_variables <- function(df) {
    df %>% mutate(
        age_category = cut(age_at_onset,
            right = FALSE, ordered_result = TRUE,
            breaks = c(0, 19, 30, 40, 50, 60, 70, 80, +Inf),
            labels = c(
                "0-18", "19-29", "30-39", "40-49",
                "50-59", "60-69", "70-79", "80+"
            )
        ),
        onset_sites = (
            bulbar_onset + spinal_onset +
                cognitive_onset + respiratory_onset
        ),
        site_of_onset = q3_as_site_of_onset(case_when(
            bulbar_onset & spinal_onset ~ "Generalized",
            onset_sites > 1 ~ "Multiple",
            spinal_onset ~ "Spinal",
            bulbar_onset ~ "Bulbar",
            respiratory_onset ~ "Respiratory",
            cognitive_onset ~ "Cognitive"
        )),
        altered_genes = (
            (c9orf72_status == "Positive") +
                (sod1_status == "Positive") +
                (fus_status == "Positive") +
                (tardbp_status == "Positive")
        ),
        causal_gene = q3_as_causal_gene(case_when(
            altered_genes > 1 ~ "Multiple",
            c9orf72_status == "Positive" ~ "C9orf72",
            sod1_status == "Positive" ~ "SOD1",
            fus_status == "Positive" ~ "FUS",
            tardbp_status == "Positive" ~ "TARDBP",
            TRUE ~ "Unknown"
        )),
        diagnosis_period = factor(if_else(!is.na(year_of_diagnosis),
            {
                period_start <- year_of_diagnosis - year_of_diagnosis %% 10
                str_glue("{period_start}-{pmin(period_start+9, 2022, na.rm = TRUE)}")
            },
            NA_character_
        ), ordered = TRUE)
    )
}

q3_as_causal_gene <- function(x) {
    factor(x, levels = c(
        "Unknown", "C9orf72", "SOD1", "FUS", "TARDBP", "Multiple"
    ))
}

q3_as_site_of_onset <- function(x) {
    factor(x, levels = c(
        "Spinal", "Bulbar", "Cognitive", "Respiratory", "Generalized", "Multiple"
    ))
}

if (file.exists(q3_output_data_path)) {
    q3_show_progress("Loading cached time to event data", {
        q3_data <- readRDS(q3_output_data_path)
        q3_data.imputed <- readRDS(q3_output_imputed_data_path)
    })
} else {
    q3_show_progress("Loading the dataset", suppressMessages({
        source("src/ext/main.r")
        source("src/ext/resp.r")
        source("src/ext/staging.r")
    }))

    q3_calculate_time_to_stage <- function(data, time, stage, values) {
        stage_col <- as_label(enquo(stage))
        time_to_stage_col <- str_glue("time_from_baseline_to_{stage_col}_")
        tibble(id = unique(data$id)) %>%
            cross_join(tibble({{ stage }} := values)) %>%
            bind_rows(data %>% select(id, {{ time }}, {{ stage }})) %>%
            slice_min({{ time }}, by = c(id, {{ stage }}), n = 1, with_ties = FALSE) %>%
            group_by(id) %>%
            arrange({{ stage }}, .by_group = TRUE) %>%
            fill({{ time }}, .direction = "up") %>%
            ungroup() %>%
            pivot_wider(
                names_from = {{ stage }},
                values_from = {{ time }},
                names_prefix = time_to_stage_col
            )
    }

    q3_calculate_age_at_origin <- function(data, origin) {
        if (origin == "birth") {
            0
        } else {
            data[[str_glue("age_at_{origin}")]]
        }
    }

    q3_analyze_time_to_event <- function(data, origin, events, duration_for, censored_for = NULL) {
        result <- tibble()
        for (o in origin) {
            odata <- data %>%
                mutate(
                    .date_of_origin = .data[[str_glue("date_of_{o}")]],
                    .age_at_origin = q3_calculate_age_at_origin(data, o)
                ) %>%
                filter(!is.na(.date_of_origin) | !is.na(.age_at_origin))
            for (e in names(events)) {
                epochs_to_event <- f_eval(events[[e]], odata)
                duration_key <- if_else(exists(e, duration_for), e, ".otherwise")
                epochs_from_origin <- f_eval(duration_for[[duration_key]], odata)
                epochs_to_loss <- epochs_from_origin
                if (!is.null(censored_for) && exists(e, censored_for)) {
                    epochs_to_loss <- pmin(
                        epochs_to_loss,
                        f_eval(censored_for[[e]], odata),
                        na.rm = TRUE
                    )
                }
                duration <- pmin(epochs_from_origin, epochs_to_event, epochs_to_loss, na.rm = TRUE)
                status <- case_when(
                    duration == epochs_to_event ~ "event",
                    duration == epochs_from_origin ~ "loss",
                    duration == epochs_to_loss ~ "censored"
                )
                result <- result %>%
                    bind_rows(tibble(
                        id = odata$id,
                        origin = o,
                        event = e,
                        status = status,
                        duration = duration,
                        time_to_event = epochs_to_event,
                        time_to_loss = pmin(epochs_from_origin, epochs_to_loss, na.rm = TRUE)
                    ))
            }
        }
        result
    }

    q3_show_progress("Calculating time to King's", {
        q3_time_to_kings <- ext_kings %>%
            q3_calculate_time_to_stage(
                time = time_from_baseline,
                stage = kings, values = 0:5
            ) %>%
            left_join(ext_baseline, by = "id") %>%
            mutate(
                across(
                    starts_with("time_from_baseline_to_kings_"),
                    ~ date_of_baseline + .x,
                    .names = "date_of__{.col}"
                ),
                across(
                    starts_with("time_from_baseline_to_kings_"),
                    ~ age_at_baseline + .x / dyears(1),
                    .names = "age_at__{.col}"
                )
            ) %>%
            rename_with(~ str_replace(.x, "__time_from_baseline_to", ""))
    })

    q3_show_progress("Calculating time to MiToS", {
        q3_time_to_mitos <- ext_mitos %>%
            q3_calculate_time_to_stage(
                time = time_from_baseline,
                stage = mitos, values = 0:5
            ) %>%
            left_join(ext_baseline, by = "id") %>%
            mutate(
                across(
                    starts_with("time_from_baseline_to_mitos_"),
                    ~ date_of_baseline + .x,
                    .names = "date_of__{.col}"
                ),
                across(
                    starts_with("time_from_baseline_to_mitos_"),
                    ~ age_at_baseline + .x / dyears(1),
                    .names = "age_at__{.col}"
                )
            ) %>%
            rename_with(~ str_replace(.x, "__time_from_baseline_to", ""))
    })

    q3_show_progress("Calculating time to walking support", {
        q3_time_to_walking_support <- ext_alsfrs %>%
            filter(time_from_baseline >= ddays(0), q8_walking <= 2) %>%
            slice_min(time_from_baseline, by = "id", n = 1, with_ties = FALSE) %>%
            select(
                id,
                age_at_walking_support = "age_at_assessment",
                date_of_walking_support = "date_of_assessment"
            )
    })

    q3_show_progress("Calculating time to respiratory onset", {
        q3_time_to_respiratory_onset <- ext_alsfrs %>%
            filter(time_from_baseline >= ddays(0), q10_dyspnea <= 3 | q11_orthopnea <= 3) %>%
            slice_min(time_from_baseline, by = "id", n = 1, with_ties = FALSE) %>%
            select(
                id,
                age_at_respiratory_onset = "age_at_assessment",
                date_of_respiratory_onset = "date_of_assessment"
            )
    })

    q3_vc_assessments <- ext_resp %>%
        left_join(
            ext_main %>% select(
                id, date_of_onset, age_at_onset, date_of_diagnosis, age_at_diagnosis
            ),
            by = "id"
        ) %>%
        mutate(
            time_from_onset = coalesce(
                date_of_assessment - date_of_onset,
                dyears(age_at_assessment - age_at_onset)
            ),
            time_from_diagnosis = coalesce(
                date_of_assessment - date_of_diagnosis,
                dyears(age_at_assessment - age_at_diagnosis)
            )
        )

    q3_vc_at_baseline <- q3_vc_assessments %>%
        filter(
            time_from_onset >= dmonths(0),
            time_from_diagnosis <= dmonths(6),
            !is.na(vc_abs) | !is.na(vc_rel)
        ) %>%
        slice_min(time_from_onset / dmonths(1), by = id, n = 1, with_ties = FALSE) %>%
        select(id, baseline_vc_abs = "vc_abs", baseline_vc_rel = "vc_rel")

    q3_show_progress("Calculating time to vital capacity decline", {
        q3_time_to_vc_lt_80 <- q3_vc_assessments %>%
            filter(fvc_rel < 80 | svc_rel < 80) %>%
            slice_min(time_from_onset, n = 1, with_ties = FALSE, by = "id") %>%
            select(
                id,
                date_of_vc_lt_80 = "date_of_assessment",
                age_at_vc_lt_80 = "age_at_assessment"
            )

        q3_time_to_vc_lt_50 <- q3_vc_assessments %>%
            filter(fvc_rel < 50 | svc_rel < 50) %>%
            slice_min(time_from_onset, n = 1, with_ties = FALSE, by = "id") %>%
            select(
                id,
                date_of_vc_lt_50 = "date_of_assessment",
                age_at_vc_lt_50 = "age_at_assessment"
            )
    })

    q3_show_progress("Calculating time to NIV by ALSFRS-R", {
        q3_time_to_niv_by_alsfrs <- ext_alsfrs %>%
            filter(time_from_baseline >= ddays(0), q12_respiratory_insufficiency <= 3) %>%
            slice_min(time_from_baseline, by = "id", n = 1, with_ties = FALSE) %>%
            select(
                id,
                age_at_niv_by_alsfrs = "age_at_assessment",
                date_of_niv_by_alsfrs = "date_of_assessment"
            )
    })

    q3_show_progress("Calculating time to NIV >23h by ALSFRS-R", {
        q3_time_to_niv_23h_by_alsfrs <- ext_alsfrs %>%
            filter(time_from_baseline >= ddays(0), q12_respiratory_insufficiency <= 1) %>%
            slice_min(time_from_baseline, by = "id", n = 1, with_ties = FALSE) %>%
            select(
                id,
                age_at_niv_23h_by_alsfrs = "age_at_assessment",
                date_of_niv_23h_by_alsfrs = "date_of_assessment"
            )
    })

    q3_show_progress("Calculating time to IMV by ALSFRS-R", {
        q3_time_to_imv_by_alsfrs <- ext_alsfrs %>%
            filter(time_from_baseline >= ddays(0), q12_respiratory_insufficiency == 0) %>%
            slice_min(time_from_baseline, by = "id", n = 1, with_ties = FALSE) %>%
            select(
                id,
                age_at_imv_by_alsfrs = "age_at_assessment",
                date_of_imv_by_alsfrs = "date_of_assessment"
            )
    })

    q3_show_progress("Calculating time to last ALSFRS-R assessment", {
        q3_time_to_last_alsfrs_assessment <- ext_alsfrs %>%
            filter(time_from_baseline >= ddays(0)) %>%
            slice_max(time_from_baseline, by = "id", n = 1, with_ties = FALSE) %>%
            select(
                id,
                age_at_last_alsfrs_assessment = "age_at_assessment",
                date_of_last_alsfrs_assessment = "date_of_assessment"
            )
    })

    q3_show_progress("Calculating time to last vital capacity assessment", {
        q3_time_to_last_vc_assessment <- q3_vc_assessments %>%
            slice_max(time_from_onset, n = 1, with_ties = FALSE, by = "id") %>%
            select(
                id,
                date_of_last_vc_assessment = "date_of_assessment",
                age_at_last_vc_assessment = "age_at_assessment"
            )
    })

    # NOTE: explicit casts are necessary because pmin seems to get confused
    #       when mixing NA days (POSIXct) and dyears(x) (lubridate) types.
    q3_show_progress("Calculating time to events", {
        q3_time_to_events <- ext_main %>%
            left_join(q3_time_to_kings, by = "id") %>%
            left_join(q3_time_to_mitos, by = "id") %>%
            left_join(q3_time_to_walking_support, by = "id") %>%
            left_join(q3_time_to_respiratory_onset, by = "id") %>%
            left_join(q3_time_to_vc_lt_80, by = "id") %>%
            left_join(q3_time_to_vc_lt_50, by = "id") %>%
            left_join(q3_time_to_niv_by_alsfrs, by = "id") %>%
            left_join(q3_time_to_niv_23h_by_alsfrs, by = "id") %>%
            left_join(q3_time_to_imv_by_alsfrs, by = "id") %>%
            left_join(q3_time_to_last_vc_assessment, by = "id") %>%
            left_join(q3_time_to_last_alsfrs_assessment, by = "id") %>%
            q3_analyze_time_to_event(
                origin = c("birth", "onset", "diagnosis"),
                events = list(
                    onset = ~ coalesce(
                        as.duration(date_of_onset - .date_of_origin),
                        dyears(age_at_onset - .age_at_origin)
                    ),
                    diagnosis = ~ coalesce(
                        as.duration(date_of_diagnosis - .date_of_origin),
                        dyears(age_at_diagnosis - .age_at_origin)
                    ),
                    walking_support = ~ coalesce(
                        as.duration(date_of_walking_support - .date_of_origin),
                        dyears(age_at_walking_support - .age_at_origin)
                    ),
                    respiratory_onset = ~ coalesce(
                        as.duration(date_of_respiratory_onset - .date_of_origin),
                        dyears(age_at_respiratory_onset - .age_at_origin)
                    ),
                    vc_lt_80 = ~ coalesce(
                        as.duration(date_of_vc_lt_80 - .date_of_origin),
                        dyears(age_at_vc_lt_50 - .age_at_origin)
                    ),
                    vc_lt_50 = ~ coalesce(
                        as.duration(date_of_vc_lt_50 - .date_of_origin),
                        dyears(age_at_vc_lt_50 - .age_at_origin)
                    ),
                    ventilatory_support = ~ pmin(
                        as.duration(date_of_niv - .date_of_origin),
                        as.duration(date_of_niv_by_alsfrs - .date_of_origin),
                        as.duration(date_of_23h_niv - .date_of_origin),
                        as.duration(date_of_niv_23h_by_alsfrs - .date_of_origin),
                        as.duration(date_of_tracheostomy - .date_of_origin),
                        as.duration(date_of_imv_by_alsfrs - .date_of_origin),
                        dyears(age_at_niv - .age_at_origin),
                        dyears(age_at_niv_by_alsfrs - .age_at_origin),
                        dyears(age_at_23h_niv - .age_at_origin),
                        dyears(age_at_niv_23h_by_alsfrs - .age_at_origin),
                        dyears(age_at_tracheostomy - .age_at_origin),
                        dyears(age_at_imv_by_alsfrs - .age_at_origin),
                        na.rm = TRUE
                    ),
                    niv = ~ pmin(
                        as.duration(date_of_niv - .date_of_origin),
                        as.duration(date_of_niv_by_alsfrs - .date_of_origin),
                        as.duration(date_of_23h_niv - .date_of_origin),
                        as.duration(date_of_niv_23h_by_alsfrs - .date_of_origin),
                        dyears(age_at_niv - .age_at_origin),
                        dyears(age_at_niv_by_alsfrs - .age_at_origin),
                        dyears(age_at_23h_niv - .age_at_origin),
                        dyears(age_at_niv_23h_by_alsfrs - .age_at_origin),
                        na.rm = TRUE
                    ),
                    niv_23h = ~ pmin(
                        as.duration(date_of_23h_niv - .date_of_origin),
                        as.duration(date_of_niv_23h_by_alsfrs - .date_of_origin),
                        dyears(age_at_23h_niv - .age_at_origin),
                        dyears(age_at_niv_23h_by_alsfrs - .age_at_origin),
                        na.rm = TRUE
                    ),
                    tracheostomy = ~ pmin(
                        as.duration(date_of_tracheostomy - .date_of_origin),
                        as.duration(date_of_imv_by_alsfrs - .date_of_origin),
                        dyears(age_at_tracheostomy - .age_at_origin),
                        dyears(age_at_imv_by_alsfrs - .age_at_origin),
                        na.rm = TRUE
                    ),
                    gastrostomy = ~ coalesce(
                        as.duration(date_of_gastrostomy - .date_of_origin),
                        dyears(age_at_gastrostomy - .age_at_origin)
                    ),
                    kings_1 = ~ coalesce(
                        as.duration(date_of_kings_1 - .date_of_origin),
                        dyears(age_at_kings_1 - .age_at_origin)
                    ),
                    kings_2 = ~ coalesce(
                        as.duration(date_of_kings_2 - .date_of_origin),
                        dyears(age_at_kings_2 - .age_at_origin)
                    ),
                    kings_3 = ~ coalesce(
                        as.duration(date_of_kings_3 - .date_of_origin),
                        dyears(age_at_kings_3 - .age_at_origin)
                    ),
                    kings_4 = ~ coalesce(
                        as.duration(date_of_kings_4 - .date_of_origin),
                        dyears(age_at_kings_4 - .age_at_origin)
                    ),
                    kings_5 = ~ coalesce(
                        as.duration(date_of_kings_5 - .date_of_origin),
                        dyears(age_at_kings_5 - .age_at_origin)
                    ),
                    mitos_1 = ~ coalesce(
                        as.duration(date_of_mitos_1 - .date_of_origin),
                        dyears(age_at_mitos_1 - .age_at_origin)
                    ),
                    mitos_2 = ~ coalesce(
                        as.duration(date_of_mitos_2 - .date_of_origin),
                        dyears(age_at_mitos_2 - .age_at_origin)
                    ),
                    mitos_3 = ~ coalesce(
                        as.duration(date_of_mitos_3 - .date_of_origin),
                        dyears(age_at_mitos_3 - .age_at_origin)
                    ),
                    mitos_4 = ~ coalesce(
                        as.duration(date_of_mitos_4 - .date_of_origin),
                        dyears(age_at_mitos_4 - .age_at_origin)
                    ),
                    mitos_5 = ~ coalesce(
                        as.duration(date_of_mitos_5 - .date_of_origin),
                        dyears(age_at_mitos_5 - .age_at_origin)
                    ),
                    death = ~ coalesce(
                        as.duration(date_of_death - .date_of_origin),
                        dyears(age_at_death - .age_at_origin)
                    )
                ),
                duration_for = list(
                    vc_lt_80 = ~ coalesce(
                        as.duration(date_of_last_vc_assessment - .date_of_origin),
                        dyears(age_at_last_vc_assessment - .age_at_origin)
                    ),
                    vc_lt_50 = ~ coalesce(
                        as.duration(date_of_last_vc_assessment - .date_of_origin),
                        dyears(age_at_last_vc_assessment - .age_at_origin)
                    ),
                    walking_support = ~ coalesce(
                        as.duration(date_of_last_alsfrs_assessment - .date_of_origin),
                        dyears(age_at_last_alsfrs_assessment - .age_at_origin)
                    ),
                    respiratory_onset = ~ coalesce(
                        as.duration(date_of_last_alsfrs_assessment - .date_of_origin),
                        dyears(age_at_last_alsfrs_assessment - .age_at_origin)
                    ),
                    death = ~ pmin(
                        dyears(age_at_transfer - .age_at_origin),
                        if_else(
                            vital_status == "Alive",
                            as.duration(date_of_transfer - .date_of_origin),
                            as.duration(NA)
                        ),
                        na.rm = TRUE
                    ),
                    .otherwise = ~ pmin(
                        dyears(age_at_death - .age_at_origin),
                        dyears(age_at_transfer - .age_at_origin),
                        if_else(
                            vital_status == "Alive",
                            as.duration(date_of_transfer - .date_of_origin),
                            as.duration(date_of_death - .date_of_origin)
                        ),
                        na.rm = TRUE
                    )
                ),
                censored_for = list(
                    death = ~ pmin(
                        as.duration(date_of_tracheostomy - .date_of_origin),
                        as.duration(date_of_imv_by_alsfrs - .date_of_origin),
                        dyears(age_at_tracheostomy - .age_at_origin),
                        dyears(age_at_imv_by_alsfrs - .age_at_origin),
                        na.rm = TRUE
                    ),
                    niv = ~ pmin(
                        as.duration(date_of_tracheostomy - .date_of_origin),
                        as.duration(date_of_imv_by_alsfrs - .date_of_origin),
                        dyears(age_at_tracheostomy - .age_at_origin),
                        dyears(age_at_imv_by_alsfrs - .age_at_origin),
                        na.rm = TRUE
                    ),
                    niv_23h = ~ pmin(
                        as.duration(date_of_tracheostomy - .date_of_origin),
                        as.duration(date_of_imv_by_alsfrs - .date_of_origin),
                        dyears(age_at_tracheostomy - .age_at_origin),
                        dyears(age_at_imv_by_alsfrs - .age_at_origin),
                        na.rm = TRUE
                    )
                )
            ) %>%
            filter(
                q3_is_valid_event_from_origin(event, origin),
                duration > as.duration(0)
            )
    })

    q3_time_to_niv <- ext_main %>%
        select(id, date_of_birth, date_of_onset, age_at_onset) %>%
        left_join(
            q3_time_to_events %>%
                q3_select_event("onset", "niv") %>%
                filter(status == "event") %>%
                select(id, time_to_niv = "duration"),
            by = "id"
        ) %>%
        transmute(
            id,
            time_from_onset = time_to_niv,
            date_of_niv = coalesce(
                date_of_onset + time_to_niv,
                date_of_birth + dyears(age_at_onset) + time_to_niv
            ),
            age_at_niv = coalesce(
                age_at_onset + time_to_niv / dyears(1),
                as.duration(date_of_onset + time_to_niv - date_of_birth) / dyears(1)
            )
        )

    q3_vc_at_niv_start <- q3_time_to_niv %>%
        rename(time_from_onset_to_niv = "time_from_onset") %>%
        left_join(
            q3_vc_assessments %>%
                rename(time_from_onset_to_assessment = "time_from_onset") %>%
                select(id, time_from_onset_to_assessment, vc_abs, vc_rel),
            by = "id"
        ) %>%
        filter(time_from_onset_to_assessment <= time_from_onset_to_niv) %>%
        slice_min(time_from_onset_to_niv - time_from_onset_to_assessment, by = id, with_ties = FALSE) %>%
        select(id, vc_at_niv_abs = "vc_abs", vc_at_niv_rel = "vc_rel")

    q3_base <- ext_main.anon %>%
        left_join(ext_baseline, by = "id") %>%
        left_join(
            q3_vc_at_baseline %>%
                select(id, baseline_vc_abs, baseline_vc_rel),
            by = "id"
        ) %>%
        transmute(
            id, site, sex, baseline_vc_rel, delta_fs, progression_category, riluzole_use,
            bulbar_onset, spinal_onset, cognitive_onset, respiratory_onset,
            age_at_onset = calculated_age_at_onset, clinical_phenotype,
            date_of_diagnosis, year_of_diagnosis, diagnostic_delay = coalesce(
                (age_at_diagnosis - age_at_onset) * 12,
                as.duration(date_of_diagnosis - date_of_onset) / dmonths(1)
            ),
            c9orf72_status, sod1_status, fus_status, tardbp_status
        )

    q3_exclude_from_imputation <- c(
        "id", "site", "date_of_diagnosis", "year_of_diagnosis"
    )

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

    q3_data <- q3_base %>%
        q3_add_derived_variables() %>%
        left_join(q3_time_to_events, by = "id")

    q3_data.imputed <- q3_base.imputed %>%
        q3_add_derived_variables() %>%
        left_join(q3_time_to_events, by = "id")

    q3_show_progress("Exporting results", {
        dir.create(q3_output_root_dir, recursive = TRUE, showWarnings = FALSE)
        q3_data %>% saveRDS(q3_output_data_path)
        q3_data.imputed %>% saveRDS(q3_output_imputed_data_path)
    })
}
