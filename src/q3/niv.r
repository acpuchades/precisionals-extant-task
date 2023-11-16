source("src/ext/main.r")
source("src/ext/resp.r")
source("src/ext/alsfrs.r")

source("src/q3/common.r")

library(writexl)

event_names <- c("diagnosis", "niv", "23h_niv", "tracheostomy", "gastrostomy", "death")

event_data <- ext_main %>%
    select(
        id, site,
        date_of_birth,
        vital_status,
        niv,
        gt_23h_niv,
        tracheostomy,
        gastrostomy,
        q3_time_of("onset"),
        q3_time_of(event_names),
    ) %>%
    mutate(
        across(
            str_c("date_of_", event_names),
            ~ coalesce(
                .x - date_of_onset,
                .x - (date_of_birth + dyears(age_at_onset))
            ),
            .names = "__duration_from_date_of_onset_to_{.col}"
        ),
        across(
            str_c("age_at_", event_names), ~ coalesce(
                dyears(.x - age_at_onset),
                dyears((date_of_onset - date_of_birth) / dyears(1) - .x)
            ),
            .names = "__duration_from_age_at_onset_to_{.col}"
        ),
        across(
            starts_with("__duration_from_date_of_onset_to_date_of_"), ~ {
                event <- str_replace(cur_column(), "^__duration_from_date_of_onset_to_date_of_", "")
                coalesce(.x, get(str_c("__duration_from_age_at_onset_to_age_at_", event)))
            },
            .names = "time_from_onset_to{.col}"
        )
    ) %>%
    rename_with(~ str_replace(
        .x,
        "time_from_onset_to__duration_from_date_of_onset_to_date_of",
        "time_from_onset_to"
    )) %>%
    select(
        id, site, niv, gt_23h_niv, tracheostomy, gastrostomy,
        str_c("time_from_onset_to_", event_names)
    ) %T>%
    write_xlsx("output/q3/time-from-onset.xlsx")

event_times_from_baseline <- ext_main %>%
    select(id, date_of_birth, q3_time_of(event_names)) %>%
    inner_join(ext_baseline %>% select(id, date_of_baseline, age_at_baseline), by = "id") %>%
    mutate(
        across(
            str_c("age_at_", event_names),
            ~ coalesce(
                dyears(.x - age_at_baseline),
                dyears(.x - (date_of_baseline - date_of_birth) / dyears(1))
            ),
            .names = "time_from_baseline__{.col}"
        ),
        across(
            str_c("date_of_", event_names),
            ~ coalesce(
                as.duration(.x - date_of_baseline),
                as.duration(.x - date_of_birth + dyears(age_at_baseline))
            ),
            .names = "time_from_baseline__{.col}"
        ),
        across(
            str_c("time_from_baseline__date_of_", event_names),
            ~ {
                event <- str_replace(cur_column(), "^time_from_baseline__date_of_", "")
                coalesce(.x, get(str_c("time_from_baseline__age_at_", event)))
            }
        )
    ) %>%
    select(-starts_with("time_from_baseline__age_at")) %>%
    rename_with(~ str_replace(.x, "^time_from_baseline__date_of_", "time_from_baseline_to_"))

event_status_by_alsfrs <- ext_alsfrs_followups %>%
    transmute(
        id, site, time_from_baseline,
        niv = q12_respiratory_insufficiency %>% between(1, 3),
        gt_23h_niv = q12_respiratory_insufficiency == 1,
        tracheostomy = q12_respiratory_insufficiency == 0
    )

fvc_followups <- ext_baseline %>%
    select(id, date_of_baseline, age_at_baseline, months_from_onset) %>%
    left_join(ext_main %>% select(id, site, date_of_birth), by = "id") %>%
    inner_join(
        ext_resp %>% select(id, date_of_assessment, age_at_assessment, vc_abs, vc_rel),
        by = "id"
    ) %>%
    mutate(
        date_of_assessment = coalesce(
            date_of_assessment,
            date_of_birth + dyears(age_at_assessment)
        ),
        date_of_baseline = coalesce(
            date_of_baseline,
            date_of_birth + dyears(age_at_baseline)
        ),
        time_from_baseline = coalesce(
            date_of_assessment - date_of_baseline,
            dmonths((age_at_assessment - age_at_baseline) * 12)
        )
    ) %>%
    select(
        id, site, time_from_baseline, vc_abs, vc_rel
    )

event_status_at_vc_assessments <- event_times_from_baseline %>%
    left_join(event_status_by_alsfrs, by = "id") %>%
    bind_rows(fvc_followups) %>%
    group_by(id, site) %>%
    arrange(time_from_baseline, .by_group = TRUE) %>%
    fill(niv, gt_23h_niv, tracheostomy, vc_abs, vc_rel) %>%
    ungroup() %>%
    mutate(
        niv = niv | (
            time_from_baseline >= time_from_baseline_to_niv &
                (is.na(time_from_baseline_to_tracheostomy) |
                    time_from_baseline < time_from_baseline_to_tracheostomy)
        ),
        gt_23h_niv = gt_23h_niv | (
            time_from_baseline >= time_from_baseline_to_23h_niv &
                (is.na(time_from_baseline_to_tracheostomy) |
                    time_from_baseline < time_from_baseline_to_tracheostomy)
        ),
        tracheostomy = tracheostomy | (time_from_baseline >= time_from_baseline_to_tracheostomy)
    )

vc_ge_70_events <- event_status_at_vc_assessments %>%
    filter(vc_rel >= 70) %>%
    summarize(
        niv = any(niv, na.rm = TRUE),
        gt_23h_niv = any(gt_23h_niv, na.rm = TRUE),
        tracheostomy = any(tracheostomy, na.rm = TRUE),
        .by = c(id, site)
    )

vc_ge_60_events <- event_status_at_vc_assessments %>%
    filter(vc_rel >= 60) %>%
    summarize(
        niv = any(niv, na.rm = TRUE),
        gt_23h_niv = any(gt_23h_niv, na.rm = TRUE),
        tracheostomy = any(tracheostomy, na.rm = TRUE),
        .by = c(id, site)
    )

vc_ge_50_events <- event_status_at_vc_assessments %>%
    filter(vc_rel >= 50) %>%
    summarize(
        niv = any(niv, na.rm = TRUE),
        gt_23h_niv = any(gt_23h_niv, na.rm = TRUE),
        tracheostomy = any(tracheostomy, na.rm = TRUE),
        .by = c(id, site)
    )

vc_lt_50_events <- event_status_at_vc_assessments %>%
    filter(vc_rel < 50) %>%
    summarize(
        niv = any(niv, na.rm = TRUE),
        gt_23h_niv = any(gt_23h_niv, na.rm = TRUE),
        tracheostomy = any(tracheostomy, na.rm = TRUE),
        .by = c(id, site)
    )


q3_summary_table <- function(...) {
    t <- table(..., useNA = "ifany")
    colnames(t) %<>% replace_na("NA")
    pt <- prop.table(t, margin = 1) * 100
    colnames(pt) <- str_c("%", colnames(pt))
    cbind(t, pt)
}

sink("output/q3/event-times-at-vc.txt")

cat("# NIV PER SITE\n\n")
q3_summary_table(event_data$site, event_data$niv) %>% print()
cat("\n\n")
cat("# NIV >23h PER SITE\n\n")
q3_summary_table(event_data$site, event_data$gt_23h_niv) %>% print()
cat("\n\n")
cat("# TRACHEOSTOMY PER SITE\n\n")
q3_summary_table(event_data$site, event_data$tracheostomy) %>% print()
cat("\n\n")
cat("# GASTROSTOMY PER SITE\n\n")
q3_summary_table(event_data$site, event_data$gastrostomy) %>% print()
cat("\n\n")

cat("# NIV PER SITE at FVC ≥70%\n\n")
q3_summary_table(vc_ge_70_events$site, vc_ge_70_events$niv) %>% print()
cat("\n\n")
cat("# NIV PER SITE at FVC ≥60%\n\n")
q3_summary_table(vc_ge_60_events$site, vc_ge_60_events$niv) %>% print()
cat("\n\n")
cat("# NIV PER SITE at FVC ≥50%\n\n")
q3_summary_table(vc_ge_50_events$site, vc_ge_50_events$niv) %>% print()
cat("\n\n")
cat("# NIV PER SITE at FVC <50%\n\n")
q3_summary_table(vc_lt_50_events$site, vc_lt_50_events$niv) %>% print()
cat("\n\n")

cat("# NIV >23h PER SITE at FVC ≥70%\n\n")
q3_summary_table(vc_ge_70_events$site, vc_ge_70_events$gt_23h_niv) %>% print()
cat("\n\n")
cat("# NIV >23h PER SITE at FVC ≥60%\n\n")
q3_summary_table(vc_ge_60_events$site, vc_ge_60_events$gt_23h_niv) %>% print()
cat("\n\n")
cat("# NIV >23h PER SITE at FVC ≥50%\n\n")
q3_summary_table(vc_ge_50_events$site, vc_ge_50_events$gt_23h_niv) %>% print()
cat("\n\n")
cat("# NIV >23h PER SITE at FVC <50%\n\n")
q3_summary_table(vc_lt_50_events$site, vc_lt_50_events$gt_23h_niv) %>% print()
cat("\n\n")

cat("# TRACHEOSTOMY PER SITE at FVC ≥70%\n\n")
q3_summary_table(vc_ge_70_events$site, vc_ge_70_events$tracheostomy) %>% print()
cat("\n\n")
cat("# TRACHEOSTOMY PER SITE at FVC ≥60%\n\n")
q3_summary_table(vc_ge_60_events$site, vc_ge_60_events$tracheostomy) %>% print()
cat("\n\n")
cat("# TRACHEOSTOMY PER SITE at FVC ≥50%\n\n")
q3_summary_table(vc_ge_50_events$site, vc_ge_50_events$tracheostomy) %>% print()
cat("\n\n")
cat("# TRACHEOSTOMY PER SITE at FVC <50%\n\n")
q3_summary_table(vc_lt_50_events$site, vc_lt_50_events$tracheostomy) %>% print()
cat("\n\n")

sink()
