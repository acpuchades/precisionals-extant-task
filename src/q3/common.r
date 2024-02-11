library(dplyr)
library(forcats)
library(stringr)

q3_output_root_dir <- "output/q3"

q3_time_of <- function(e) {
    all_of(c(str_glue("date_of_{e}"), str_glue("age_at_{e}")))
}

q3_str_restore_allcaps <- function(s) {
    s %>%
        str_replace_all("Fus", "FUS") %>%
        str_replace_all("Mitos", "MiToS") %>%
        str_replace_all("Niv", "NIV") %>%
        str_replace_all("Sod1", "SOD1") %>%
        str_replace_all("Tardbp", "TARDBP")
}

q3_str_to_sentence <- function(s) {
    s %>%
        str_to_sentence() %>%
        q3_str_restore_allcaps()
}

q3_str_to_title <- function(s) {
    s %>%
        str_to_title() %>%
        q3_str_restore_allcaps()
}

q3_is_valid_event_from_origin <- function(event, origin) {
    case_when(
        event == origin ~ FALSE,
        event != origin ~ case_match(
            origin,
            "onset" ~ event != "birth",
            "diagnosis" ~ !(event %in% c("birth", "onset")),
            .default = TRUE
        )
    )
}

q3_select_event <- function(data, origin, event, epoch = dmonths(1)) {
    data %>%
        filter(
            .data$origin == .env$origin,
            .data$event == .env$event
        ) %>%
        select(-origin, -event, -time_to_event, -time_to_loss)
}

q3_show_progress <- function(m, f) {
    message(m, "…", appendLF = FALSE)
    result <- f
    message("\r", m, "… done.")
    invisible(result)
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

q3_save_plot <- function(plt, path) {
    png(path)
    plot(plt)
    dev.off()
}
