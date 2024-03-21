library(dplyr)
library(forcats)
library(ggplot2)
library(ggsci)
library(stringr)

source("src/ext/common.r")

ext_source("src/ext/alsfrs.r")
ext_source("src/q3/helpers.r")

q3_output_root_dir <- "output/q3"

scale_fill_custom <- scale_fill_frontiers
scale_colour_custom <- scale_color_custom <- scale_color_frontiers

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
            "birth" ~ event == "onset",
            "onset" ~ event != "birth",
            "diagnosis" ~ !(event %in% c("birth", "onset")),
            .default = TRUE
        )
    )
}

q3_as_survival_data <- function(data, unit = "years", censor_after = NULL) {
    data %<>%
        mutate(
            status = as.integer(status == "event"),
            time = duration / lubridate::duration(1, unit)
        ) %>%
        filter(time > 0)

    if (!is.null(censor_after)) {
        data %<>% mutate(
            status = if_else(status == 1 & duration <= censor_after, 1, 0),
            time = min(time, censor_after / lubridate::duration(1, unit))
        )
    }

    data
}

q3_select_event <- function(x, ...) {
    UseMethod("q3_select_event", x)
}

q3_select_event.default <- function(data, origin, event, epoch = dmonths(1), censor_after_epochs = NULL) {
    data %<>%
        filter(.data$origin == .env$origin, .data$event == .env$event) %>%
        select(-origin, -event, -time_to_event, -time_to_loss)

    if (!is.null(censor_after_epochs)) {
        data %<>% mutate(
            status = as.integer((status == 1) & (time <= censor_after_epochs)),
            time = pmin(time, censor_after_epochs)
        )
    }

    data
}

q3_select_event.mids <- function(mids, origin, event, event_required = FALSE, censor_after_epochs = NULL) {
    data <- complete(mids, action = "long", include = TRUE) %>%
        rename(
            time = str_glue("time_{origin}_{event}"),
            status = str_glue("status_{origin}_{event}")
        )

    if (!is.null(censor_after_epochs)) {
        data <- data %>% mutate(
            status = as.numeric((status == 1) & (time <= censor_after_epochs)),
            time = pmin(time, censor_after_epochs)
        )
    }

    if (event_required) {
        data <- mutate(data, status = na_if(status, status != 1))
    }

    as.mids(data)
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
            generalized_onset ~ "Generalized",
            onset_sites > 1 ~ "Generalized",
            bulbar_onset ~ "Bulbar",
            spinal_onset ~ "Spinal",
            cognitive_onset ~ "Cognitive",
            respiratory_onset ~ "Respiratory",
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
        ), ordered = TRUE),
        progression_category = ext_as_progression_category(
            case_when(
                delta_fs < ext_baseline_deltafs_p25 ~ "Slow",
                delta_fs %>% between(
                    ext_baseline_deltafs_p25,
                    ext_baseline_deltafs_p75
                ) ~ "Intermediate",
                delta_fs > ext_baseline_deltafs_p75 ~ "Fast"
            )
        )
    )
}

q3_summary_table <- function(...) {
    t <- table(...)
    mt <- margin.table(t, margin = 1)
    pt <- round(prop.table(t, margin = 1) * 100, 2)
    colnames(pt) <- str_c(colnames(pt), " (%)")
    cbind(n = mt, t, pt)
}

q3_plot_coxph <- function(model, ...) {
    zph <- cox.zph(model)
    nfigs <- nrow(zph$table) - 1
    nrows <- round(sqrt(nfigs))
    ncols <- nrows + as.integer(nrows^2 < nfigs)
    par(mfrow = c(nrows, ncols))
    for (i in 1:(nrow(zph$table) - 1)) {
        plot(zph[i], col = "red", ...)
    }
    par(mfrow = c(1, 1))
    zph
}

q3_print_object <- function(x, path) {
    sink(path)
    if (first(class(x)) == "list") {
        for (key in names(x)) {
            if (key != "") {
                cat(str_glue("# {key}\n"))
                cat("\n\n")
            }
            print(x[[key]])
            cat("\n")
        }
    } else {
        print(x)
    }
    sink()
}

q3_save_plot <- function(plt, path) {
    png(path)
    plot(plt)
    dev.off()
}
