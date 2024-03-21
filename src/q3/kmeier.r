suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggsurvfit)
    library(lubridate)
    library(magrittr)
    library(progress)
    library(rlang)
    library(stringr)
    library(survival)
    library(xfun)
})

source("src/ext/common.r")

ext_source("src/ext/resp.r")
ext_source("src/q3/timetoevent.r")

q3_survplots_output_dpi <- 300
q3_survplots_output_width <- 10
q3_survplots_output_height <- 10
q3_survplots_output_format <- "png"

q3_origin_labels <- list(
    birth = "birth",
    onset = "onset",
    diagnosis = "diagnosis"
)

q3_event_labels <- list(
    onset = "onset",
    diagnosis = "diagnosis",
    walking_support = "walking support",
    respiratory_onset = "respiratory onset",
    vc_lt_80 = "decline in vital capacity (<80%)",
    vc_lt_50 = "decline in vital capacity (<50%)",
    ventilatory_support = "ventilatory support",
    niv = "NIV",
    niv_23h = "NIV >23h",
    tracheostomy = "tracheostomy",
    gastrostomy = "gastrostomy",
    kings_1 = "King's 1",
    kings_2 = "King's 2",
    kings_3 = "King's 3",
    kings_4 = "King's 4",
    kings_5 = "King's 5",
    mitos_1 = "MiToS 1",
    mitos_2 = "MiToS 2",
    mitos_3 = "MiToS 3",
    mitos_4 = "MiToS 4",
    mitos_5 = "MiToS 5",
    death = "death"
)

q3_group_labels <- list(
    site = "cohort",
    sex = "sex",
    age_category = "age at onset",
    diagnosis_period = "diagnosis period",
    clinical_phenotype = "clinical phenotype",
    progression_category = "progression category",
    c9orf72_status = "C9orf72 status",
    sod1_status = "SOD1 status",
    tardbp_status = "TARDBP status",
    fus_status = "FUS status",
    causal_gene = "causal gene",
    site_of_onset = "site of onset",
    riluzole_use = "riluzole use"
)

q3_clinical_milestones <- names(q3_event_labels)[
    !(names(q3_event_labels) %in% c("onset", "diagnosis"))
]

q3_plots <- list(
    list(
        origins = "birth",
        events = "onset",
        event_required = TRUE,
        censor_after = dyears(100),
        output_name = "kmeier/@overall/time-from-{origin}-to-{event}"
    ),
    list(
        risk_tables = TRUE,
        origins = "birth",
        events = "onset",
        event_required = TRUE,
        censor_after = dyears(100),
        output_name = "kmeier.risktables/@overall/time-from-{origin}-to-{event}"
    ),
    list(
        origins = "birth",
        events = "onset",
        event_required = TRUE,
        censor_after = dyears(100),
        groups = names(q3_group_labels),
        output_name = "kmeier/{group}/time-from-{origin}-to-{event}"
    ),
    list(
        risk_tables = TRUE,
        origins = "birth",
        events = "onset",
        event_required = TRUE,
        censor_after = dyears(100),
        groups = names(q3_group_labels),
        output_name = "kmeier.risktables/{group}/time-from-{origin}-to-{event}"
    ),
    list(
        origins = "onset",
        events = "diagnosis",
        event_required = TRUE,
        censor_after = dyears(10),
        output_name = "kmeier/@overall/time-from-{origin}-to-{event}"
    ),
    list(
        risk_tables = TRUE,
        origins = "onset",
        events = "diagnosis",
        event_required = TRUE,
        censor_after = dyears(10),
        output_name = "kmeier.risktables/@overall/time-from-{origin}-to-{event}"
    ),
    list(
        origins = "onset",
        events = "diagnosis",
        event_required = TRUE,
        censor_after = dyears(10),
        groups = names(q3_group_labels),
        output_name = "kmeier/{group}/time-from-{origin}-to-{event}"
    ),
    list(
        risk_tables = TRUE,
        origins = "onset",
        events = "diagnosis",
        event_required = TRUE,
        censor_after = dyears(10),
        groups = names(q3_group_labels),
        output_name = "kmeier.risktables/{group}/time-from-{origin}-to-{event}"
    ),
    list(
        origins = c("onset", "diagnosis"),
        events = q3_clinical_milestones,
        censor_after = dyears(10),
        output_name = "kmeier/@overall/time-from-{origin}-to-{event}"
    ),
    list(
        risk_tables = TRUE,
        origins = c("onset", "diagnosis"),
        events = q3_clinical_milestones,
        censor_after = dyears(10),
        output_name = "kmeier.risktables/@overall/time-from-{origin}-to-{event}"
    ),
    list(
        origins = c("onset", "diagnosis"),
        events = q3_clinical_milestones,
        groups = names(q3_group_labels),
        censor_after = dyears(10),
        output_name = "kmeier/{group}/time-from-{origin}-to-{event}"
    ),
    list(
        risk_tables = TRUE,
        origins = c("onset", "diagnosis"),
        events = q3_clinical_milestones,
        groups = names(q3_group_labels),
        censor_after = dyears(10),
        output_name = "kmeier.risktables/{group}/time-from-{origin}-to-{event}"
    )
)

q3_survival_plot_count <- function(plots) {
    count <- 0
    for (p in plots) {
        n_events <- length(p$events)
        n_origins <- length(p$origins)
        n_groups <- length(p$groups %||% list(NULL))
        n_plots <- n_events * n_origins * n_groups
        count <- count + n_plots
    }
    count
}

q3_trim_survival_groups <- function(data, group) {
    if (group == "causal_gene") {
        data %>% filter(causal_gene != "Multiple")
    } else if (group == "clinical_phenotype") {
        data %>% filter(clinical_phenotype %in% c(
            "ALS", "ALS/FTD", "FTD",
            "PBP", "PLS", "PMA",
            "LMN-Predominant", "UMN-Predominant",
            "Flail-Arm", "Flail-Leg"
        ))
    } else if (group == "diagnosis_period") {
        data %>% filter(diagnosis_period %in% c(
            "1990-1999", "2000-2009", "2010-2019", "2020-2022"
        ))
    } else {
        data
    }
}

q3_analyze_survival <- function(data, origin, event, group = NULL, event_required = FALSE, censor_after = NULL, risk_tables = FALSE, unit = "years") {
    origin_lbl <- q3_origin_labels[[origin]]
    event_lbl <- q3_event_labels[[event]]
    title <- q3_str_to_title(event_lbl)
    xlab <- str_glue("Time from {origin_lbl}, {unit}")

    data %<>% filter(.data$origin == .env$origin, .data$event == .env$event)
    if (!is.null(censor_after)) {
        data %<>% mutate(
            status = if_else((status == "event") & (duration <= censor_after), "event", "censored"),
            duration = pmin(duration, censor_after)
        )
    }
    data %<>% mutate(duration = duration / lubridate::duration(1, unit))

    if (event_required) {
        data %<>% filter(status == "event")
    }

    if (is.null(group)) {
        km_fit <- survfit2(Surv(duration, status == "event") ~ 1, data)
        km_plot <- ggsurvfit(km_fit) + add_quantile()
    } else {
        group_lbl <- q3_group_labels[[group]]
        data <- q3_trim_survival_groups(data, group)
        km_fit <- survfit2(as.formula(
            str_glue("Surv(duration, status == 'event') ~ {group}")
        ), data)
        km_plot <- ggsurvfit(km_fit) +
            add_legend_title(q3_str_to_sentence(group_lbl))
    }

    km_plot <- km_plot +
        scale_ggsurvfit() +
        add_confidence_interval() +
        labs(title = title, x = xlab, y = NULL)

    if (risk_tables) {
        km_plot <- km_plot + add_risktable()
    }

    list("fit" = km_fit, "plot" = km_plot)
}

q3_save_plot <- function(plot, path, width = q3_survplots_output_width,
                         height = q3_survplots_output_height,
                         dpi = q3_survplots_output_dpi) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    ggsave(path, plot, width = width, height = height, dpi = dpi)
}

progress_bar <- progress::progress_bar$new(
    format = "Exporting [:bar] :current/:total (:percent)",
    total = q3_survival_plot_count(q3_plots)
)

q3_data <- q3_data %>% mutate(across(
    c(c9orf72_status, sod1_status, tardbp_status, fus_status), ~ case_when(
        .x == "Negative" & causal_gene != "Unknown" ~ "Negative (known gene)",
        .x == "Negative" & causal_gene == "Unknown" ~ "Negative (unknown gene)",
        TRUE ~ .x
    )
))

progress_bar$tick(0)
for (p in q3_plots) {
    for (origin in p$origins) {
        for (event in p$events) {
            if (!q3_is_valid_event_from_origin(event, origin)) {
                skip_groups <- length(p$groups %||% list(NULL))
                skip_facets <- length(p$facets %||% list(NULL))
                progress_bar$tick(skip_groups * skip_facets)
                next
            }

            for (group in p$groups %||% list(NULL)) {
                results <- q3_analyze_survival(
                    q3_data, origin, event, group,
                    censor_after = p$censor_after,
                    unit = p$unit %||% "years",
                    risk_tables = p$risk_tables %||% FALSE,
                    event_required = p$event_required %||% FALSE
                )
                output_path <- file.path("output", "q3", str_glue(p$output_name))
                q3_save_plot(results$plot, output_path %>% with_ext(q3_survplots_output_format))
                q3_print_object(results$fit, output_path %>% with_ext("txt"))
                progress_bar$tick()
            }
        }
    }
}
