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

source("src/ext/resp.r")

source("src/q3/common.r")
source("src/q3/timetoevent.r")

q3_survplots_output_width <- 7
q3_survplots_output_height <- 7
q3_survplots_output_dpi <- 300
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
    vc_decline = "decline in vital capacity (<80%)",
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
    site = "site",
    sex = "sex",
    age_category = "age at onset",
    clinical_phenotype = "clinical phenotype",
    progression_category = "progression category",
    c9orf72_status = "C9orf72 status",
    sod1_status = "SOD1 status",
    tardbp_status = "TARDBP status",
    fus_status = "FUS status",
    causal_gene = "causal gene",
    site_of_onset = "site of onset"
)

q3_plots <- list(
    list(
        origins = "birth",
        events = "onset",
        output_name = "@overall/time-from-{origin}-to-{event}"
    ),
    list(
        origins = "birth",
        events = "onset",
        groups = names(q3_group_labels),
        output_name = "{group}/time-from-{origin}-to-{event}"
    ),
    list(
        origins = "onset",
        events = names(q3_event_labels),
        output_name = "@overall/time-from-{origin}-to-{event}"
    ),
    list(
        origins = "onset",
        events = names(q3_event_labels),
        groups = names(q3_group_labels),
        output_name = "{group}/time-from-{origin}-to-{event}"
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
            "ALS", "PBP", "PLS", "PMA",
            "LMN-Predominant", "UMN-Predominant"
        ))
    } else if (group == "site_of_onset") {
        data %>% filter(site_of_onset %in% c(
            "Bulbar", "Respiratory", "Spinal"
        ))
    } else {
        data
    }
}

q3_make_survival_plot <- function(data, origin, event, group = NULL, unit = "years") {
    origin_lbl <- q3_origin_labels[[origin]]
    event_lbl <- q3_event_labels[[event]]
    title <- q3_str_to_title(event_lbl)
    xlab <- str_glue("Time from {origin_lbl}, {unit}")

    data %<>%
        filter(.data$origin == .env$origin, .data$event == .env$event) %>%
        mutate(duration = duration / lubridate::duration(1, unit))

    if (is.null(group)) {
        km_fit <- survfit2(Surv(duration, status == "event") ~ 1, data)
        km_plot <- ggsurvfit(km_fit) + add_quantile()
    } else {
        group_lbl <- q3_group_labels[[group]]
        data %<>% q3_trim_survival_groups(group)
        km_fit <- survfit2(as.formula(
            str_glue("Surv(duration, status == 'event') ~ {group}")
        ), data)

        km_plot <- ggsurvfit(km_fit) +
            add_pvalue("annotation") +
            add_legend_title(q3_str_to_sentence(group_lbl))
    }

    km_plot +
        scale_ggsurvfit() +
        add_confidence_interval() +
        labs(title = title, x = xlab)
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
        epoch_unit <- p$unit %||% "years"
        for (event in p$events) {
            if (!q3_is_valid_event_from_origin(event, origin)) {
                skip_groups <- length(p$groups %||% list(NULL))
                skip_facets <- length(p$facets %||% list(NULL))
                progress_bar$tick(skip_groups * skip_facets)
                next
            }

            for (group in p$groups %||% list(NULL)) {
                plot <- q3_make_survival_plot(q3_data, origin, event, group, epoch_unit)
                output_fname <- str_glue(p$output_name) %>% with_ext(q3_survplots_output_format)
                q3_save_plot(plot, file.path("output", "q3", output_fname))
                progress_bar$tick()
            }
        }
    }
}
