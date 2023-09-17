library(magrittr)
library(survival)
library(coxme)

source("src/q3/timetoevent.r")

q3_events <- unique(q3_data$event)
q3_periods <- list(
    list(origin = "birth", events = c("onset"), unit = "years"),
    list(origin = "onset", events = q3_events[q3_events != "onset"], unit = "months")
)

q3_factors <- c(
    "site",
    "sex",
    "age_at_onset",
    "site_of_onset",
    "clinical_phenotype",
    "delta_fs",
    "c9orf72_status",
    "fus_status",
    "sod1_status",
    "tardbp_status"
)

q3_as_survival_data <- function(data, ..., unit = "months") {
    data %>%
        mutate(
            status = if_else(status == "event", 1, 0),
            time = duration / lubridate::duration(1, unit)
        ) %>%
        select(time, status, ...)
}

sink("output/q3/coxreg.txt")
for (period in q3_periods) {
    origin <- period$origin
    for (event in period$events) {
        data <- q3_data %>%
            q3_select_event(origin, event) %>%
            mutate(
                across(where(is.factor), fct_drop),
                sex = sex %>% fct_relevel("Male"),
                site_of_onset = site_of_onset %>% fct_relevel("Spinal"),
                clinical_phenotype = clinical_phenotype %>% fct_relevel("ALS")
            ) %>%
            q3_as_survival_data(unit = period$unit, all_of(q3_factors))

        if (event == "onset") {
            formula <- Surv(time, status) ~ sex + site_of_onset + clinical_phenotype +
                # c9orf72_status + fus_status + sod1_status + tardbp_status +
                delta_fs + (1 | site)
        } else {
            formula <- Surv(time, status) ~ sex + age_at_onset + site_of_onset + clinical_phenotype +
                # c9orf72_status + fus_status + sod1_status + tardbp_status +
                delta_fs + (1 | site)
        }

        model <- coxme(formula, data)
        cat(str_glue("# {origin} -> {event}\n\n\n"))
        print(model)
        cat("\n")
    }
}
sink()
