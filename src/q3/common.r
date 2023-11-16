library(dplyr)
library(forcats)
library(stringr)

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

q3_save_plot <- function(plt, path) {
    png(path)
    plot(plt)
    dev.off()
}
