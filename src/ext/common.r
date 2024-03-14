suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(stringr)
    library(tidyr)
})

ext_source <- function(path) {
    if (!exists("ext_source_level__")) {
        ext_source_level__ <- 0
        ext_source_interactive__ <- FALSE
    } else {
        ext_source_level__ <- ext_source_level__ + 1
    }

    source(path)

    ext_source_level__ <- ext_source_level__ - 1
    if (ext_source_level__ == 0) {
        rm(ext_source_level__)
        rm(ext_source_interactive__)
    }
}

ext_interactive <- function(f) {
    if (!exists("ext_source_interactive__") || ext_source_interactive__ == TRUE) {
        eval(f)
    }
}

ext_load_data <- function(path, ...) {
    data_dir <- Sys.getenv("PALS_EXTANT_DATADIR", unset = "./data")
    read_excel(file.path(data_dir, path), na = c(
        "Missing", "N/A", "NA", "Unknown"
    ), ...) %>% ext_normalize_names()
}

ext_normalize_names <- function(xs) {
    xs %>%
        rename_with(~ .x %>%
            str_to_lower() %>%
            str_replace_all("<", "lt_") %>%
            str_replace_all(">", "gt_") %>%
            str_replace_all("[^A-Za-z0-9]+", "_") %>%
            str_replace_all("(^_)|(_$)", "") %>%
            str_replace_all("^([^A-Za-z])", "x\\1"))
}

ext_rows_update <- function(data, t, ...) {
    rows_update(data, t, ..., unmatched = "ignore")
}

ext_rows_delete <- function(data, t, ...) {
    rows_delete(data, t, ..., unmatched = "ignore")
}

ext_parse_boolean <- function(x) {
    x %>% case_match(
        c("Yes", "yes") ~ TRUE,
        c("No", "no") ~ FALSE
    )
}

ext_apply_corrections <- function(main, corrected) {
    corrected_cols <- corrected %>%
        select(starts_with("corrected_")) %>%
        colnames() %>%
        str_replace("corrected_", "")

    main %>%
        left_join(corrected, by = "id") %>%
        mutate(across(all_of(corrected_cols), ~ (function(col, prev) {
            corrected <- get(str_c("corrected_", col))
            coalesce(corrected, prev)
        })(cur_column(), .x))) %>%
        select(-starts_with("corrected_"))
}

ext_as_sex <- function(x) {
    factor(x, levels = c("Male", "Female"))
}

ext_as_gene_status <- function(x) {
    factor(x, levels = c("Negative", "Positive"))
}

ext_logical_to_factor <- function(x, when_true = "Yes", when_false = "No") {
    factor(if_else(x, when_true, when_false), levels = c(when_false, when_true))
}

ext_as_clinical_phenotype <- function(x) {
    factor(x, levels = c(
        "ALS", "ALS/FTD", "FTD",
        "PLS", "PMA", "PBP",
        "UMN-Predominant", "LMN-Predominant",
        "Flail-Arm", "Flail-Leg"
    ))
}
