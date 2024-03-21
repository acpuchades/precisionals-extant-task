library(writexl)

source("src/ext/common.r")

ext_source("src/q3/helpers.r")
ext_source("src/q3/impute2.r")

q3_calculate_cohort_numeric_stats <- function(x, ...) {
  UseMethod("q3_calculate_cohort_numeric_stats", x)
}

q3_calculate_cohort_categorical_stats <- function(x, ...) {
  UseMethod("q3_calculate_cohort_categorical_stats", x)
}

q3_calculate_cohort_numeric_stats.data.frame <- function(df) {
  df %>%
    select(where(is.numeric)) %>%
    pivot_longer(everything()) %>%
    summarize(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n()),
      .by = name
    ) %>%
    arrange(name)
}

q3_calculate_cohort_categorical_stats.data.frame <- function(df) {
  df %>%
    select(where(is.factor)) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(everything()) %>%
    summarize(n = n(), .by = c(name, value)) %>%
    mutate(p = n / sum(n), .by = name) %>%
    arrange(name, value)
}

# REF: https://github.com/rwnahhas/RMPH_Resources/blob/main/Functions_rmph.R
q3_calculate_cohort_numeric_stats.mids <- function(mids, ...) {
  stats <- NULL
  impdat <- complete(mids, action = "long")
  num.cols <- colnames(mids$data %>% select(where(is.numeric)))
  for (col in num.cols) {
    col_mean <- tapply(impdat[[col]], impdat$.imp, mean, na.rm = TRUE)
    col_varmean <- tapply(impdat[[col]], impdat$.imp, var.mean, na.rm = TRUE)
    col_sd <- tapply(impdat[[col]], impdat$.imp, sd, na.rm = TRUE)
    pooled <- pool.scalar(
      Q = col_mean,
      U = col_varmean,
      n = nrow(mids$data)
    )

    stats <- bind_rows(stats, tibble(
      name = col, mean = pooled$qbar, se = sqrt(pooled$t), sd = mean(col_sd)
    ))
  }
  stats
}

# REF: https://github.com/rwnahhas/RMPH_Resources/blob/main/Functions_rmph.R
q3_calculate_cohort_categorical_stats.mids <- function(mids, ...) {
  stats <- NULL
  cat.cols <- colnames(mids$data %>% select(where(is.factor)))
  for (col in cat.cols) {
    pmat <- with(mids, table(get(col)))$analyses %>%
      unlist() %>%
      matrix(nrow = mids$m, byrow = TRUE)
    n <- apply(pmat, MARGIN = 2, FUN = mean, na.rm = TRUE)
    p <- n / sum(n)
    stats <- bind_rows(stats, tibble(
      name = col, value = levels(mids$data[[col]]), n = n, p = p
    ))
  }
  stats
}

q3_cohort.ref <- q3_base %>%
  q3_add_derived_variables()

q3_cohort.target <- q3_base_recent %>%
  q3_add_derived_variables()

q3_cohort.imputed <- q3_data_recent_w.mids %>%
  complete(action = "long", include = TRUE) %>%
  select(all_of(colnames(q3_base)), .imp) %>%
  q3_add_derived_variables() %>%
  as.mids()

ext_interactive({
  q3_cohorts_root_dir <- file.path(q3_output_root_dir, "cohorts")
  dir.create(q3_cohorts_root_dir, showWarnings = FALSE, recursive = TRUE)

  write_xlsx(list(
    "Numeric" = bind_rows(
      q3_calculate_cohort_numeric_stats(q3_cohort.ref) %>% mutate(cohort = "ref", .before = everything()),
      q3_calculate_cohort_numeric_stats(q3_cohort.target) %>% mutate(cohort = "target", .before = everything()),
      q3_calculate_cohort_numeric_stats(q3_cohort.imputed) %>% mutate(cohort = "target.imputed", .before = everything())
    ) %>% arrange(name, cohort),
    "Categorical" = bind_rows(
      q3_calculate_cohort_categorical_stats(q3_cohort.ref) %>% mutate(cohort = "ref", .before = everything()),
      q3_calculate_cohort_categorical_stats(q3_cohort.target) %>% mutate(cohort = "target", .before = everything()),
      q3_calculate_cohort_categorical_stats(q3_cohort.imputed) %>% mutate(cohort = "target.imputed", .before = everything())
    ) %>% arrange(name, value, cohort)
  ), file.path(q3_cohorts_root_dir, "cohorts-stats.xlsx"))
})
