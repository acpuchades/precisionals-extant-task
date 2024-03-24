library(tibble)
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

q3_compare_cohort_numeric_fields <- function(x, ...) {
  UseMethod("q3_compare_cohort_numeric_fields", x)
}

q3_compare_cohort_numeric_fields.data.frame <- function(target, ref, ref.group = "ref", target.group = "target", exclude = NULL) {
  bind_rows(
    target %>% select(where(is.numeric), -all_of(exclude)) %>% mutate(group = target.group, .before = everything()),
    ref %>% select(where(is.numeric), -all_of(exclude)) %>% mutate(group = ref.group, .before = everything())
  ) %>%
    pivot_longer(-group) %>%
    reframe(
      lm(value ~ group) %>%
        tidy(conf.int = TRUE) %>%
        filter(term == str_c("group", target.group)) %>%
        select(-term),
      .by = name
    ) %>%
    mutate(target = target.group, ref = ref.group, .before = everything())
}

q3_compare_cohort_numeric_fields.mids <- function(target.mids, ref, target.group = "target", ref.group = "ref", exclude = NULL) {
  stats <- NULL
  num.cols <- colnames(target.mids$data %>% select(where(is.numeric)))
  groups <- cbind(target.mids, group = rep(target.group, nrow(complete(target.mids))))
  groups <- rbind(groups, mutate(ref, group = ref.group))
  for (col in num.cols) {
    results <- with(groups, lm(reformulate("group", col)))
    stats <- bind_rows(stats, bind_cols(
      summary(pool(results), conf.int = TRUE) %>%
        filter(term == paste0("group", target.group)) %>%
        select(-term, -df),
      name = col
    ))
  }

  mutate(stats, target = target.group, ref = ref.group, .before = everything())
}

q3_compare_cohort_categorical_fields <- function(x, ...) {
  UseMethod("q3_compare_cohort_categorical_fields", x)
}

q3_compare_cohort_categorical_fields.data.frame <- function(target, ref, target.group = "target", ref.group = "ref", exclude = NULL) {
  bind_rows(
    target %>%
      select(-where(is.numeric), -all_of(exclude)) %>%
      mutate(group = target.group, .before = everything()),
    ref %>%
      select(-where(is.numeric), -all_of(exclude)) %>%
      mutate(group = ref.group, .before = everything())
  ) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(-group) %>%
    reframe(
      chisq.test(table(group, value)) %>% tidy(),
      .by = name
    ) %>%
    mutate(target = target.group, ref = ref.group, .before = everything()) %>%
    relocate(method, .after = name)
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
    "Numeric Tests" = bind_rows(
      q3_compare_cohort_numeric_fields(q3_cohort.target, q3_cohort.ref, target.group = "target", ref.group = "ref", exclude = "diagnosis_period"),
      q3_compare_cohort_numeric_fields(q3_cohort.imputed, q3_cohort.target, target.group = "target.imputed", ref.group = "target", exclude = c("riluzole_use", "diagnosis_period"))
    ) %>% arrange(name, target),
    "Categorical" = bind_rows(
      q3_calculate_cohort_categorical_stats(q3_cohort.ref) %>% mutate(cohort = "ref", .before = everything()),
      q3_calculate_cohort_categorical_stats(q3_cohort.target) %>% mutate(cohort = "target", .before = everything()),
      q3_calculate_cohort_categorical_stats(q3_cohort.imputed) %>% mutate(cohort = "target.imputed", .before = everything())
    ) %>% arrange(name, value, cohort),
    "Categorical Tests" = bind_rows(
      q3_compare_cohort_categorical_fields(q3_cohort.target, q3_cohort.ref, target.group = "target", ref.group = "ref", exclude = c("id", "diagnosis_period")),
    )
  ), file.path(q3_cohorts_root_dir, "cohorts-stats.xlsx"))
})
