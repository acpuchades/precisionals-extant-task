library(dplyr)

source("src/ext/main.r")
source("src/ext/alsfrs.r")
source("src/ext/resp.r")
source("src/q3/timetoevent.r")

time_of_niv_as_reported <- ext_main.anon %>%
  select(id, niv, date_of_niv, age_at_niv)

niv_status_as_reported <- ext_main.anon %>%
  select(id, date_of_assessment = "date_of_niv", age_at_assessment = "age_at_niv", niv) %>%
  left_join(ext_baseline %>% select(id, date_of_baseline, age_at_baseline), by = "id") %>%
  mutate(.after = id, time_from_baseline = coalesce(
    as.duration(date_of_assessment - date_of_baseline),
    dyears(age_at_assessment - age_at_baseline)
  ))

niv_status_by_alsfrs <- ext_alsfrs_followups %>%
  mutate(niv = q12_respiratory_insufficiency %>% between(1, 3)) %>%
  select(id, time_from_baseline, date_of_assessment, age_at_assessment, niv)

time_of_niv_by_alsfrs <- niv_status_by_alsfrs %>%
  summarize(niv = any(niv, na.rm = TRUE), .by = id) %>%
  left_join(
    niv_status_by_alsfrs %>%
      filter(niv == TRUE) %>%
      slice_min(time_from_baseline, by = id, n = 1, with_ties = FALSE) %>%
      select(id, date_of_niv = "date_of_assessment", age_at_niv = "age_at_assessment"),
    by = "id"
  )

niv_assessments <- niv_status_as_reported %>%
  bind_rows(niv_status_by_alsfrs)

time_of_niv <- time_of_niv_as_reported %>%
  bind_rows(time_of_niv_by_alsfrs) %>%
  summarize(
    .by = id,
    niv = any(niv, na.rm = TRUE),
    age_at_niv = if_else(sum(!is.na(age_at_niv)) != 0, suppressWarnings(min(age_at_niv, na.rm = TRUE)), NA),
    date_of_niv = if_else(sum(!is.na(date_of_niv)) != 0, suppressWarnings(min(date_of_niv, na.rm = TRUE)), NA)
  )

vc_assessments <- ext_resp %>%
  select(id, date_of_assessment, age_at_assessment, vc_rel) %>%
  drop_na(vc_rel) %>%
  left_join(
    ext_baseline %>% select(id, date_of_baseline, age_at_baseline),
    by = "id"
  ) %>%
  mutate(
    .after = id, time_from_baseline = coalesce(
      as.duration(date_of_assessment - date_of_baseline),
      dyears(age_at_assessment - age_at_baseline)
    )
  ) %>%
  filter(time_from_baseline >= ddays(0)) %>%
  select(-date_of_baseline, -age_at_baseline)

vc_and_niv_assessments <- vc_assessments %>%
  mutate(time_of_last_vc_assessment = time_from_baseline) %>%
  bind_rows(
    niv_assessments %>%
      semi_join(vc_assessments, by = "id") %>%
      mutate(time_of_last_niv_assessment = time_from_baseline)
  ) %>%
  group_by(id) %>%
  arrange(time_from_baseline, .by_group = TRUE) %>%
  fill(vc_rel, time_of_last_vc_assessment, niv, time_of_last_niv_assessment) %>%
  ungroup() %>%
  mutate(
    vc_rel = if_else(
      (time_from_baseline - time_of_last_vc_assessment) <= dmonths(6),
      vc_rel, NA
    )
  )

vc_at_niv_start <- vc_and_niv_assessments %>%
  filter(niv == TRUE) %>%
  slice_min(time_from_baseline, by = id, n = 1, with_ties = FALSE) %>%
  filter(time_from_baseline > ddays(0)) %>%
  select(id, vc_at_niv = "vc_rel") %>%
  mutate(vc_at_niv_interval = factor(case_when(
    vc_at_niv < 50 ~ "<50",
    vc_at_niv %>% between(50, 60) ~ "51-60",
    vc_at_niv %>% between(60, 70) ~ "61-70",
    vc_at_niv > 70 ~ ">70"
  ), ordered = TRUE, levels = c("<50", "51-60", "61-70", ">70")))

vc_summary <- vc_assessments %>% summarize(
  vc_min = min(vc_rel, na.rm = TRUE),
  vc_max = max(vc_rel, na.rm = TRUE),
  .by = id
)

patients_info <- q3_base %>%
  q3_add_derived_variables() %>%
  transmute(id, site, diagnosis_period) %>%
  left_join(ext_main %>% select(id, vital_status), by = "id") %>%
  left_join(time_of_niv, by = "id") %>%
  left_join(vc_summary, by = "id") %>%
  left_join(vc_at_niv_start, by = "id")

patients_info.vc_at_niv <- patients_info %>% filter(niv == TRUE)
