library(forcats)
library(ggplot2)
library(stringr)
library(tibble)

source("src/ext/resp.r")
source("src/ext/alsfrs.r")

q3_summary_table <- function(...) {
    t <- table(...)
    mt <- margin.table(t, margin = 1)
    pt <- round(prop.table(t, margin = 1) * 100, 2)
    colnames(pt) <- str_c(colnames(pt), " (%)")
    cbind(n = mt, t, pt)
}

patients_vstatus <- ext_main %>%
    transmute(site, id, is_alive = vital_status == "Alive")

niv_status_as_reported <- ext_main %>%
    select(site, id, date_of_assessment = "date_of_niv", age_at_assessment = "age_at_niv") %>%
    mutate(niv_carrier = !is.na(date_of_assessment) | !is.na(age_at_assessment)) %>%
    left_join(ext_baseline %>% select(id, date_of_baseline, age_at_baseline), by = "id") %>%
    mutate(.after = id, time_from_baseline = coalesce(
        date_of_assessment - date_of_baseline,
        dyears(age_at_assessment - age_at_baseline),
        ddays(0)
    )) %>%
    select(-date_of_baseline, -age_at_baseline)

niv_status_by_alsfrs <- ext_alsfrs_followups %>%
    select(site, id, time_from_baseline, date_of_assessment, age_at_assessment, q12_respiratory_insufficiency) %>%
    mutate(niv_carrier = q12_respiratory_insufficiency %>% between(1, 3), .keep = "unused") %>%
    arrange(id, time_from_baseline)

vc_assessments_from_baseline <- ext_resp %>%
    select(site, id, date_of_assessment, age_at_assessment, vc_rel) %>%
    left_join(
        ext_baseline %>% select(id, date_of_baseline, age_at_baseline),
        by = "id"
    ) %>%
    mutate(
        .after = id, time_from_baseline = coalesce(
            date_of_assessment - date_of_baseline,
            dyears(age_at_assessment - age_at_baseline)
        )
    ) %>%
    filter(time_from_baseline >= ddays(0)) %>%
    select(-date_of_baseline, -age_at_baseline) %>%
    drop_na(vc_rel)

vc_and_niv_status_from_baseline <- niv_status_as_reported %>%
    bind_rows(niv_status_by_alsfrs) %>%
    semi_join(vc_assessments_from_baseline, by = "id") %>%
    bind_rows(vc_assessments_from_baseline) %>%
    group_by(id) %>%
    arrange(time_from_baseline, .by_group = TRUE) %>%
    fill(niv_carrier, vc_rel) %>%
    ungroup()

vc_at_niv_start <- vc_and_niv_status_from_baseline %>%
    filter(niv_carrier == TRUE) %>%
    slice_min(time_from_baseline, by = id, n = 1, with_ties = FALSE) %>%
    select(site, id, vc_at_niv = "vc_rel")

vc_niv_summary <- patients_vstatus %>%
    inner_join(
        vc_and_niv_status_from_baseline %>% select(id, vc_rel),
        by = "id"
    ) %>%
    summarize(
        vc_min = min(vc_rel, na.rm = TRUE),
        vc_max = max(vc_rel, na.rm = TRUE),
        .by = c("site", "id")
    ) %>%
    left_join(vc_at_niv_start %>% select(id, vc_at_niv), by = "id") %>%
    mutate(vc_at_niv_interval = factor(case_when(
        vc_at_niv < 50 ~ "<50",
        vc_at_niv %>% between(50, 60) ~ "51-60",
        vc_at_niv %>% between(60, 70) ~ "61-70",
        vc_at_niv > 70 ~ ">70"
    ), ordered = TRUE, levels = c("<50", "51-60", "61-70", ">70")))

sink("output/q3/niv/niv-per-site.table.txt")
cat("# NIV STATUS PER SITE AS REPORTED (N/A removed)\n\n")
ext_main %$%
    q3_summary_table(site, niv) %>%
    print()
cat("\n\n")
cat("# NIV STATUS PER SITE AS REPORTED (N/A included)\n\n")
ext_main %$%
    q3_summary_table(site, niv, useNA = "ifany") %>%
    print()
sink()

sink("output/q3/niv/niv-per-period.table.txt")
cat("# NIV STATUS PER PERIOD AS REPORTED (N/A removed)\n\n")
ext_main %$%
    q3_summary_table(diagnosis_period, niv) %>%
    print()
cat("\n\n")
cat("# NIV STATUS PER PERIOD AS REPORTED (N/A included)\n\n")
ext_main %$%
    q3_summary_table(diagnosis_period, niv, useNA = "ifany") %>%
    print()
sink()

sink("output/q3/niv/vc-at-niv.table.txt")
vc_niv_summary %$%
    q3_summary_table(site, vc_at_niv_interval) %>%
    print()
sink()

ggplot(vc_niv_summary, aes(sample = vc_at_niv)) +
    geom_qq_line() +
    geom_qq() +
    facet_wrap(~site) +
    theme_bw()
ggsave("output/q3/niv/vc-at-niv.qqplot.png")

sink("output/q3/niv/vc-at-niv.aov.txt")
aov(vc_at_niv ~ site, vc_niv_summary) %>%
    summary() %>%
    print()
sink()

ggplot(vc_niv_summary, aes(vc_at_niv, fill = site)) +
    geom_density(alpha = 0.3) +
    labs(x = "Vital capacity (%)", y = "Density", fill = "Site") +
    theme_bw() +
    theme(legend.position = "none")
ggsave("output/q3/niv/vc-at-niv-per-site.density.png")

ggplot(vc_niv_summary, aes(x = site, y = vc_at_niv, fill = site)) +
    geom_boxplot() +
    labs(title = "VC at NIV", x = "Site", y = "VC (%)") +
    theme_bw() +
    theme(legend.position = "none")
ggsave("output/q3/niv/vc-at-niv-per-site.boxplot.png")
