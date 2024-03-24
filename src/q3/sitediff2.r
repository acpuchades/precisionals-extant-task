library(ggplot2)
library(writexl)

source("src/ext/common.r")

ext_source("src/q3/vc_niv.r")

ext_interactive({
  q3_sitediff2_output_dir <- file.path(q3_output_root_dir, "sitediff2")
  dir.create(q3_sitediff2_output_dir, showWarnings = FALSE, recursive = TRUE)

  sink(file.path(q3_sitediff2_output_dir, "niv-per-site.txt"))
  cat("# NIV STATUS PER SITE (2010-2022)\n\n")
  patients_info.current %$%
    q3_summary_table(site, niv, useNA = "ifany") %>%
    print()
  sink()

  sink(file.path(q3_sitediff2_output_dir, "vc-at-niv-per-site.txt"))
  cat("# VC AT NIV START PER SITE (2010-2022)\n\n")
  patients_info.current.niv %$%
    q3_summary_table(site, vc_at_niv_interval, useNA = "ifany") %>%
    print()
  cat("\n\n")

  cat("# ANOVA: VC AT NIV START PER SITE (2010-2022)\n\n")
  aov(vc_at_niv ~ site, data = patients_info.current.niv) %>%
    summary() %>%
    print()
  sink()

  patients_info.current.niv %>%
    drop_na(vc_at_niv) %>%
    mutate(site = fct_drop(site)) %>%
    ggplot(aes(sample = vc_at_niv)) +
    geom_qq_line() +
    geom_qq() +
    facet_wrap(~site) +
    theme_bw()
  ggsave(file.path(q3_sitediff2_output_dir, "vc-at-niv.qqplot.png"))

  ext_main %>%
    drop_na(age_at_onset, year_of_diagnosis) %>%
    ggplot(aes(year_of_diagnosis, age_at_onset)) +
    geom_jitter(alpha = .3, size = .5, width = .5, height = .5) +
    geom_smooth(aes(color = site), method = "lm") +
    theme_bw() +
    facet_wrap(~site, scales = "free") +
    labs(
      title = "Year of diagnosis vs Age at onset among sites",
      x = "Year of diagnosis", y = "Age at onset"
    ) +
    theme(legend.position = "none")
  ggsave(file.path(q3_sitediff2_output_dir, "age_at_onset-vs-year_of_diagnosis-per-site.png"))

  patients_info.niv %>%
    drop_na(vc_at_niv) %>%
    ggplot(aes(vc_at_niv, fill = site)) +
    geom_density(alpha = .3) +
    scale_fill_custom(drop = FALSE, breaks = str_c("Cohort ", c(1, 3, 4, 5, 7, 8, 9))) +
    labs(title = "Vital Capacity at NIV (Entire Cohort)", x = "Vital capacity (%)", y = "Density", fill = NULL) +
    theme_bw()
  ggsave(file.path(q3_sitediff2_output_dir, "vc-at-niv-overall.density.png"))

  patients_info.niv %>%
    drop_na(vc_at_niv) %>%
    ggplot(aes(x = site, y = vc_at_niv, fill = site)) +
    geom_boxplot() +
    scale_fill_custom(drop = FALSE) +
    labs(title = "Vital Capacity at NIV (Entire Cohort)", x = NULL, y = "VC (%)") +
    theme_bw() +
    theme(legend.position = "none")
  ggsave(file.path(q3_sitediff2_output_dir, "vc-at-niv-overall.boxplot.png"))

  patients_info.current.niv %>%
    drop_na(vc_at_niv) %>%
    ggplot(aes(vc_at_niv, fill = site)) +
    geom_density(alpha = .3) +
    scale_fill_custom(drop = FALSE, breaks = str_c("Cohort ", c(1, 5, 7:9))) +
    labs(title = "Vital Capacity at NIV (2010 – 2022)", x = "Vital capacity (%)", y = "Density", fill = NULL) +
    theme_bw()
  ggsave(file.path(q3_sitediff2_output_dir, "vc-at-niv-after-2010.density.png"))

  patients_info.current.niv %>%
    drop_na(vc_at_niv) %>%
    ggplot(aes(x = site, y = vc_at_niv, fill = site)) +
    geom_boxplot() +
    scale_fill_custom(drop = FALSE) +
    labs(title = "Vital Capacity at NIV (2010 – 2022)", x = NULL, y = "VC (%)") +
    theme_bw() +
    theme(legend.position = "none")
  ggsave(file.path(q3_sitediff2_output_dir, "vc-at-niv-after-2010.boxplot.png"))
})
