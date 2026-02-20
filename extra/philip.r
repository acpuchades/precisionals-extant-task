library(ggsurvfit)

source("src/q3/timetoevent.r")

cohort.2010 <- q3_data %>% filter(diagnosis_period == "2010-2019")

onset_death.2010 <- cohort.2010 %>%
  filter(origin == "onset", event == "death") %>%
  mutate(
    status = if_else((status == "event") & (duration <= dyears(10)), 1, 0),
    time = pmin(duration, dyears(10)) / dyears(1)
  )

onset_niv.2010 <- cohort.2010 %>%
  filter(origin == "onset", event == "niv") %>%
  mutate(
    status = if_else((status == "event") & (duration <= dyears(10)), 1, 0),
    time = pmin(duration, dyears(10)) / dyears(1)
  )

onset_death.2010.km <- survfit2(Surv(time, status) ~ 1, onset_death.2010)
onset_niv.2010.km <- survfit2(Surv(time, status) ~ 1, onset_niv.2010)

ggsurvfit(onset_death.2010.km) + ggtitle("Death") + add_risktable() + scale_x_continuous(breaks = 0:10)
ggsave("ppc-survival-onset-death-2010_2019.png")

ggsurvfit(onset_niv.2010.km) + ggtitle("NIV") + add_risktable() + scale_x_continuous(breaks = 0:10)
ggsave("ppc-survival-onset-niv-2010_2019.png")

sink("ppc-survival-onset-niv-2010_2019.txt")
print(summary(onset_niv.2010.km))
sink()

vc_tours <- ext_resp %>%
  filter(id %>% str_starts("FRA")) %>%
  select(id, site, fvc_abs, fvc_rel, svc_abs, svc_rel, vc_rel, vc_abs)
