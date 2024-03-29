library("tidyverse")


# Is the age at diagnosis different at each site?
ext_main %>%
  filter(site != "Karolinska") %>% # Removing for now to investigate why
  ggplot(aes(x = calculated_age_at_diagnosis, fill = site)) +
  geom_histogram() +
  facet_wrap(~site)

# Karolinsk is missing this - checked with Éanna, should be removed for age of diag anal

ext_main %>%
  filter(site == "Karolinska")
select(ext_main$calculated_age_at_diagnosis)
unique(ext_main$site)

# I almost prefer it without the facet wrap actually, clearer image that not that different
ext_main %>%
  ggplot(aes(x = calculated_age_at_diagnosis, fill = site)) +
  geom_histogram()

# With c9, is it different?
ext_main %>%
  filter(site != "Karolinsk") %>%
  filter(., c9orf72_stat_BOOL == TRUE) %>%
  ggplot(aes(x = calculated_age_at_diagnosis, fill = site)) +
  geom_histogram() +
  facet_wrap(~site)

# Violins work better
ext_main %>%
  filter(site != "Karolinska") %>%
  ggplot(aes(x = site, y = calculated_age_at_diagnosis, fill = site)) +
  geom_boxplot()


# Age of onset
#karolinksa has date of onset and dob - needs calc

karo <- ext_main %>%
  filter(site == "Karolinska") %>% 
  select(id, date_of_birth, date_of_onset) %>% 
  mutate(karo_age_onset = as.numeric(interval(date_of_birth, date_of_onset) / dyears(1))
         )
ext_main <- ext_main %>% 
  left_join(karo, by = c("id", "date_of_birth", "date_of_onset")) %>%
  mutate(new_age_of_onset = coalesce(calculated_age_at_onset, karo_age_onset))

ext_main %>%
  ggplot(aes(x = new_age_of_onset, fill = site)) +
  geom_histogram() +
  facet_wrap(~site)

ext_main %>%
  filter(date_of_diagnosis) %>%
  ggplot(aes(x = site, y = calculated_age_at_onset, fill = site)) +
  geom_violin()


library(nortest)
ad.test(ext_main$new_age_of_onset)
#non-normal
kruskal.test(new_age_of_onset ~ site, data = ext_main)

library(FSA)
dunnTest(new_age_of_onset ~ site, 
          data = ext_main,
          method = "bonferroni")

#dif with genes?

headline_stats <- ext_main %>%
  summarise(n = n(),
            spinal = mean(spinal_onset)*100,
            bulbar = mean(bulbar_onset)*100,
            respiratory = mean(respiratory_onset)*100,
            cognitive = mean(cognitive_onset)*100,
            median_age_of_onset = median(new_age_of_onset, na.rm = TRUE),
            IQR_onset = IQR(new_age_of_onset, na.rm = TRUE),
            median_age_of_death = median(age_at_death, na.rm = TRUE),
            IQR_death = IQR(age_at_death, na.rm = TRUE))

overall_tested_onset <- ext_main %>% 
  filter(c9orf72_test_BOOL == T | sod1_test_BOOL ==T | fus_tested_BOOL == T | tardbp_tested_BOOL == T) %>% 
  summarise(n = n(),
            spinal = mean(spinal_onset)*100,
            bulbar = mean(bulbar_onset)*100,
            respiratory = mean(respiratory_onset)*100,
            cognitive = mean(cognitive_onset)*100,
            median_age_of_onset = median(new_age_of_onset, na.rm = TRUE),
            IQR_onset = IQR(new_age_of_onset, na.rm = TRUE),
            median_age_of_death = median(age_at_death, na.rm = TRUE),
            IQR_death = IQR(age_at_death, na.rm = TRUE))

fus_onset <- ext_main %>% 
  filter(fus_tested_BOOL == T & fus_status_BOOL == F) %>% 
  summarise(n = n(),
            spinal = mean(spinal_onset)*100,
            bulbar = mean(bulbar_onset)*100,
            respiratory = mean(respiratory_onset)*100,
            cognitive = mean(cognitive_onset)*100,
            median_age_of_onset = median(new_age_of_onset, na.rm = TRUE),
            IQR_onset = IQR(new_age_of_onset, na.rm = TRUE),
            median_age_of_death = median(age_at_death, na.rm = TRUE),
            IQR_death = IQR(age_at_death, na.rm = TRUE))
tardbp_onset <- ext_main %>% 
  filter(tardbp_tested_BOOL == T & tardbp_status_BOOL == F) %>% 
  summarise(n = n(),
            spinal = mean(spinal_onset)*100,
            bulbar = mean(bulbar_onset)*100,
            respiratory = mean(respiratory_onset)*100,
            cognitive = mean(cognitive_onset)*100,
            median_age_of_onset = median(new_age_of_onset, na.rm = TRUE),
            IQR_onset = IQR(new_age_of_onset, na.rm = TRUE),
            median_age_of_death = median(age_at_death, na.rm = TRUE),
            IQR_death = IQR(age_at_death, na.rm = TRUE))

sod1_onset <- ext_main %>% 
  filter(sod1_test_BOOL == T & sod1_stat_BOOL == F) %>% 
  summarise(n = n(),
            spinal = mean(spinal_onset)*100,
            bulbar = mean(bulbar_onset)*100,
            respiratory = mean(respiratory_onset)*100,
            cognitive = mean(cognitive_onset)*100,
            median_age_of_onset = median(new_age_of_onset, na.rm = TRUE),
            IQR_onset = IQR(new_age_of_onset, na.rm = TRUE),
            median_age_of_death = median(age_at_death, na.rm = TRUE),
            IQR_death = IQR(age_at_death, na.rm = TRUE))

tested <- ext_main %>% 
  filter(c9orf72_test_BOOL == T | sod1_test_BOOL ==T | fus_tested_BOOL == T | tardbp_tested_BOOL == T)
wilcox.test(new_age_of_onset ~ c9orf72_stat_BOOL, data = tested) #yes
wilcox.test(new_age_of_onset ~ sod1_stat_BOOL, data = tested) #yes
wilcox.test(new_age_of_onset ~ fus_status_BOOL, data = tested) #yes
wilcox.test(new_age_of_onset ~ tardbp_status_BOOL, data = tested) #yes



median(ext_main$new_age_of_onset, na.rm = TRUE)
IQR(ext_main$new_age_of_onset, na.rm = TRUE)

