suppressPackageStartupMessages({
    library(lubridate)
    library(rlang)
})

source("src/ext/common.r")

# Load source datasets and rename fields
ext_main <- ext_load_data(
    "P-ALS_Ext_V2_Main_Data_File.xlsx",
    col_types = c(
        "text", # ID
        "text", # Site
        "text", # Sex
        "date", # Date of Birth
        "text", # Year/Year and Month of Birth
        "numeric", # Age
        "text", # First Symptom
        "text", # Site of Onset
        "text", # Site of Onset 2
        "text", # Site of Onset 3
        "text", # Site of Onset 4
        "text", # Side of Onset
        "date", # Date of Onset
        "date", # Month of Onset
        "numeric", # Age at Onset
        "numeric", # Calculated Age at Onset
        "text", # Diagnosis
        "text", # Diagnosis 2
        "text", # Diagnosis 3
        "text", # Other Diagnosis
        "text", # Motor Neuron Predominance
        "date", # Date of Diagnosis
        "date", # Month of Diagnosis
        "numeric", # Age at Diagnosis
        "numeric", # Calculated Age at Diagnosis
        "text", # Vital Status
        "date", # Date of Death
        "numeric", # Age at Death
        "numeric", # Calculated Age at Death
        "text", # Tracheostomy
        "date", # Date of Tracheostomy
        "numeric", # Age at Tracheostomy
        "text", # >23h NIV
        "date", # Date of 23h NIV
        "numeric", # Age at >23 h NIV
        "date", # Date of Last Follow Up
        "numeric", # Age at Last Follow-up (if alive)
        "date", # Date of Transfer
        "numeric", # Age at Transfer (if alive)
        "text", # Non-invasive Ventilation
        "date", # Date of Non-invasive Ventialtion
        "numeric", # Age at Non-invasive Ventilation
        "text", # Gastrostomy
        "date", # Date of Gastrostomy
        "numeric", # Age at Gastrostomy
        "text", # C9orf72 Tested
        "text", # C9orf72 Status
        "text", # Commercial Result
        "text", # SOD1 Tested
        "text", # SOD1 Status
        "text", # FUS Tested
        "text", # FUS Status
        "text", # TARDBP Tested
        "text", # TARDBP Status
        "text", # Riluzole Use
        "date", # Riluzole Start Date
        "numeric", # Riluzole Start Age
        "text", # Rilzole Stopped
        "date", # Riluzole Stop Date
        "text", # Edaravone Use
        "text", # Edaravone Stopped
        "date", # Edaravone Stop Date
        "text" # Current Working Status
    )
) %>%
    rename_with(~ str_replace(.x, "non_invasive_venti(la|al)tion$", "niv")) %>%
    rename_with(~ str_replace(.x, "rilzole", "riluzole")) %>%
    rename_with(~ str_replace(.x, "_if_alive$", "")) %>%
    rename_with(~ str_replace(.x, "_gt_23h_niv", "_23h_niv")) %>%
    rows_update(tibble(id = "FRA-0560", date_of_niv = dmy("15/10/2016")), by = "id")

# Apply Leuven corrections
ext_leuven_corrections_path <- "PALS_Leuven_Corrected_Age_Death_or_Transfer_2024-01-18_cleaning.xlsx"
if (file.exists(file.path("data", ext_leuven_corrections_path))) {
    ext_leuven_corrections <- ext_load_data(ext_leuven_corrections_path) %>%
        rename(
            id = "pals_id",
            corrected_age_at_transfer = "precision_age_at_death_transfer_if_still_alive",
        ) %>%
        inner_join(
            ext_main %>% select(id, vital_status, age_at_transfer, age_at_death),
            by = "id"
        ) %>%
        mutate(
            corrected_vital_status = if_else(
                !is.na(corrected_age_at_death) & vital_status != "Deceased",
                "Deceased", NA_character_
            )
        ) %>%
        select(id, starts_with("corrected_"))

    ext_main <- ext_main %>% ext_apply_corrections(ext_leuven_corrections)
}

# Make data fields uniform and add calculated fields
ext_main <- ext_main %>%
    mutate(
        across(ends_with("_tested"), ext_parse_boolean),
        across(c(gastrostomy, niv, tracheostomy), ext_parse_boolean),
        sex = ext_as_sex(case_match(
            sex,
            c("Man", "Male") ~ "Male",
            c("Woman", "Female") ~ "Female",
        )),
        riluzole_use = case_match(
            riluzole_use,
            "Yes" ~ TRUE,
            "No" ~ FALSE
        ),
        c9orf72_status = ext_as_gene_status(case_match(
            c9orf72_status,
            "Positive" ~ "Positive",
            c("Intermediate", "Negative") ~ "Negative"
        )),
        sod1_status = ext_as_gene_status(case_when(
            sod1_status == "Unknown effect" ~ NA_character_,
            TRUE ~ sod1_status
        )),
        fus_status = ext_as_gene_status(fus_status),
        tardbp_status = ext_as_gene_status(tardbp_status),
        date_of_birth = coalesce(
            date_of_birth,
            make_date(str_extract(year_year_and_month_of_birth, "\\d{4}"), 1, 1),
            make_date(
                str_extract(year_year_and_month_of_birth, "(\\d{4})-\\d{1,2}", group = 1),
                str_extract(year_year_and_month_of_birth, "\\d{4}-(\\d{1,2})", group = 1),
                1
            )
        ),
        across(starts_with("date_of_") & -date_of_birth,
            ~ (.x - date_of_birth) / dyears(1),
            .names = "calculated_age_from_{.col}"
        ),
        calculated_age_from_date_of_transfer = if_else(
            vital_status == "Alive", calculated_age_from_date_of_transfer, NA_real_
        ),
        # BUG: for some reason this doesn't work...
        #
        # across(starts_with("age_at_"), \(x) {
        #    event <- str_extract(cur_column(), "age_at_([a-z0-9_]+)", group = 1)
        #    calculated_age_col <- str_glue("calculated_age_from_date_of_{event}")
        #    coalesce(x, .data[[calculated_age_col]])
        # }, .names = "coalesced_{.col}"),
        age_at_onset = coalesce(age_at_onset, calculated_age_at_onset, calculated_age_from_date_of_onset),
        age_at_death = coalesce(age_at_death, calculated_age_at_death, calculated_age_from_date_of_death),
        age_at_diagnosis = coalesce(age_at_diagnosis, calculated_age_at_diagnosis, calculated_age_from_date_of_diagnosis),
        age_at_gastrostomy = coalesce(age_at_gastrostomy, calculated_age_from_date_of_gastrostomy),
        age_at_23h_niv = coalesce(age_at_23h_niv, calculated_age_from_date_of_23h_niv),
        age_at_niv = coalesce(age_at_niv, calculated_age_from_date_of_niv),
        age_at_tracheostomy = coalesce(age_at_tracheostomy, calculated_age_from_date_of_tracheostomy),
        age_at_transfer = if_else(vital_status == "Alive", coalesce(
            age_at_transfer, calculated_age_from_date_of_transfer
        ), NA_real_),
        age_at_last_follow_up = coalesce(
            age_at_last_follow_up,
            calculated_age_from_date_of_last_follow_up
        ),
        year_of_diagnosis = coalesce(
            year(date_of_diagnosis),
            year(date_of_birth + dyears(age_at_diagnosis))
        ),
        bulbar_onset = (site_of_onset %in% c(
            "Bulbar", "Bulbaire",
            "Bulbar and Cognitive/Behavioural",
            "Bulbar and Spinal", "Bulbar and Spinal",
            "Bulbar and Thoracic/Respiratory",
            "Cognitive/Behavioural and Bulbar",
            "PBP"
        )) | coalesce(diagnosis == "PBP", FALSE),
        spinal_onset = site_of_onset %in% c(
            "Arms", "Cervical",
            "Cognitive/Behavioural and Spinal",
            "Flail-Arm", "Flail-Leg",
            "Hemiplegic", "Lower limb",
            "Membre inférieur distal D",
            "Membre inférieur distal G",
            "Membre inférieur distal Bilat",
            "Membre inférieur proximal D",
            "Membre inférieur proximal G",
            "Membre inférieur proximal Bilat",
            "Membre supérieur distal D",
            "Membre supérieur distal G",
            "Membre supérieur distal Bilat",
            "Membre supérieur proximal D",
            "Membre supérieur proximal G",
            "Membre supérieur proximal Bilat",
            "Monomyelic", "Neck",
            "Pseudopolyneuritic",
            "Spinal", "Bulbar and Spinal",
            "Spinal and Cognitive/Behavioural",
            "Thoracic/Respiratory and Spinal",
            "Trunk", "trunk", "Upper limb"
        ),
        respiratory_onset = site_of_onset %in% c(
            "Bulbar and Thoracic/Respiratory",
            "Respiratory", "Respiratoire",
            "respiratory",
            "Thoracic/respiratory",
            "Thoracic/Respiratory",
            "Thoracic/Respiratory and Spinal"
        ),
        cognitive_onset = (site_of_onset %in% c(
            "Bulbar and Cognitive/Behavioural",
            "Cognitive",
            "Cognitive/Behavioural",
            "Cognitive/Behavioural and Bulbar",
            "Cognitive/Behavioural and Spinal",
            "Cognitive impairment",
            "FTD"
        )) | coalesce(diagnosis == "FTD", FALSE),
        generalized_onset = site_of_onset %in% c(
            "Generalized",
            "Generalised",
            "Mixed Presentation"
        ),
        side_of_onset = case_when(
            side_of_onset == "Right" ~ "R",
            side_of_onset == "Left" ~ "L",
            side_of_onset == "Both sides" ~ "B",
            str_ends(site_of_onset, " D") ~ "R",
            str_ends(site_of_onset, " G") ~ "L",
            str_ends(site_of_onset, " Bilat") ~ "B"
        ),
        clinical_phenotype = ext_as_clinical_phenotype(coalesce(
            case_match(
                diagnosis,
                "ALS" ~ case_match(
                    motor_neuron_predominance,
                    "LMN" ~ "LMN-Predominant",
                    "UMN" ~ "UMN-Predominant",
                    .default = "ALS"
                ),
                "ALS/FTD" ~ "ALS/FTD",
                "FTD" ~ "FTD",
                "Flail arm" ~ "Flail-Arm",
                "Flail leg" ~ "Flail-Leg",
                "LMN Predominant ALS" ~ "LMN-Predominant",
                "UMN Predominant ALS" ~ "UMN-Predominant",
                "PBP" ~ "PBP",
                "PMA" ~ "PMA",
                c("PLS", "PLS/ALS") ~ "PLS"
            ),
            case_match(
                site_of_onset,
                "PMA" ~ "PMA",
                "PBP" ~ "PBP",
                "Flail-Arm" ~ "Flail-Arm",
                "Flail-Leg" ~ "Flail-Leg",
                c("PLS", "PLS/ALS", "Suspected PLS") ~ "PLS",
                .default = if_else(site == "Bellvitge", "ALS", NA)
            )
        )),
        diagnosis_period = factor(if_else(!is.na(year_of_diagnosis),
            {
                period_start <- year_of_diagnosis - year_of_diagnosis %% 10
                str_glue("{period_start}-{pmin(period_start+9, 2022, na.rm = TRUE)}")
            },
            NA_character_
        ), ordered = TRUE)
    )

ext_sod1_info <- read_xlsx("data/SOD1 raw data collated.xlsx") %>%
    select(
        patient_id = "pals_id",
        sod1_position = "Pos",
        sod1_zygosity = "Zygosity",
        sod1_variant_p = "Corrected name",
        sod1_significance = "Variant significance",
        sod1_clinvar_name = "Name"
    ) %>%
    mutate(
        across(sod1_position, as.integer),
        across(!where(is.numeric), ~ na_if(.x, "-")),
    )

ext_main <- ext_main %>% left_join(
    ext_sod1_info,
    by = c(id = "patient_id")
)

ext_main.anon <- ext_main %>%
    mutate(site = factor(site, labels = str_c("Cohort ", 1:n_distinct(site)))) %>%
    arrange(site, id)

ext_interactive({
    library(writexl)

    ext_main.anon %>% write_xlsx("output/entire-cohort.xlsx")

    followup_times <- ext_main.anon %>%
        transmute(
            id, site,
            vital_status,
            calculated_age_from_date_of_transfer,
            calculated_age_at_death,
            age_at_followup_start = calculated_age_at_onset,
            age_at_followup_end = coalesce(
                calculated_age_at_death,
                calculated_age_from_date_of_transfer
            ),
            followup_duration = age_at_followup_end - age_at_followup_start
        )

    write_xlsx(list(
        "From onset" = followup_times %>%
            bind_rows(followup_times %>% mutate(site = "Total")) %>%
            summarize(
                total_patient_years = sum(followup_duration, na.rm = TRUE),
                mean = mean(followup_duration, na.rm = TRUE),
                median = median(followup_duration, na.rm = TRUE),
                sd = sd(followup_duration, na.rm = TRUE),
                iqr = IQR(followup_duration, na.rm = TRUE),
                .by = site
            ),
        "Per patient" = followup_times
    ), "output/followup-times.xlsx")
})
