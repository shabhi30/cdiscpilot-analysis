# CDISC pilot alzheimer's trial — complete biostatistical analysis
#loading packages
install.packages(c("haven","tidyverse","gtsummary","survival","broom","mmrm","flextable"))
install.packages("webshot2")
webshot2::install_chromium()
library(haven)
library(tidyverse)
library(gtsummary)
library(survival)
library(broom)
library(mmrm)
library(flextable)
#setup
dir.create("data",    showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)
dir.create("programs",showWarnings = FALSE)
setwd("~/cdiscpilot-analysis")
# publicly available CDISC ADaM pilot datasets — real FDA submission data .xpt = SAS transport format, 
base_url <- paste0(
  "https://github.com/cdisc-org/sdtm-adam-pilot-project/raw/master/",
  "updated-pilot-submission-package/900172/m5/datasets/",
  "cdiscpilot01/analysis/adam/datasets/"
)
download.file(paste0(base_url, "adsl.xpt"),     "data/adsl.xpt",     mode = "wb")
download.file(paste0(base_url, "adtte.xpt"),    "data/adtte.xpt",    mode = "wb")
download.file(paste0(base_url, "adqsadas.xpt"), "data/adqsadas.xpt", mode = "wb")
download.file(paste0(base_url, "adae.xpt"),     "data/adae.xpt",     mode = "wb")
#read data
# ADSL  — one row per subject, demographics and baseline values
# ADTTE — one row per subject per event, time-to-event data
# ADQS  — one row per subject per visit per scale item, ADAS-Cog scores
# ADAE  — one row per adverse event per subject
adsl  <- read_xpt("data/adsl.xpt")
adtte <- read_xpt("data/adtte.xpt")
adqs  <- read_xpt("data/adqsadas.xpt")
adae  <- read_xpt("data/adae.xpt")
cat("ADSL:", nrow(adsl),  "subjects\n")
cat("ADTTE:", nrow(adtte), "rows\n")
cat("ADQS:", nrow(adqs),  "rows\n")
cat("ADAE:", nrow(adae),  "rows\n")
#baseline demographics 
# ITT population = every randomized patient — required by ICH E9 guidelines as it can confirm arms were balanced at baseline after randomization

adsl_itt <- adsl %>%
  filter(ITTFL == "Y")

table1 <- adsl_itt %>%
  select(TRT01P, AGE, SEX, RACE, BMIBL, MMSETOT, DURDIS, EDUCLVL) %>%
  tbl_summary(
    by = TRT01P,
    label = list(
      AGE     ~ "Age (years)",
      SEX     ~ "Sex",
      RACE    ~ "Race",
      BMIBL   ~ "BMI at baseline (kg/m²)",
      MMSETOT ~ "MMSE total score at baseline",
      DURDIS  ~ "Duration of disease (months)",
      EDUCLVL ~ "Years of education"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_labels() %>%
  modify_caption("**Table 1. Baseline Demographics and Clinical Characteristics (ITT Population)**")

table1

table1 %>%
  as_flex_table() %>%
  flextable::save_as_image(path = "output/Table1_Demographics.png", zoom = 2)

cat("Table 1 saved\n")
#figure 1 - Kaplan Meier model for time to dermatologic event - dermatologic events are the primary safety signal of interest
# CNSR flipped to EVENT because: Raw data: 0 = event happened, 1 = censored
# R survival functions expect: 1 = event happened, 0 = censored
#Placebo set as reference level so legend and model show Drug vs Placebo

adtte_km <- adtte %>%
  filter(PARAM == "Time to First Dermatologic Event", SAFFL == "Y") %>%
  mutate(
    EVENT = 1 - CNSR,
    TRTP  = factor(TRTP, levels = c(
      "Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"
    ))
  )

# Fit Kaplan-Meier — separate curve per treatment arm
km_fit <- survfit(Surv(AVAL, EVENT) ~ TRTP, data = adtte_km)

# Log rank test — are the three curves statistically different or not
logrank_test <- survdiff(Surv(AVAL, EVENT) ~ TRTP, data = adtte_km)
logrank_test

# Convert model to tidy dataframe for ggplot
km_data <- broom::tidy(km_fit)
km_data$strata <- gsub("TRTP=", "", km_data$strata)
km_data$strata <- factor(km_data$strata, levels = c(
  "Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"
))

km_plot <- ggplot(km_data, aes(x = time, y = estimate, color = strata)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = strata),
    alpha = 0.15, color = NA
  ) +
  scale_color_manual(values = c(
    "Placebo"              = "#2E75B6",
    "Xanomeline Low Dose"  = "#70AD47",
    "Xanomeline High Dose" = "#FF0000"
  )) +
  scale_fill_manual(values = c(
    "Placebo"              = "#2E75B6",
    "Xanomeline Low Dose"  = "#70AD47",
    "Xanomeline High Dose" = "#FF0000"
  )) +
  labs(
    x        = "Time (Days)",
    y        = "Probability of Remaining Event-Free",
    title    = "Time to First Dermatologic Event by Treatment Arm",
    subtitle = "Safety Population  |  Log-rank p = 8 x 10^-14",
    color    = "Treatment",
    fill     = "Treatment"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    plot.subtitle   = element_text(size = 11, color = "darkred"),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold")
  )

km_plot

ggsave(
  filename = "output/Figure1_KaplanMeier.png",
  plot = km_plot, width = 10, height = 7, dpi = 300
)

cat("Figure 1 saved\n")

# Table 2 for MMRM primary efficacy analysis
# Endpoint: ADAS-Cog(11) change from baseline — range 0-70, higher = worse. MMRM used because measurements from the same patient across visits are correlated — regular regression would ignore this, inflate sample size, and produce falsely small p-values (serious regulatory problem)

# Model specification per ICH E9 and FDA guidance:
# CHG= change from baseline (outcome — not absolute score)
# TRT01P = treatment main effect
# AVISIT = visit main effect (captures natural disease progression)
# TRT01P:AVISIT  = interaction — treatment effect at each visit separately
# BASE = baseline score covariate (removes starting point noise)
#  us() = unstructured covariance (FDA recommended — no assumptions)
# Data cleaning decisions:
# DTYPE == ""  - keeps only observed records, drops LOCF imputed rows
# ANL01FL == "Y" resolves duplicates using pre-specified analysis flag
# TRT01P from adsl — planned treatment for ITT, not actual (TRTP)
# Baseline excluded from outcome rows — it is the BASE covariate instead
adqs_mmrm <- adqs %>%
  filter(PARAM == "Adas-Cog(11) Subscore") %>%
  filter(ITTFL == "Y") %>%
  filter(AVISIT %in% c("Week 8", "Week 16", "Week 24")) %>%
  filter(DTYPE == "") %>%
  filter(ANL01FL == "Y") %>%
  left_join(adsl %>% select(USUBJID, TRT01P), by = "USUBJID") %>%
  mutate(
    TRT01P  = factor(TRT01P, levels = c(
      "Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"
    )),
    AVISIT  = factor(AVISIT, levels = c("Week 8", "Week 16", "Week 24")),
    USUBJID = factor(USUBJID)
  )

mmrm_fit <- mmrm(
  formula = CHG ~ TRT01P + AVISIT + TRT01P:AVISIT + BASE +
    us(AVISIT | USUBJID),
  data = adqs_mmrm
)

summary(mmrm_fit)
# Extract results — key numbers are TRT01P:AVISIT interaction at Week 24
# Negative estimate = drug arm declined less than placebo = efficacy signal
broom::tidy(mmrm_fit) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  rename(
    Term = term, Estimate = estimate, Std_Error = std.error,
    Statistic = statistic, P_value = p.value, DF = df
  ) %>%
  mutate(Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01  ~ "**",
    P_value < 0.05  ~ "*",
    P_value < 0.10  ~ ".",
    TRUE            ~ "ns"
  )) %>%
  flextable() %>%
  autofit() %>%
  save_as_image(path = "output/Table2_MMRM_Results.png", zoom = 2)

cat("Table 2 saved\n")

# safety profile - AE
# Safety population — patients who received at least one dose (SAFFL)
# TRTEMFL = "Y" filters to treatment-emergent AEs only
# (AEs that started after first dose — pre-existing conditions excluded)
# Denominators from ADSL not ADAE — patients with zero AEs don't appear
# in ADAE so using ADAE as denominator would inflate percentages
# Drug-related = POSSIBLE or PROBABLE per ICH E2A guideline
# "RELATED" does not exist in this dataset — verified with count(AEREL)
# AOCCPFL = "Y" counts each patient once per preferred term
# Standard — safety tables count patients with the event, not event episodes

adae_te <- adae %>% filter(TRTEMFL == "Y")

# Denominators — total Safety population per arm
n_patients <- adsl %>%
  filter(SAFFL == "Y") %>%
  count(TRT01A, name = "N") %>%
  rename(TRTA = TRT01A)

# table 3- overall AE summary

overall_summary <- adae_te %>%
  group_by(TRTA) %>%
  summarise(
    any_ae      = n_distinct(USUBJID),
    any_sae     = n_distinct(USUBJID[AESER == "Y"]),
    any_related = n_distinct(USUBJID[AEREL %in% c("POSSIBLE", "PROBABLE")]),
    any_severe  = n_distinct(USUBJID[AESEV == "SEVERE"]),
    any_death   = n_distinct(USUBJID[AESDTH == "Y"]),
    any_hosp    = n_distinct(USUBJID[AESHOSP == "Y"]),
    .groups = "drop"
  ) %>%
  left_join(n_patients, by = "TRTA") %>%
  mutate(
    fmt_any_ae  = paste0(any_ae,      " (", round(any_ae      / N * 100, 1), "%)"),
    fmt_sae     = paste0(any_sae,     " (", round(any_sae     / N * 100, 1), "%)"),
    fmt_related = paste0(any_related, " (", round(any_related / N * 100, 1), "%)"),
    fmt_severe  = paste0(any_severe,  " (", round(any_severe  / N * 100, 1), "%)"),
    fmt_death   = paste0(any_death,   " (", round(any_death   / N * 100, 1), "%)"),
    fmt_hosp    = paste0(any_hosp,    " (", round(any_hosp    / N * 100, 1), "%)")
  )

overall_summary %>%
  select(
    `Treatment Arm`           = TRTA,
    `N`                       = N,
    `Any TEAE`                = fmt_any_ae,
    `Any serious AE`          = fmt_sae,
    `Any drug-related AE`     = fmt_related,
    `Any severe AE`           = fmt_severe,
    `Any AE leading to death` = fmt_death,
    `Any AE leading to hosp.` = fmt_hosp
  ) %>%
  flextable() %>%
  set_caption("Table 3. Overall Summary of Treatment-Emergent Adverse Events (Safety Population)") %>%
  bold(part = "header") %>%
  bg(bg = "#D6E4F0", part = "header") %>%
  autofit() %>%
   save_as_image(path = "output/Table3_AE_Overall_Summary.png", zoom = 2)

cat("Table 3 saved\n")

# table 4: AEs by Body System and Preferred Term (>= 5% in any arm)

# build wide table — one column per treatment arm
ae_by_system <- adae_te %>%
  filter(AOCCPFL == "Y") %>%
  group_by(TRTA, AEBODSYS, AEDECOD) %>%
  summarise(n_pts = n_distinct(USUBJID), .groups = "drop") %>%
  left_join(n_patients, by = "TRTA") %>%
  mutate(pct = round(n_pts / N * 100, 1), fmt = paste0(n_pts, " (", pct, "%)")) %>%
  select(TRTA, AEBODSYS, AEDECOD, fmt) %>%
  pivot_wider(names_from = TRTA, values_from = fmt, values_fill = "0 (0.0%)")

# calculate max percentage across arms for threshold filtering
ae_max_pct <- adae_te %>%
  filter(AOCCPFL == "Y") %>%
  group_by(TRTA, AEBODSYS, AEDECOD) %>%
  summarise(n_pts = n_distinct(USUBJID), .groups = "drop") %>%
  left_join(n_patients, by = "TRTA") %>%
  mutate(pct = round(n_pts / N * 100, 1)) %>%
  group_by(AEBODSYS, AEDECOD) %>%
  summarise(max_pct = max(pct), .groups = "drop")

# filter to >= 5% threshold and save
ae_by_system %>%
  left_join(ae_max_pct, by = c("AEBODSYS", "AEDECOD")) %>%
  filter(max_pct >= 5) %>%
  select(-max_pct) %>%
  arrange(AEBODSYS, AEDECOD) %>%
  rename(`Body System` = AEBODSYS, `Preferred Term` = AEDECOD) %>%
  flextable() %>%
  set_caption("Table 4. Treatment-Emergent Adverse Events >= 5% in Any Arm by Body System and Preferred Term (Safety Population)") %>%
  bold(part = "header") %>%
  bg(bg = "#D6E4F0", part = "header") %>%
  autofit() %>%
  save_as_image(path = "output/Table4_AE_ByBodySystem.png", zoom = 2)

cat("Table 4 saved\n")


#verify output
cat("\n=== Deliverables ===\n")
list.files("output/")



