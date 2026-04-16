Phase 3 Alzheimer's Trial — Biostatistical Analysis
CDISC ADaM Pilot Dataset | R | Demographics | MMRM | Kaplan-Meier | Safety analysis

Background
Alzheimer's disease is a progressive neurodegenerative disorder affecting millions globally. This project replicates the primary biostatistical analyses from a Phase 3 randomized controlled trial testing Xanomeline, a muscarinic receptor agonist delivered via transdermal skin patch against placebo in patients with mild-to-moderate Alzheimer's disease.
The dataset used is the publicly available CDISC ADaM Pilot Project dataset is a real FDA submission dataset widely used in the pharmaceutical industry for training and tool development. It contains data from 254 subjects across three treatment arms: Placebo, Xanomeline Low Dose (54mg), and Xanomeline High Dose (81mg).

Objectives

Confirm baseline comparability across treatment arms (Table 1)
Analyze time to first dermatologic event by treatment arm (Kaplan-Meier)
Estimate the treatment effect on cognitive decline using ADAS-Cog(11) (MMRM)


Dataset

| Dataset | Description | Rows | Key Variables |
|---------|-------------|------|---------------|
| ADSL | Subject-level demographics | 254 | AGE, SEX, RACE, TRT01P, MMSETOT |
| ADTTE | Time to event | 254 | AVAL, CNSR, PARAM, TRTP |
| ADQSADAS | ADAS-Cog questionnaire scores | ~3,000 | AVAL, CHG, BASE, AVISIT, PARAM |

Methods
Table 1 — Baseline Demographics
Descriptive statistics were computed for the Intent-to-Treat (ITT) population. Continuous variables are presented as mean (SD) and categorical variables as n (%). Between-arm differences were tested using ANOVA for continuous variables and chi-square for categorical variables.
Kaplan-Meier Survival Analysis
Time to first dermatologic event was analyzed in the safety population using the Kaplan-Meier estimator. Between-arm differences were tested using the log-rank test. Confidence intervals are shown as shaded bands around each curve.
MMRM — Primary Efficacy Analysis
Change from baseline in ADAS-Cog(11) total score was analyzed using a Mixed Effects Model for Repeated Measures (MMRM) with the following specification:
CHG ~ TRT01P + AVISIT + TRT01P:AVISIT + BASE + us(AVISIT | USUBJID)
Fixed effects: Treatment arm, visit, treatment-by-visit interaction, baseline ADAS-Cog score
Covariance structure: Unstructured (FDA recommended)
Estimation method: REML
Population: ITT
Primary endpoint: Treatment difference at Week 24

Results
Table 1 — Baseline Demographics
Baseline demographic and clinical characteristics were well balanced across treatment arms, confirming the integrity of randomization. Mean age was approximately 75 years, mean baseline MMSE was approximately 21 points (mild-to-moderate Alzheimer's), and the three arms were comparable on sex, race, BMI, and disease duration.
Kaplan-Meier — Dermatologic Safety Signal
Xanomeline was associated with significantly higher rates of dermatologic events compared to placebo. Placebo: 29/86 (34%), Xanomeline Low Dose: 62/84 (74%), Xanomeline High Dose: 61/84 (73%). 
Log-rank p = 8 x 10^-14.
Xanomeline was associated with significantly higher rates of dermatologic events compared to placebo. Events occurred earlier in the High Dose arm, confirming a dose-dependent skin tolerability signal driven by the transdermal delivery mechanism.

MMRM — Efficacy Results
Low Dose vs Placebo at Week 24: -1.65 ADAS-Cog points (p = 0.075, trending toward significance)
High Dose vs Placebo at Week 24: -1.04 ADAS-Cog points (p = 0.285, not significant)
Negative estimates indicate less cognitive decline (improvement) relative to placebo. ADAS-Cog is scored where higher = worse, so negative = better.
Neither result reached conventional statistical significance (p < 0.05). However, the consistent directional trend for both doses and the clinically meaningful effect size for Low Dose (−1.65 points) suggest a potential efficacy signal that may be detectable in a larger adequately powered trial.

Conclusions
Xanomeline showed a consistent but statistically underpowered efficacy signal in slowing cognitive decline in Alzheimer's disease. The Low Dose arm demonstrated a numerically larger treatment benefit than the High Dose arm — a counterintuitive finding explained by the tolerability data: 73% of High Dose patients experienced dermatologic events, with events occurring faster and earlier, likely affecting compliance and cognitive test performance.

The integrated story:

Randomization was successful (Table 1 — balanced baseline)
The drug caused significantly more skin reactions than placebo, dose-dependently (Kaplan-Meier)
The efficacy signal trends positive but the sample size (n=254) was insufficient to confirm it statistically (MMRM)
The Low Dose offers a better benefit-risk profile than High Dose

These findings are consistent with the historical development of Xanomeline. The molecule was later reformulated as KarXT (Xanomeline + Trospium) combining Xanomeline with a peripheral blocker that eliminates the tolerability burden while preserving central efficacy and received FDA approval in 2024 for schizophrenia under the brand name Cobenfy.

Repository Structure
cdiscpilot-analysis/
├── data/                        ← ADaM datasets (.xpt files, not uploaded to GitHub)
├── programs/
│   └── cdiscpilot_analysis.R    ← complete annotated analysis script
├── output/
│   ├── Table1_Demographics.docx
│   ├── Figure1_KaplanMeier.png
│   └── Table2_MMRM_Results.docx
└── README.md

How to Reproduce
1. Clone this repository
2. Open cdiscpilot_analysis.R in RStudio
3. Run the script from top to bottom
All packages are listed at the top — install once, load each session
Data downloads automatically from the public CDISC GitHub repository

Tools & Packages
- R 4.4+ — analysis environment
- haven 2.5+ — reads SAS .xpt files
- tidyverse 2.0+ — data manipulation and plotting
- gtsummary 2.5+ — publication-quality Table 1
- survival 3.8+ — Kaplan-Meier estimation
- mmrm — MMRM model fitting
- flextable 0.9+ — exports tables to Word

Key Statistical Concepts Demonstrated

CDISC ADaM data standards — working with ADSL, ADTTE, ADQS datasets
ITT vs Safety population — appropriate population selection per endpoint type
Kaplan-Meier estimator — non-parametric survival analysis with censoring
Log-rank test — statistical comparison of survival curves
MMRM — handling repeated measures with unstructured covariance (FDA standard)
Baseline adjustment — covariate adjustment to increase precision and power
Analysis flags — using ANL01FL and DTYPE to select pre-specified records
Clinical interpretation — connecting statistical results to drug development decisions


Author
Abhishek Chirag Shah
Senior Associate — Clinical Data & Analytics, Eli Lilly
MS Health Informatics, University of Michigan, Ann Arbor


Dataset: CDISC ADaM Pilot Project — publicly available at https://github.com/cdisc-org/sdtm-adam-pilot-project
