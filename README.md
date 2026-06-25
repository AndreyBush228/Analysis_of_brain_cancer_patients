# Brain Cancer Clinical Analytics

**Анализ клинических данных пациентов с глиомами**

An end-to-end data analytics and machine learning project focused on glioma and glioblastoma patient data from the Chinese Glioma Genome Atlas (CGGA) database.

The project combines Python-based data processing, exploratory data analysis, survival analysis, Cox proportional hazards regression, and interactive Yandex DataLens dashboards to analyze survival outcomes, genetic mutations, and prognostic factors in brain cancer patients.

---

## Project Overview

This project analyzes clinical and genomic data from brain cancer patients from three perspectives:

- **Clinical Data Analysis** — patient demographics, tumor histology, WHO grade, and treatment history (radio- and chemotherapy).

- **Survival Analysis** — overall survival (OS) distribution, Kaplan-Meier survival curves, and group comparisons by IDH status, tumor grade, and molecular subtypes.

- **Prognostic Factor Analysis** — Cox proportional hazards regression to identify key risk factors affecting patient survival, including IDH mutation status, tumor grade, age, gender, and genetic markers (TP53, PTEN, EGFR, ATRX, EZH2).

- **Genetic Mutation Analysis** — mutation frequency analysis for key glioma-associated genes (IDH1, IDH2, TP53, PTEN, EGFR, ATRX, EZH2) and their association with patient outcomes.

- **Yandex DataLens Dashboard** — interactive business intelligence dashboard for exploring key clinical indicators, survival metrics, and mutation patterns.

The goal of this project is to demonstrate a complete analytics workflow: from raw data to clean datasets, statistical models, and a final dashboard suitable for portfolio presentation.

---

## Key Findings

| Factor | Hazard Ratio | 95% CI | p-value | Effect |
|--------|--------------|--------|---------|--------|
| **WHO Grade** | 2.59 | 2.00 - 3.34 | <0.001 | ⬆️ Increases risk |
| **IDH Mutation** | 0.45 | 0.29 - 0.71 | <0.001 | ⬇️ Decreases risk |
| **Radiotherapy** | 0.47 | 0.32 - 0.69 | <0.001 | ⬇️ Decreases risk |
| **Chemotherapy** | 0.63 | 0.43 - 0.92 | 0.017 | ⬇️ Decreases risk |

- **Concordance Index: 0.80** — excellent model performance
- IDH-mutant patients have **~55% lower risk** of death
- Grade 4 tumors increase risk **~2.6x** compared to lower grades

---

## Tools & Technologies

| Category | Tools |
|----------|-------|
| **Data Processing** | Python, Pandas, NumPy |
| **Statistical Analysis** | SciPy, Lifelines (Cox Regression) |
| **Visualization** | Matplotlib, Seaborn |
| **Machine Learning** | Scikit-learn |
| **Dashboard** | Yandex DataLens |
| **Version Control** | Git, GitHub |

---

## Project Structure
