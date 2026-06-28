# Brain Cancer Clinical Analytics

**English version** | [**Русская версия**](README_RU.md)

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

## Dashboard Preview

### Survival by Age Group

<img width="1731" height="837" alt="Снимок экрана 2026-06-28 223016" src="https://github.com/user-attachments/assets/31df450f-c08d-4355-be3e-695536eb7a56" />

### Distribution by Molecular Subtypes

<img width="1708" height="445" alt="Снимок экрана 2026-06-28 223556" src="https://github.com/user-attachments/assets/07c434a0-61ba-4c64-b4eb-bb401ea0dc48" />

### Mutation Analysis

<img width="1730" height="600" alt="Снимок экрана 2026-06-28 223757" src="https://github.com/user-attachments/assets/da86a9ed-3162-4fa4-a26a-fbc8cf389bba" />
<img width="1727" height="633" alt="Снимок экрана 2026-06-28 223945" src="https://github.com/user-attachments/assets/c5675064-5111-438d-91d9-c37cfdf3657d" />

### Effect of Treatment

<img width="1732" height="661" alt="Снимок экрана 2026-06-28 224118" src="https://github.com/user-attachments/assets/d0481258-4c7e-487a-9e34-03145d17530d" />


---

## Chart Preview

### Survival Distribution

![Survival Distribution](results/figures/survival_distribution.png)

### Mutation Frequency

![Mutation Frequency](results/figures/mutation_frequency.png)

### Correlation Matrix

![Correlation Matrix](results/figures/correlation_matrix.png)

## Dataset

The project uses clinical and genomic data from the **Chinese Glioma Genome Atlas (CGGA)** database.

**Source:** CGGA mRNAseq 325 dataset  
[https://www.cgga.org.cn/]([https://www.cgga.org.cn/](https://zenodo.org/records/8190378))

Raw data files are not included in this repository because they are stored in `data/`, which is excluded from Git tracking.

### Main source file:

- `CCGA_clinical_mRNAseq325.csv`

### Dataset description:

- **Patients:** 325 glioma patients
- **Features:** 24 clinical and genomic variables
- **Key variables:**
  - Demographics: Age, Gender
  - Clinical: WHO Grade, Histology, TCGA subtype
  - Treatment: Radiotherapy, Chemotherapy
  - Outcomes: Overall Survival (OS), Censoring status
  - Genetics: IDH1/IDH2, TP53, PTEN, EGFR, ATRX, EZH2 mutations

## Project Structure

```
Brain/
├── data/
│   ├── processed/
│   │   └── brain_cancer_clean.csv
│   └── datalens/
│       ├── patient_listings.csv
│       ├── survival_by_grade.csv
│       └── ...
├── results/
│   ├── figures/
│   │   ├── age_distribution.png
│   │   └── ...
│   └── tables/
│       ├── dataset_overview.csv
│       └── ...
├── scripts/
│   ├── config.py
│   ├── main.py
│   ├── 01_check_data.py
│   ├── 02_basic_eda.py
│   └── ...
├── screenshots/
│   └── ...
├── requirements.txt
├── README.md
├── README_RU.md
└── .gitignore
```

---

## Python Pipeline

The project is organized as a sequential pipeline.

### 1. Data Quality Check

```bash
python research/01_check_data.py
```

Checks source datasets, column structure, missing values and basic dataset statistics.

**Outputs:**
- `results/tables/dataset_overview.csv`
- `results/tables/numeric_summary.csv`
- `results/tables/categorical_summary.csv`
- `results/tables/correlation_matrix.csv`

---

### 2. Exploratory Data Analysis (EDA)

```bash
python research/02_basic_eda.py
```

Generates exploratory charts and summary tables for clinical and genomic data.

**Outputs:**
- `results/figures/age_distribution.png`
- `results/figures/survival_distribution.png`
- `results/figures/mutation_frequency.png`
- `results/figures/survival_by_idh.png`
- `results/figures/survival_by_grade.png`
- `results/figures/correlation_matrix.png`
- `results/figures/patient_status_pie.png`
- `results/figures/age_vs_survival.png`
- `results/figures/survival_by_subtype.png`
- `results/figures/dashboard_overview.png`

**Tables:**
- `results/tables/survival_stats.csv`
- `results/tables/patient_status.csv`
- `results/tables/mutation_stats.csv`
- `results/tables/gender_stats.csv`
- `results/tables/grade_stats.csv`
- `results/tables/idh_stats.csv`
- `results/tables/subtype_stats.csv`

---

### 3. Data Cleaning and Feature Engineering

```bash
python research/03_prepare_clean_data.py
```

Creates cleaned datasets and engineered features, including:

- Age groups
- Overall Survival in years
- Survival status labels
- IDH status labels
- Grade labels
- Mutation count
- Any mutation indicator
- High grade indicator
- Treatment indicators (radiotherapy, chemotherapy)
- IDH1/IDH2 mutation indicators

**Outputs:**
- `data/processed/brain_cancer_clean.csv`

---

### 4. Yandex DataLens Table Preparation

```bash
python research/04_prepare_yandex_tables.py
```

Creates analytical CSV tables used in the Yandex DataLens dashboard.

**Outputs:**
- `patient_listings.csv` — detailed patient-level data
- `survival_by_grade.csv` — survival metrics by WHO grade
- `survival_by_idh.csv` — survival metrics by IDH status
- `survival_by_gender.csv` — survival metrics by gender
- `survival_by_subtype.csv` — survival metrics by TCGA subtype
- `survival_by_mutation.csv` — survival metrics by mutation status
- `survival_by_age_group.csv` — survival metrics by age group
- `mutation_frequency.csv` — mutation frequency by gene
- `clinical_kpi.csv` — key clinical indicators

**Output location:** `data/datalens/`

---

### 5. Cox Proportional Hazards Regression

```bash
python research/05_cox_analysis.py
```

Performs Cox regression to identify key prognostic factors affecting patient survival.

**Features tested:**
- Age
- Gender
- Grade
- IDH mutation
- TP53
- PTEN
- EGFR
- ATRX
- EZH2
- Radiotherapy
- Chemotherapy
- Mutation count

**Outputs:**
- `results/tables/cox_regression_results.csv`
- `results/figures/cox_regression_results.png`
- `results/figures/survival_curves_cox.png`
- `results/figures/forest_plot.png`

---

### Cox Regression Results

| Factor | Hazard Ratio | 95% CI | p-value | Effect |
|--------|--------------|--------|---------|--------|
| **WHO Grade** | 2.59 | 2.00 - 3.34 | <0.001 |  Increases risk by 158.5% |
| **IDH Mutation** | 0.45 | 0.29 - 0.71 | <0.001 |  Decreases risk by 54.7% |
| **Radiotherapy** | 0.47 | 0.32 - 0.69 | <0.001 |  Decreases risk by 52.9% |
| **Chemotherapy** | 0.63 | 0.43 - 0.92 | 0.017 |  Decreases risk by 37.5% |

**Concordance Index: 0.80** — excellent model performance

### Clinical Observations

- **IDH mutation is the strongest protective factor**: IDH-mutant patients have ~55% lower risk of death compared to IDH-wildtype patients.
- **WHO Grade is the strongest risk factor**: Grade 4 tumors increase risk ~2.6x compared to lower grades.
- **Standard treatments significantly improve prognosis**: Both radiotherapy and chemotherapy show strong protective effects.
- **Age and gender are not independent prognostic factors** after accounting for other variables.

---

## Yandex DataLens Dashboard

The Yandex DataLens dashboard contains five pages:

### Executive Overview
High-level clinical summary with patient demographics (age distribution, gender), survival overview (median OS, patient status), key clinical KPIs and mutation frequency.

### Survival Analysis
Detailed survival analysis by WHO Grade, IDH status, TCGA subtype, age group, gender and mutation status with median overall survival metrics.

### Genetic Analysis
Genetic landscape overview showing mutation frequency by gene, mutation count distribution, correlation between mutations and survival outcomes by mutation status.

### Patient Demographics
Patient population overview with age distribution, gender ratio, WHO Grade distribution, histological types and TCGA molecular subtypes distribution.

### Treatment Analysis
Treatment impact analysis showing radiotherapy and chemotherapy distribution, survival comparison by treatment status and treatment patterns across different patient groups.

The dashboard file is available here:

```[
https://707.su/bsWq
```

## Key Clinical Questions

This project helps answer questions such as:

- How does IDH mutation status affect patient survival outcomes?
- What is the impact of WHO grade on overall survival?
- Which genetic mutations (TP53, PTEN, EGFR, ATRX, EZH2) are most common in glioma patients?
- How do age and gender influence survival in glioma patients?
- Do radiotherapy and chemotherapy improve survival outcomes?
- Which factors are the strongest independent prognostic factors in glioma?
- How does survival differ across TCGA molecular subtypes?
- What is the relationship between mutation count and prognosis?

---

## Important Limitations

- **Sample size**: 325 patients is relatively small for subgroup analyses.
- **Data completeness**: Some variables (e.g., TP53, PTEN) have significant missing data.
- **Treatment data**: Information on specific treatment protocols is limited.
- **External validation**: Findings should be validated on external cohorts.
- **Patient cohort**: Data comes from a single center in China; results may not generalize to other populations.
- **Proportional hazards**: The Cox model assumes proportional hazards, which may not hold for all covariates.

The analysis should be interpreted as exploratory and hypothesis-generating, not as definitive clinical guidance.

---

## How to Run the Project

### 1. Clone the repository

```bash
git clone https://github.com/AndreyBush228/Analysis_of_brain_cancer_patients.git
cd Brain
```

### 2. Create and activate a virtual environment

```bash
python -m venv .venv
```

**Windows:**
```bash
.venv\Scripts\activate
```

**macOS / Linux:**
```bash
source .venv/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Download the dataset

Download the dataset from CGGA and place the CSV file into:

```
data/
```

Required file:
- `CCGA_clinical_mRNAseq325.csv`

### 5. Run the pipeline

```bash
python research/01_check_data.py
python research/02_basic_eda.py
python research/03_prepare_clean_data.py
python research/04_prepare_yandex_tables.py
python research/05_cox_analysis.py
```

Or use the main entry point:

```bash
python scripts/main.py
```

### 6. Import tables into Yandex DataLens

```
[data/datalens/](https://707.su/bsWq)
```

---

## Project Status

**Completed:**
- Data quality checks
- Exploratory data analysis
- Data cleaning
- Feature engineering
- Cox proportional hazards regression
- Yandex DataLens tables
- Dashboard-ready CSV exports
- Dashboard screenshots

**In Progress:**
- Yandex DataLens dashboard design
- External validation (planned)

---

## Author

Created by Andrey Bushuev (@AndreyBush228) as a portfolio data analytics project.
