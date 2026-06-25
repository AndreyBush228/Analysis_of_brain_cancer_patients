from pathlib import Path
import numpy as np
import pandas as pd

project_dir = Path(__file__).resolve().parents[1]

processed_data_dir = project_dir / "data" / "processed"
datalens_data_dir = project_dir / "data" / "datalens"
reports_tables_dir = project_dir / "results" / "tables"

datalens_data_dir.mkdir(parents=True, exist_ok=True)
reports_tables_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(processed_data_dir / "brain_cancer_clean.csv")

print("=" * 80)
print("Data for Yandex DataLens loaded")
print(f"Shape: {df.shape[0]} rows, {df.shape[1]} columns")
print()

df['SampleID'] = df['SampleID'].astype(str)

categorical_cols = ['TCGA subtype', 'Histology']
for col in categorical_cols:
    if col in df.columns:
        df[col] = df[col].fillna('Unknown')

df['age_group'] = pd.cut(
    df['Age'],
    bins=[0, 30, 40, 50, 60, 100],
    labels=['<30', '30-40', '40-50', '50-60', '60+'],
    include_lowest=True
)

df['os_years'] = (df['OS'] / 365.25).round(2)

df['log_os'] = np.log1p(df['OS'])

df['survival_status'] = df['Censor (1=dead,0=alive)'].map({0: 'Alive', 1: 'Dead'}).fillna('Unknown')

df['idh_status'] = df['IDH(DNA&RNA)'].map({0: 'Wildtype', 1: 'Mutant'}).fillna('Unknown')

df['grade_label'] = df['Grade'].map({2: 'Grade 2', 3: 'Grade 3', 4: 'Grade 4'}).fillna('Unknown')

mutation_cols = ['IDH(DNA&RNA)', 'TP53', 'PTEN', 'EGFR', 'ATRX', 'EZH2']
mutation_cols_existing = []
for col in mutation_cols:
    if col in df.columns:
        if df[col].dtype in ['int64', 'float64']:
            mutation_cols_existing.append(col)
        else:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

if mutation_cols_existing:
    mutation_count = df[mutation_cols_existing].sum(axis=1)
    df['mutation_count'] = mutation_count
else:
    df['mutation_count'] = 0

mutation_frequency = {}
for col in mutation_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        count = df[col].sum()
        pct = (count / len(df)) * 100 if len(df) > 0 else 0
        mutation_frequency[col.replace('(DNA&RNA)', '')] = round(pct, 2)

df['has_any_mutation'] = (df['mutation_count'] > 0).astype(int)
df['has_any_mutation_label'] = df['has_any_mutation'].map({0: 'No mutations', 1: 'Has mutation'})

df['is_high_grade'] = (df['Grade'] >= 4).astype(int)
df['is_high_grade_label'] = df['is_high_grade'].map({0: 'Low grade (2-3)', 1: 'High grade (4)'})

df['gender_label'] = df['Gender'].map({0: 'Male', 1: 'Female'})

df['treatment_radio'] = df['Radio'].map({0: 'No', 1: 'Yes'}).fillna('Unknown')
df['treatment_chemo'] = df['Chemo'].map({0: 'No', 1: 'Yes'}).fillna('Unknown')

if 'IDH1-DNA' in df.columns and 'IDH2-DNA' in df.columns:
    df['idh1_mutation'] = pd.to_numeric(df['IDH1-DNA'], errors='coerce').fillna(0).astype(int)
    df['idh2_mutation'] = pd.to_numeric(df['IDH2-DNA'], errors='coerce').fillna(0).astype(int)
    df['idh1_or_idh2_mutation'] = ((df['idh1_mutation'] == 1) | (df['idh2_mutation'] == 1)).astype(int)

patient_listings = df.copy()
patient_columns = [
    'SampleID',
    'Gender',
    'gender_label',
    'Age',
    'age_group',
    'Histology',
    'Grade',
    'grade_label',
    'TCGA subtype',
    'OS',
    'os_years',
    'log_os',
    'Censor (1=dead,0=alive)',
    'survival_status',
    'Radio',
    'treatment_radio',
    'Chemo',
    'treatment_chemo',
    'IDH(DNA&RNA)',
    'idh_status',
    'IDH1-DNA',
    'idh1_mutation',
    'IDH2-DNA',
    'idh2_mutation',
    'idh1_or_idh2_mutation',
    'TP53',
    'PTEN',
    'EGFR',
    'ATRX',
    'EZH2',
    'mutation_count',
    'has_any_mutation',
    'has_any_mutation_label',
    'is_high_grade',
    'is_high_grade_label'
]
existing_columns = [col for col in patient_columns if col in patient_listings.columns]
patient_listings = patient_listings[existing_columns]


def create_survival_table(df, group_col, group_name_map=None):
    result = (
        df
        .groupby(group_col, as_index=False)
        .agg(
            patient_count=('SampleID', 'count'),
            median_os=('OS', 'median'),
            mean_os=('OS', 'mean'),
            min_os=('OS', 'min'),
            max_os=('OS', 'max'),
            deaths=('Censor (1=dead,0=alive)', 'sum'),
            alive=('Censor (1=dead,0=alive)', lambda x: (x == 0).sum())
        )
    )
    result['death_rate_pct'] = (result['deaths'] / result['patient_count'] * 100).round(2)
    result['alive_rate_pct'] = (result['alive'] / result['patient_count'] * 100).round(2)
    result['median_os_years'] = (result['median_os'] / 365.25).round(2)

    if group_name_map:
        result[group_col] = result[group_col].map(group_name_map).fillna(result[group_col])

    return result


survival_by_grade = create_survival_table(df, 'Grade')
survival_by_grade['grade_label'] = survival_by_grade['Grade'].map({2: 'Grade 2', 3: 'Grade 3', 4: 'Grade 4'})

survival_by_idh = create_survival_table(df, 'IDH(DNA&RNA)', {0: 'Wildtype', 1: 'Mutant'})
survival_by_idh.rename(columns={'IDH(DNA&RNA)': 'idh_status'}, inplace=True)

survival_by_gender = create_survival_table(df, 'Gender', {0: 'Male', 1: 'Female'})
survival_by_gender.rename(columns={'Gender': 'gender_label'}, inplace=True)

survival_by_subtype = create_survival_table(df, 'TCGA subtype')

survival_by_mutation = create_survival_table(df, 'has_any_mutation', {0: 'No mutations', 1: 'Has mutation'})
survival_by_mutation.rename(columns={'has_any_mutation': 'mutation_status'}, inplace=True)

survival_by_age_group = create_survival_table(df, 'age_group')

mutation_frequency_df = pd.DataFrame([
    {'gene': gene, 'mutation_frequency_pct': pct}
    for gene, pct in mutation_frequency.items()
])
if len(mutation_frequency_df) > 0:
    mutation_frequency_df = mutation_frequency_df.sort_values('mutation_frequency_pct', ascending=False)

clinical_kpi = pd.DataFrame([
    {'metric': 'total_patients', 'value': len(df)},
    {'metric': 'deaths', 'value': df['Censor (1=dead,0=alive)'].sum()},
    {'metric': 'alive', 'value': (df['Censor (1=dead,0=alive)'] == 0).sum()},
    {'metric': 'death_rate_pct', 'value': (df['Censor (1=dead,0=alive)'].sum() / len(df) * 100).round(2)},
    {'metric': 'median_survival_days', 'value': df['OS'].median()},
    {'metric': 'median_survival_years', 'value': (df['OS'].median() / 365.25).round(2)},
    {'metric': 'mean_survival_days', 'value': df['OS'].mean().round(0)},
    {'metric': 'mean_survival_years', 'value': (df['OS'].mean() / 365.25).round(2)},
    {'metric': 'min_survival_days', 'value': df['OS'].min()},
    {'metric': 'max_survival_days', 'value': df['OS'].max()},
    {'metric': 'male_patients', 'value': (df['Gender'] == 0).sum()},
    {'metric': 'female_patients', 'value': (df['Gender'] == 1).sum()},
    {'metric': 'female_pct', 'value': ((df['Gender'] == 1).sum() / len(df) * 100).round(2)},
    {'metric': 'idh_mutant_patients', 'value': df['IDH(DNA&RNA)'].sum() if 'IDH(DNA&RNA)' in df.columns else 0},
    {'metric': 'idh_mutant_pct',
     'value': (df['IDH(DNA&RNA)'].sum() / len(df) * 100).round(2) if 'IDH(DNA&RNA)' in df.columns else 0},
    {'metric': 'grade_4_patients', 'value': (df['Grade'] == 4).sum()},
    {'metric': 'grade_4_pct', 'value': ((df['Grade'] == 4).sum() / len(df) * 100).round(2)},
    {'metric': 'grade_2_patients', 'value': (df['Grade'] == 2).sum()},
    {'metric': 'grade_2_pct', 'value': ((df['Grade'] == 2).sum() / len(df) * 100).round(2)},
    {'metric': 'grade_3_patients', 'value': (df['Grade'] == 3).sum()},
    {'metric': 'grade_3_pct', 'value': ((df['Grade'] == 3).sum() / len(df) * 100).round(2)},
    {'metric': 'patients_with_any_mutation', 'value': df['has_any_mutation'].sum()},
    {'metric': 'patients_with_any_mutation_pct', 'value': (df['has_any_mutation'].sum() / len(df) * 100).round(2)},
    {'metric': 'avg_mutation_count', 'value': df['mutation_count'].mean().round(2)}
])

tables_for_rounding = [
    patient_listings,
    survival_by_grade,
    survival_by_idh,
    survival_by_gender,
    survival_by_subtype,
    survival_by_mutation,
    survival_by_age_group,
    mutation_frequency_df,
    clinical_kpi,
]

for table in tables_for_rounding:
    if len(table) > 0:
        numeric_columns = table.select_dtypes(include=[np.number]).columns
        if len(numeric_columns) > 0:
            table[numeric_columns] = table[numeric_columns].round(2)

patient_listings.to_csv(datalens_data_dir / "patient_listings.csv", index=False)

survival_by_grade.to_csv(datalens_data_dir / "survival_by_grade.csv", index=False)
survival_by_idh.to_csv(datalens_data_dir / "survival_by_idh.csv", index=False)
survival_by_gender.to_csv(datalens_data_dir / "survival_by_gender.csv", index=False)
survival_by_subtype.to_csv(datalens_data_dir / "survival_by_subtype.csv", index=False)
survival_by_mutation.to_csv(datalens_data_dir / "survival_by_mutation.csv", index=False)
survival_by_age_group.to_csv(datalens_data_dir / "survival_by_age_group.csv", index=False)

if len(mutation_frequency_df) > 0:
    mutation_frequency_df.to_csv(datalens_data_dir / "mutation_frequency.csv", index=False)
clinical_kpi.to_csv(datalens_data_dir / "clinical_kpi.csv", index=False)

datalens_tables_report = pd.DataFrame([
    {"table_name": "patient_listings.csv", "rows": len(patient_listings), "columns": patient_listings.shape[1]},
    {"table_name": "survival_by_grade.csv", "rows": len(survival_by_grade), "columns": survival_by_grade.shape[1]},
    {"table_name": "survival_by_idh.csv", "rows": len(survival_by_idh), "columns": survival_by_idh.shape[1]},
    {"table_name": "survival_by_gender.csv", "rows": len(survival_by_gender), "columns": survival_by_gender.shape[1]},
    {"table_name": "survival_by_subtype.csv", "rows": len(survival_by_subtype),
     "columns": survival_by_subtype.shape[1]},
    {"table_name": "survival_by_mutation.csv", "rows": len(survival_by_mutation),
     "columns": survival_by_mutation.shape[1]},
    {"table_name": "survival_by_age_group.csv", "rows": len(survival_by_age_group),
     "columns": survival_by_age_group.shape[1]},
    {"table_name": "mutation_frequency.csv", "rows": len(mutation_frequency_df),
     "columns": mutation_frequency_df.shape[1]},
    {"table_name": "clinical_kpi.csv", "rows": len(clinical_kpi), "columns": clinical_kpi.shape[1]},
])

datalens_tables_report.to_csv(
    reports_tables_dir / "datalens_tables_report.csv",
    index=False,
)

print("=" * 80)
print("Yandex DataLens tables prepared")
print()
print(datalens_tables_report)
print()
print(f"DataLens data saved to: {datalens_data_dir}")
print(f"Report saved to: {reports_tables_dir / 'datalens_tables_report.csv'}")