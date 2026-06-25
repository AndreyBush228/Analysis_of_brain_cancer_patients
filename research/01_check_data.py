from pathlib import Path
import pandas as pd
import numpy as np

project_dir = Path(__file__).resolve().parents[1]

raw_data_dir = project_dir / "data"
reports_tables_dir = project_dir / "results" / "tables"

reports_tables_dir.mkdir(parents=True, exist_ok=True)

files = {
    "brain_cancer_clinical": "CCGA_clinical_mRNAseq325.csv",
}

overview_rows = []

for dataset_name, file_name in files.items():
    file_path = raw_data_dir / file_name

    if not file_path.exists():
        print(f"File not found: {file_path}")
        print("Check that the CSV file is located in data/")
        continue

    df = pd.read_csv(file_path)

    print("=" * 80)
    print(f"DATASET: {dataset_name}")
    print("=" * 80)
    print(f"File: {file_name}")
    print(f"Shape: {df.shape[0]} rows, {df.shape[1]} columns")
    print()

    rows_count = df.shape[0]
    columns_count = df.shape[1]
    duplicates_count = df.duplicated().sum()
    missing_values_count = df.isna().sum().sum()
    unique_ids = df['SampleID'].nunique() if 'SampleID' in df.columns else 'N/A'

    overview_rows.append(
        {
            "dataset": dataset_name,
            "file_name": file_name,
            "rows": rows_count,
            "columns": columns_count,
            "duplicates": duplicates_count,
            "missing_values_total": missing_values_count,
            "unique_patients": unique_ids,
        }
    )

    print("COLUMNS:")
    print(df.columns.tolist())
    print()

    print("DATA TYPES:")
    print(df.dtypes)
    print()

    print("TOP 15 COLUMNS BY MISSING VALUES:")
    missing = df.isna().sum().sort_values(ascending=False)
    print(missing.head(15))
    print()

    print("STATISTICS FOR NUMERIC COLUMNS:")
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if numeric_cols:
        print(df[numeric_cols].describe())
    else:
        print("No numeric columns")
    print()

    print("TARGET VARIABLES ANALYSIS:")
    if 'OS' in df.columns:
        os_data = df['OS'].dropna()
        print(f"\nOverall Survival (OS):")
        print(f"  Median: {os_data.median():.0f} days")
        print(f"  Mean: {os_data.mean():.0f} days")
        print(f"  Range: {os_data.min()} - {os_data.max()} days")

    if 'Censor (1=dead,0=alive)' in df.columns:
        status = df['Censor (1=dead,0=alive)']
        deaths = status.sum()
        alive = len(status) - deaths - status.isna().sum()
        print(f"\nPatient status:")
        print(f"  Dead: {deaths} ({deaths / len(status) * 100:.1f}%)")
        print(f"  Alive: {alive} ({alive / len(status) * 100:.1f}%)")

    if 'Gender' in df.columns:
        print(f"\nGender distribution:")
        print(df['Gender'].value_counts())

    if 'Grade' in df.columns:
        print(f"\nGrade distribution:")
        print(df['Grade'].value_counts().sort_index())

    print("\nMUTATION FREQUENCY:")
    mutation_cols = ['IDH(DNA&RNA)', 'TP53', 'PTEN', 'EGFR', 'ATRX', 'IDH1-DNA', 'IDH2-DNA']
    for col in mutation_cols:
        if col in df.columns:
            if df[col].dtype in ['int64', 'float64']:
                count = df[col].sum()
            else:
                count = (df[col] == 1).sum() if col in df.columns else 0
            if isinstance(count, (int, float)):
                pct = (count / len(df)) * 100 if len(df) > 0 else 0
                print(f"  {col}: {count} ({pct:.1f}%)")

    print("\nFIRST 5 ROWS:")
    print(df.head())
    print()

    if duplicates_count > 0:
        print(f"Duplicates found: {duplicates_count}")
    else:
        print("No duplicates found")
    print()

    if 'SampleID' in df.columns:
        print(f"Unique SampleID: {df['SampleID'].nunique()} out of {len(df)}")
        if df['SampleID'].nunique() == len(df):
            print("All SampleIDs are unique")
        else:
            print("Duplicate SampleIDs exist")

overview = pd.DataFrame(overview_rows)
output_path = reports_tables_dir / "dataset_overview.csv"
overview.to_csv(output_path, index=False)

print("=" * 80)
print("DATASET OVERVIEW")
print("=" * 80)
print(overview.to_string(index=False))
print()
print(f"Report saved to: {output_path}")

print("\n" + "=" * 80)
print("CORRELATION ANALYSIS")
print("=" * 80)

all_numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

interesting_cols = ['Age', 'OS']
if 'Grade' in all_numeric_cols:
    interesting_cols.append('Grade')

mutation_cols = ['IDH(DNA&RNA)', 'TP53', 'PTEN', 'EGFR', 'ATRX']
for col in mutation_cols:
    if col in all_numeric_cols:
        interesting_cols.append(col)

corr_cols = [col for col in interesting_cols if col in df.columns]

print(f"\nColumns used for correlation: {corr_cols}")

if len(corr_cols) > 1:
    df_corr = df[corr_cols].copy()

    print(f"\nMissing values before filling:")
    print(df_corr.isna().sum())

    for col in df_corr.columns:
        if df_corr[col].isna().sum() > 0:
            mean_val = df_corr[col].mean(skipna=True)
            df_corr[col] = df_corr[col].fillna(mean_val)
            print(f"  {col}: filled {df_corr[col].isna().sum()} missing values with {mean_val:.2f}")

    correlation_matrix = df_corr.corr()

    print("\nCorrelation matrix (key features):")
    print(correlation_matrix.round(3))

    corr_path = reports_tables_dir / "correlation_matrix.csv"
    correlation_matrix.to_csv(corr_path)
    print(f"\nCorrelation matrix saved to: {corr_path}")
else:
    print("Not enough numeric columns for correlation analysis")

print("\n" + "=" * 80)
print("SURVIVAL ANALYSIS BY GROUPS")
print("=" * 80)

if 'IDH(DNA&RNA)' in df.columns and 'OS' in df.columns:
    if df['IDH(DNA&RNA)'].dtype in ['int64', 'float64']:
        print("\nSurvival by IDH status:")
        idh_groups = df.groupby('IDH(DNA&RNA)')['OS'].agg(['median', 'mean', 'count'])
        idh_groups.columns = ['Median (days)', 'Mean (days)', 'Count']
        print(idh_groups)

if 'Grade' in df.columns and 'OS' in df.columns:
    if df['Grade'].dtype in ['int64', 'float64']:
        print("\nSurvival by Grade:")
        grade_groups = df.groupby('Grade')['OS'].agg(['median', 'mean', 'count'])
        grade_groups.columns = ['Median (days)', 'Mean (days)', 'Count']
        print(grade_groups)

if 'Gender' in df.columns and 'OS' in df.columns:
    print("\nSurvival by Gender:")
    gender_groups = df.groupby('Gender')['OS'].agg(['median', 'mean', 'count'])
    gender_groups.columns = ['Median (days)', 'Mean (days)', 'Count']
    print(gender_groups)