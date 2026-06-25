from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

project_dir = Path(__file__).resolve().parents[1]

data_dir = project_dir / "data"
reports_dir = project_dir / "results"
figures_dir = reports_dir / "figures"
tables_dir = reports_dir / "tables"

figures_dir.mkdir(parents=True, exist_ok=True)
tables_dir.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("DATA LOADING")
print("=" * 80)

df = pd.read_csv(data_dir / "CCGA_clinical_mRNAseq325.csv")

print(f"\nData loaded successfully")
print(f"Shape: {df.shape[0]} rows, {df.shape[1]} columns")
print()

if 'date' in df.columns:
    df['date'] = pd.to_datetime(df['date'], errors='coerce')

print("\n" + "=" * 80)
print("GENERAL STATISTICS")
print("=" * 80)

numeric_cols = df.select_dtypes(include=[np.number]).columns
df[numeric_cols].describe(include="all").to_csv(
    tables_dir / "numeric_summary.csv"
)
print(f"Numeric statistics saved to {tables_dir / 'numeric_summary.csv'}")

categorical_cols = df.select_dtypes(include=['object', 'string']).columns
if len(categorical_cols) > 0:
    df[categorical_cols].describe(include="all").to_csv(
        tables_dir / "categorical_summary.csv"
    )
    print(f"Categorical statistics saved to {tables_dir / 'categorical_summary.csv'}")

print("\n" + "=" * 80)
print("SURVIVAL STATISTICS")
print("=" * 80)

if 'OS' in df.columns:
    os_data = df['OS'].dropna()
    survival_stats = {
        'metric': ['count', 'mean', 'median', 'min', 'max', 'std'],
        'value': [
            len(os_data),
            os_data.mean(),
            os_data.median(),
            os_data.min(),
            os_data.max(),
            os_data.std()
        ]
    }
    survival_df = pd.DataFrame(survival_stats)
    print("\nOverall Survival (OS):")
    print(survival_df.to_string(index=False))
    survival_df.to_csv(tables_dir / "survival_stats.csv", index=False)

if 'Censor (1=dead,0=alive)' in df.columns:
    status = df['Censor (1=dead,0=alive)']
    status_stats = {
        'status': ['alive', 'dead', 'unknown'],
        'count': [
            (status == 0).sum(),
            (status == 1).sum(),
            status.isna().sum()
        ]
    }
    status_df = pd.DataFrame(status_stats)
    print("\nPatient status:")
    print(status_df.to_string(index=False))
    status_df.to_csv(tables_dir / "patient_status.csv", index=False)

print("\n" + "=" * 80)
print("MUTATION STATISTICS")
print("=" * 80)

mutation_cols = ['IDH(DNA&RNA)', 'TP53', 'PTEN', 'EGFR', 'ATRX', 'IDH1-DNA', 'IDH2-DNA']
mutation_stats = []

for col in mutation_cols:
    if col in df.columns:
        if df[col].dtype in ['int64', 'float64']:
            count = df[col].sum()
            if count > 0:
                pct = (count / len(df)) * 100
                mutation_stats.append({
                    'gene': col.replace('(DNA&RNA)', ''),
                    'mutated_count': count,
                    'mutated_percent': round(pct, 2),
                    'wild_type_count': len(df) - count,
                    'wild_type_percent': round(100 - pct, 2)
                })

if mutation_stats:
    mutation_df = pd.DataFrame(mutation_stats)
    print("\nMutation frequency:")
    print(mutation_df.to_string(index=False))
    mutation_df.to_csv(tables_dir / "mutation_stats.csv", index=False)
else:
    print("Mutation columns not found")

print("\n" + "=" * 80)
print("GROUP STATISTICS")
print("=" * 80)

if 'Gender' in df.columns:
    gender_stats = df.groupby('Gender').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nStatistics by gender:")
    print(gender_stats)
    gender_stats.to_csv(tables_dir / "gender_stats.csv")

if 'Grade' in df.columns:
    grade_stats = df.groupby('Grade').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nStatistics by grade:")
    print(grade_stats)
    grade_stats.to_csv(tables_dir / "grade_stats.csv")

if 'IDH(DNA&RNA)' in df.columns:
    idh_stats = df.groupby('IDH(DNA&RNA)').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nStatistics by IDH status:")
    print(idh_stats)
    idh_stats.to_csv(tables_dir / "idh_stats.csv")

if 'TCGA subtype' in df.columns:
    subtype_stats = df.groupby('TCGA subtype').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nStatistics by TCGA subtype:")
    print(subtype_stats)
    subtype_stats.to_csv(tables_dir / "subtype_stats.csv")

print("\n" + "=" * 80)
print("PLOTTING")
print("=" * 80)

plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

if 'Age' in df.columns:
    age_data = df['Age'].dropna()
    plt.figure(figsize=(10, 6))
    plt.hist(age_data, bins=25, edgecolor='black', alpha=0.7)
    plt.axvline(age_data.median(), color='red', linestyle='--',
                label=f'Median: {age_data.median():.0f}')
    plt.axvline(age_data.mean(), color='blue', linestyle='--',
                label=f'Mean: {age_data.mean():.0f}')
    plt.title('Age distribution', fontsize=14)
    plt.xlabel('Age (years)')
    plt.ylabel('Number of patients')
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "age_distribution.png", dpi=300)
    plt.close()
    print("age_distribution.png")

if 'OS' in df.columns:
    os_data = df['OS'].dropna()
    plt.figure(figsize=(10, 6))
    plt.hist(os_data, bins=30, edgecolor='black', alpha=0.7)
    plt.axvline(os_data.median(), color='red', linestyle='--',
                label=f'Median: {os_data.median():.0f}')
    plt.axvline(os_data.mean(), color='blue', linestyle='--',
                label=f'Mean: {os_data.mean():.0f}')
    plt.title('Survival distribution', fontsize=14)
    plt.xlabel('Time (days)')
    plt.ylabel('Number of patients')
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_distribution.png", dpi=300)
    plt.close()
    print("survival_distribution.png")

if mutation_stats:
    mutation_df = pd.DataFrame(mutation_stats)
    plt.figure(figsize=(10, 6))
    bars = plt.bar(mutation_df['gene'], mutation_df['mutated_percent'],
                   color='lightgreen', edgecolor='black')
    plt.title('Mutation frequency', fontsize=14)
    plt.xlabel('Gene')
    plt.ylabel('Mutated patients (%)')
    plt.ylim(0, max(mutation_df['mutated_percent']) * 1.2)
    for bar, pct in zip(bars, mutation_df['mutated_percent']):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                 f'{pct:.1f}%', ha='center', va='bottom')
    plt.tight_layout()
    plt.savefig(figures_dir / "mutation_frequency.png", dpi=300)
    plt.close()
    print("mutation_frequency.png")

if 'IDH(DNA&RNA)' in df.columns and 'OS' in df.columns:
    df_idh = df[['IDH(DNA&RNA)', 'OS']].dropna()
    df_idh['IDH'] = df_idh['IDH(DNA&RNA)'].map({0: 'IDH-wildtype', 1: 'IDH-mutant'})
    plt.figure(figsize=(8, 6))
    df_idh.boxplot(column='OS', by='IDH')
    plt.title('Survival by IDH status', fontsize=14)
    plt.suptitle('')
    plt.xlabel('IDH status')
    plt.ylabel('Survival (days)')
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_by_idh.png", dpi=300)
    plt.close()
    print("survival_by_idh.png")

if 'Grade' in df.columns and 'OS' in df.columns:
    df_grade = df[['Grade', 'OS']].dropna()
    plt.figure(figsize=(8, 6))
    df_grade.boxplot(column='OS', by='Grade')
    plt.title('Survival by grade', fontsize=14)
    plt.suptitle('')
    plt.xlabel('Grade')
    plt.ylabel('Survival (days)')
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_by_grade.png", dpi=300)
    plt.close()
    print("survival_by_grade.png")

numeric_for_corr = df.select_dtypes(include=[np.number]).columns
if len(numeric_for_corr) > 1:
    corr_data = df[numeric_for_corr].copy()
    for col in corr_data.columns:
        median_val = corr_data[col].median()
        corr_data[col] = corr_data[col].fillna(median_val)
    plt.figure(figsize=(12, 10))
    corr_matrix = corr_data.corr()
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, mask=mask, annot=True, cmap='coolwarm',
                center=0, fmt='.2f', square=True)
    plt.title('Correlation matrix', fontsize=14)
    plt.tight_layout()
    plt.savefig(figures_dir / "correlation_matrix.png", dpi=300)
    plt.close()
    print("correlation_matrix.png")
    corr_matrix.to_csv(tables_dir / "correlation_matrix.csv")
    print(f"correlation_matrix.csv saved to {tables_dir}")

if 'Censor (1=dead,0=alive)' in df.columns:
    status = df['Censor (1=dead,0=alive)']
    status_counts = {
        'Alive': (status == 0).sum(),
        'Dead': (status == 1).sum(),
        'Unknown': status.isna().sum()
    }
    status_counts = {k: v for k, v in status_counts.items() if v > 0}
    plt.figure(figsize=(8, 8))
    plt.pie(status_counts.values(), labels=status_counts.keys(),
            autopct='%1.1f%%', startangle=90)
    plt.title('Patient status', fontsize=14)
    plt.tight_layout()
    plt.savefig(figures_dir / "patient_status_pie.png", dpi=300)
    plt.close()
    print("patient_status_pie.png")

if 'Age' in df.columns and 'OS' in df.columns:
    df_age_os = df[['Age', 'OS']].dropna()
    plt.figure(figsize=(10, 6))
    plt.scatter(df_age_os['Age'], df_age_os['OS'], alpha=0.6)
    plt.title('Age vs Survival', fontsize=14)
    plt.xlabel('Age (years)')
    plt.ylabel('Survival (days)')
    z = np.polyfit(df_age_os['Age'], df_age_os['OS'], 1)
    p = np.poly1d(z)
    plt.plot(df_age_os['Age'].sort_values(),
             p(df_age_os['Age'].sort_values()),
             color='red', linestyle='--', label='Trend')
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "age_vs_survival.png", dpi=300)
    plt.close()
    print("age_vs_survival.png")

if 'TCGA subtype' in df.columns and 'OS' in df.columns:
    df_subtype = df[['TCGA subtype', 'OS']].dropna()
    subtype_order = df_subtype.groupby('TCGA subtype')['OS'].median().sort_values(ascending=False).index
    plt.figure(figsize=(10, 6))
    df_subtype.boxplot(column='OS', by='TCGA subtype')
    plt.title('Survival by TCGA subtype', fontsize=14)
    plt.suptitle('')
    plt.xlabel('TCGA subtype')
    plt.ylabel('Survival (days)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_by_subtype.png", dpi=300)
    plt.close()
    print("survival_by_subtype.png")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

if 'Age' in df.columns:
    ax = axes[0, 0]
    df['Age'].dropna().hist(ax=ax, bins=20, edgecolor='black', alpha=0.7)
    ax.set_title('Age distribution')
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Count')
    ax.axvline(df['Age'].median(), color='red', linestyle='--',
               label=f'Median: {df["Age"].median():.0f}')
    ax.legend()

if 'OS' in df.columns:
    ax = axes[0, 1]
    df['OS'].dropna().hist(ax=ax, bins=30, edgecolor='black', alpha=0.7)
    ax.set_title('Survival distribution')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Count')
    ax.axvline(df['OS'].median(), color='red', linestyle='--',
               label=f'Median: {df["OS"].median():.0f}')
    ax.legend()

if mutation_stats:
    ax = axes[1, 0]
    bars = ax.bar(mutation_df['gene'], mutation_df['mutated_percent'],
                  color='lightgreen', edgecolor='black')
    ax.set_title('Mutation frequency')
    ax.set_ylabel('Mutated (%)')
    for bar, pct in zip(bars, mutation_df['mutated_percent']):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f'{pct:.1f}%', ha='center', va='bottom')

if 'Gender' in df.columns:
    ax = axes[1, 1]
    gender_counts = df['Gender'].value_counts()
    ax.pie(gender_counts.values, labels=[f'Male ({gender_counts.get("M", 0)})',
                                         f'Female ({gender_counts.get("F", 0)})'],
           autopct='%1.1f%%', startangle=90)
    ax.set_title('Gender distribution')

plt.suptitle('Dashboard Overview', fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(figures_dir / "dashboard_overview.png", dpi=300, bbox_inches='tight')
plt.close()
print("dashboard_overview.png")
