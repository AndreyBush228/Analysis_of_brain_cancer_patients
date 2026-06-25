from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from sklearn.preprocessing import StandardScaler
import warnings

warnings.filterwarnings('ignore')

project_dir = Path(__file__).resolve().parents[1]

processed_data_dir = project_dir / "data" / "processed"
figures_dir = project_dir / "results" / "figures"
tables_dir = project_dir / "results" / "tables"

figures_dir.mkdir(parents=True, exist_ok=True)
tables_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(processed_data_dir / "brain_cancer_clean.csv")

print("=" * 80)
print("COX PROPORTIONAL HAZARDS REGRESSION ANALYSIS")
print("=" * 80)
print(f"Data loaded: {df.shape[0]} patients, {df.shape[1]} features")
print()

cox_df = df.copy()

cox_df = cox_df.dropna(subset=['OS', 'Censor (1=dead,0=alive)'])

old_names = {
    'IDH(DNA&RNA)': 'IDH_mutation',
    'Censor (1=dead,0=alive)': 'status'
}

for old, new in old_names.items():
    if old in cox_df.columns:
        cox_df[new] = cox_df[old]
        cox_df = cox_df.drop(columns=[old])

duration_col = 'OS'
event_col = 'status'

features = [
    'Age',
    'Gender',
    'Grade',
    'IDH_mutation',
    'TP53',
    'PTEN',
    'EGFR',
    'ATRX',
    'EZH2',
    'Radio',
    'Chemo',
    'mutation_count'
]

existing_features = [col for col in features if col in cox_df.columns]
print(f"Features used in Cox model: {len(existing_features)}")
print(existing_features)
print()

numeric_cols = ['Age', 'Grade', 'mutation_count']
for col in numeric_cols:
    if col in cox_df.columns:
        cox_df[col] = pd.to_numeric(cox_df[col], errors='coerce')

for col in ['Gender', 'IDH_mutation', 'TP53', 'PTEN', 'EGFR', 'ATRX', 'EZH2', 'Radio', 'Chemo']:
    if col in cox_df.columns:
        cox_df[col] = pd.to_numeric(cox_df[col], errors='coerce').fillna(0)

scaler = StandardScaler()
cox_df_scaled = cox_df.copy()
for col in ['Age', 'Grade', 'mutation_count']:
    if col in cox_df_scaled.columns:
        cox_df_scaled[col] = scaler.fit_transform(cox_df_scaled[[col]])

cph = CoxPHFitter()
formula = ' + '.join(existing_features)

try:
    cph.fit(cox_df_scaled, duration_col=duration_col, event_col=event_col, formula=formula)
    print("Cox regression completed successfully")
except Exception as e:
    print(f"Error with full model: {e}")
    print("Trying without mutation_count...")
    features_no_mutation = [f for f in existing_features if f != 'mutation_count']
    formula_no_mutation = ' + '.join(features_no_mutation)
    cph.fit(cox_df_scaled, duration_col=duration_col, event_col=event_col, formula=formula_no_mutation)

print()
print("=" * 80)
print("COX REGRESSION RESULTS")
print("=" * 80)

cph.print_summary()

summary_df = cph.summary
summary_df.to_csv(tables_dir / "cox_regression_results.csv")
print(f"\nResults saved to: {tables_dir / 'cox_regression_results.csv'}")

significant = summary_df[summary_df['p'] < 0.05]
print("\n" + "=" * 80)
print(f"SIGNIFICANT FEATURES (p < 0.05): {len(significant)}")
print("=" * 80)

if len(significant) > 0:
    for feature in significant.index:
        hr = significant.loc[feature, 'exp(coef)']
        p = significant.loc[feature, 'p']
        ci_lower = significant.loc[feature, 'exp(coef) lower 95%']
        ci_upper = significant.loc[feature, 'exp(coef) upper 95%']
        direction = "INCREASES risk" if hr > 1 else "DECREASES risk"
        print(f"\n{feature}:")
        print(f"  HR = {hr:.3f} (95% CI: {ci_lower:.3f} - {ci_upper:.3f})")
        print(f"  p = {p:.4f}")
        print(f"  {direction} by {abs(hr - 1) * 100:.1f}%")
else:
    print("No significant features found (p < 0.05)")

print("\n" + "=" * 80)
print("VISUALIZATION")
print("=" * 80)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

cph.plot(ax=axes[0])
axes[0].set_title('Coefficients in Cox Regression', fontsize=14)
axes[0].axvline(x=0, color='red', linestyle='--', alpha=0.5)
axes[0].set_xlabel('log(HR)')
axes[0].tick_params(axis='x', rotation=45)

cph.plot(hazard_ratios=True, ax=axes[1])
axes[1].set_title('Hazard Ratios (HR) with 95% CI', fontsize=14)
axes[1].axvline(x=1, color='red', linestyle='--', alpha=0.5)
axes[1].set_xlabel('Hazard Ratio')
axes[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig(figures_dir / "cox_regression_results.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"Coefficients plot saved to: {figures_dir / 'cox_regression_results.png'}")

print("\n" + "=" * 80)
print("SURVIVAL CURVES FOR KEY FACTORS")
print("=" * 80)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

factors = [
    ('IDH_mutation', [0, 1], 'IDH status', ['Wildtype', 'Mutant']),
    ('Grade', [2, 3, 4], 'Grade', ['Grade 2', 'Grade 3', 'Grade 4']),
    ('Gender', [0, 1], 'Gender', ['Male', 'Female']),
    ('TP53', [0, 1], 'TP53 mutation', ['Wildtype', 'Mutant'])
]

for i, (factor, values, title, labels) in enumerate(factors):
    if factor in cox_df_scaled.columns:
        ax = axes[i // 2, i % 2]
        try:
            cph.plot_partial_effects_on_outcome(factor, values=values, cmap='coolwarm', ax=ax)
            ax.set_title(f'Survival by {title}', fontsize=12)
            ax.set_xlabel('Time (days)')
            ax.set_ylabel('Survival probability')
            ax.legend(labels)
        except Exception as e:
            ax.text(0.5, 0.5, f'Cannot plot {factor}', ha='center', va='center')
            print(f"Warning: Could not plot {factor}: {e}")

plt.tight_layout()
plt.savefig(figures_dir / "survival_curves_cox.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"Survival curves saved to: {figures_dir / 'survival_curves_cox.png'}")

print("\n" + "=" * 80)
print("PROPORTIONAL HAZARDS TEST")
print("=" * 80)

try:
    test_results = proportional_hazard_test(cph, cox_df_scaled, time_transform='rank')
    test_results.print_summary()
except Exception as e:
    print(f"Could not perform test: {e}")

print("\n" + "=" * 80)
print("FOREST PLOT FOR KEY FEATURES")
print("=" * 80)

plt.figure(figsize=(10, 8))

selected_features = []
for feature in summary_df.index:
    if summary_df.loc[feature, 'p'] < 0.1:
        selected_features.append(feature)

if len(selected_features) > 0:
    hr = summary_df.loc[selected_features, 'exp(coef)']
    ci_lower = summary_df.loc[selected_features, 'exp(coef) lower 95%']
    ci_upper = summary_df.loc[selected_features, 'exp(coef) upper 95%']
    p_values = summary_df.loc[selected_features, 'p']

    y_pos = np.arange(len(selected_features))

    plt.errorbar(hr, y_pos, xerr=[hr - ci_lower, ci_upper - hr],
                 fmt='o', capsize=5, color='black', ecolor='gray')

    plt.axvline(x=1, color='red', linestyle='--', alpha=0.5)

    for i, (feature, h, p) in enumerate(zip(selected_features, hr, p_values)):
        color = 'red' if h > 1 else 'blue'
        plt.scatter(h, i, color=color, s=80)
        label = f"{feature}: HR={h:.2f}, p={p:.3f}"
        if p < 0.05:
            label += " *"
        plt.text(h + 0.05, i, label, va='center', fontsize=9)

    plt.yticks(y_pos, selected_features)
    plt.xlabel('Hazard Ratio (HR)')
    plt.title('Forest Plot of Hazard Ratios', fontsize=14)
    plt.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(figures_dir / "forest_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Forest plot saved to: {figures_dir / 'forest_plot.png'}")

print(f"\nResults saved to:")
print(f"  - {tables_dir / 'cox_regression_results.csv'}")
print(f"\nFigures saved to:")
print(f"  - {figures_dir / 'cox_regression_results.png'}")
print(f"  - {figures_dir / 'survival_curves_cox.png'}")
print(f"  - {figures_dir / 'forest_plot.png'}")
