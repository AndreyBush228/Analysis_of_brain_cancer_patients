import pandas as pd
from research.config import tables_dir
from research.data_units import load_data, get_numeric_cols, get_categorical_cols

df = load_data()

print("=" * 80)
print(f"DATASET OVERVIEW")
print("=" * 80)
print(f"Shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")
print(f"Missing values: {df.isna().sum().sum()}")

numeric_cols = get_numeric_cols(df)
categorical_cols = get_categorical_cols(df)

df[numeric_cols].describe().to_csv(tables_dir / "numeric_summary.csv")
df[categorical_cols].describe().to_csv(tables_dir / "categorical_summary.csv")

print(f"Reports saved to {tables_dir}")