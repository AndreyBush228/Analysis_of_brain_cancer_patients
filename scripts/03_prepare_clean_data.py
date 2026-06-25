from research.config import processed_data_dir, tables_dir
from research.data_units import load_data, clean_missing_values, encode_categorical

df = load_data()
df = clean_missing_values(df)
df = encode_categorical(df)

df.to_csv(processed_data_dir / "brain_cancer_clean.csv", index=False)