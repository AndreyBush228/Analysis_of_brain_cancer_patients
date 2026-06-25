import matplotlib.pyplot as plt
import seaborn as sns
from research.config import figures_dir, tables_dir
from research.data_units import load_data, get_numeric_cols

df = load_data()

plt.figure(figsize=(10, 6))
df['Age'].hist(bins=25)
plt.savefig(figures_dir / "age_distribution.png")

numeric_cols = get_numeric_cols(df)
corr_matrix = df[numeric_cols].corr()
corr_matrix.to_csv(tables_dir / "correlation_matrix.csv")