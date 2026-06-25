from pathlib import Path

project_dir = Path(__file__).resolve().parents[1]

data_dir = project_dir / "data"
processed_data_dir = project_dir / "data" / "processed"
results_dir = project_dir / "results"
figures_dir = results_dir / "figures"
tables_dir = results_dir / "tables"
reports_dir = results_dir

DATA_FILE = data_dir / "CCGA_clinical_mRNAseq325.csv"

for dir_path in [data_dir, processed_data_dir, figures_dir, tables_dir]:
    dir_path.mkdir(parents=True, exist_ok=True)