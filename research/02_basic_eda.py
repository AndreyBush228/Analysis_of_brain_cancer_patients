from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

# Пути
project_dir = Path(__file__).resolve().parents[1]

# Папки для данных и результатов
data_dir = project_dir / "data"
reports_dir = project_dir / "results"
figures_dir = reports_dir / "figures"
tables_dir = reports_dir / "tables"

# Создаем папки
figures_dir.mkdir(parents=True, exist_ok=True)
tables_dir.mkdir(parents=True, exist_ok=True)

# Загружаем данные
print("=" * 80)
print("📊 ЗАГРУЗКА ДАННЫХ")
print("=" * 80)

df = pd.read_csv(data_dir / "CCGA_clinical_mRNAseq325.csv")

print(f"\n✅ Данные загружены")
print(f"Размер: {df.shape[0]} строк, {df.shape[1]} столбцов")
print()

# Преобразуем даты (если есть)
# В вашем файле нет дат, но оставим для совместимости
if 'date' in df.columns:
    df['date'] = pd.to_datetime(df['date'], errors='coerce')

# 1. АНАЛИЗ ПРОПУСКОВ
print("=" * 80)
print("🔍 АНАЛИЗ ПРОПУСКОВ")
print("=" * 80)

missing_rows = []

for column in df.columns:
    missing_count = df[column].isna().sum()
    missing_percent = missing_count / len(df) * 100

    missing_rows.append({
        "column": column,
        "missing_count": missing_count,
        "missing_percent": round(missing_percent, 2),
    })

missing_report = pd.DataFrame(missing_rows)
missing_report = missing_report.sort_values(
    by=["missing_percent", "missing_count"],
    ascending=False,
)

print("\nТоп-10 столбцов с пропусками:")
print(missing_report.head(10).to_string(index=False))

# Сохраняем отчет о пропусках
missing_report.to_csv(tables_dir / "missing_values_report.csv", index=False)
print(f"\n✅ Отчет о пропусках сохранен в {tables_dir / 'missing_values_report.csv'}")

# 2. ОБЩАЯ СТАТИСТИКА
print("\n" + "=" * 80)
print("📊 ОБЩАЯ СТАТИСТИКА")
print("=" * 80)

# Числовые столбцы
numeric_cols = df.select_dtypes(include=[np.number]).columns
df[numeric_cols].describe(include="all").to_csv(
    tables_dir / "numeric_summary.csv"
)
print(f"✅ Статистика числовых данных сохранена в {tables_dir / 'numeric_summary.csv'}")

# Категориальные столбцы
categorical_cols = df.select_dtypes(include=['object']).columns
if len(categorical_cols) > 0:
    df[categorical_cols].describe(include="all").to_csv(
        tables_dir / "categorical_summary.csv"
    )
    print(f"✅ Статистика категориальных данных сохранена в {tables_dir / 'categorical_summary.csv'}")

# 3. СТАТИСТИКА ПО ВЫЖИВАЕМОСТИ
print("\n" + "=" * 80)
print("🎯 СТАТИСТИКА ВЫЖИВАЕМОСТИ")
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
    print("\nОбщая выживаемость (OS):")
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
    print("\nСтатус пациентов:")
    print(status_df.to_string(index=False))
    status_df.to_csv(tables_dir / "patient_status.csv", index=False)

# 4. СТАТИСТИКА ПО МУТАЦИЯМ
print("\n" + "=" * 80)
print("🧬 СТАТИСТИКА МУТАЦИЙ")
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
    print("\nЧастота мутаций:")
    print(mutation_df.to_string(index=False))
    mutation_df.to_csv(tables_dir / "mutation_stats.csv", index=False)
else:
    print("⚠️ Столбцы с мутациями не найдены")

# 5. СТАТИСТИКА ПО ГРУППАМ
print("\n" + "=" * 80)
print("📊 СТАТИСТИКА ПО ГРУППАМ")
print("=" * 80)

# По полу
if 'Gender' in df.columns:
    gender_stats = df.groupby('Gender').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nСтатистика по полу:")
    print(gender_stats)
    gender_stats.to_csv(tables_dir / "gender_stats.csv")

# По Grade
if 'Grade' in df.columns:
    grade_stats = df.groupby('Grade').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nСтатистика по степени злокачественности:")
    print(grade_stats)
    grade_stats.to_csv(tables_dir / "grade_stats.csv")

# По IDH статусу
if 'IDH(DNA&RNA)' in df.columns:
    idh_stats = df.groupby('IDH(DNA&RNA)').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nСтатистика по IDH статусу:")
    print(idh_stats)
    idh_stats.to_csv(tables_dir / "idh_stats.csv")

# По TCGA подтипу
if 'TCGA subtype' in df.columns:
    subtype_stats = df.groupby('TCGA subtype').agg({
        'OS': ['count', 'mean', 'median'],
        'Age': ['mean', 'median'] if 'Age' in df.columns else []
    })
    print("\nСтатистика по TCGA подтипам:")
    print(subtype_stats)
    subtype_stats.to_csv(tables_dir / "subtype_stats.csv")

# 6. ВИЗУАЛИЗАЦИЯ
print("\n" + "=" * 80)
print("📈 ПОСТРОЕНИЕ ГРАФИКОВ")
print("=" * 80)

# Настройка стиля
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# 1. Распределение возраста
if 'Age' in df.columns:
    age_data = df['Age'].dropna()

    plt.figure(figsize=(10, 6))
    plt.hist(age_data, bins=25, edgecolor='black', alpha=0.7)
    plt.axvline(age_data.median(), color='red', linestyle='--',
                label=f'Медиана: {age_data.median():.0f}')
    plt.axvline(age_data.mean(), color='blue', linestyle='--',
                label=f'Среднее: {age_data.mean():.0f}')
    plt.title('Распределение возраста пациентов', fontsize=14)
    plt.xlabel('Возраст (лет)')
    plt.ylabel('Количество пациентов')
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "age_distribution.png", dpi=300)
    plt.close()
    print("✅ age_distribution.png")

# 2. Распределение выживаемости
if 'OS' in df.columns:
    os_data = df['OS'].dropna()

    plt.figure(figsize=(10, 6))
    plt.hist(os_data, bins=30, edgecolor='black', alpha=0.7)
    plt.axvline(os_data.median(), color='red', linestyle='--',
                label=f'Медиана: {os_data.median():.0f}')
    plt.axvline(os_data.mean(), color='blue', linestyle='--',
                label=f'Среднее: {os_data.mean():.0f}')
    plt.title('Распределение выживаемости', fontsize=14)
    plt.xlabel('Время (дни)')
    plt.ylabel('Количество пациентов')
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_distribution.png", dpi=300)
    plt.close()
    print("✅ survival_distribution.png")

# 3. Частота мутаций (bar plot)
if mutation_stats:
    mutation_df = pd.DataFrame(mutation_stats)

    plt.figure(figsize=(10, 6))
    bars = plt.bar(mutation_df['gene'], mutation_df['mutated_percent'],
                   color='lightgreen', edgecolor='black')
    plt.title('Частота мутаций', fontsize=14)
    plt.xlabel('Ген')
    plt.ylabel('Мутировавшие пациенты (%)')
    plt.ylim(0, max(mutation_df['mutated_percent']) * 1.2)

    for bar, pct in zip(bars, mutation_df['mutated_percent']):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                 f'{pct:.1f}%', ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(figures_dir / "mutation_frequency.png", dpi=300)
    plt.close()
    print("✅ mutation_frequency.png")

# 4. Boxplot выживаемости по IDH
if 'IDH(DNA&RNA)' in df.columns and 'OS' in df.columns:
    df_idh = df[['IDH(DNA&RNA)', 'OS']].dropna()
    df_idh['IDH'] = df_idh['IDH(DNA&RNA)'].map({0: 'IDH-дикий', 1: 'IDH-мутантный'})

    plt.figure(figsize=(8, 6))
    df_idh.boxplot(column='OS', by='IDH')
    plt.title('Выживаемость по IDH статусу', fontsize=14)
    plt.suptitle('')  # Убираем автоматический заголовок
    plt.xlabel('IDH статус')
    plt.ylabel('Выживаемость (дни)')
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_by_idh.png", dpi=300)
    plt.close()
    print("✅ survival_by_idh.png")

# 5. Boxplot выживаемости по Grade
if 'Grade' in df.columns and 'OS' in df.columns:
    df_grade = df[['Grade', 'OS']].dropna()

    plt.figure(figsize=(8, 6))
    df_grade.boxplot(column='OS', by='Grade')
    plt.title('Выживаемость по степени злокачественности', fontsize=14)
    plt.suptitle('')
    plt.xlabel('Grade')
    plt.ylabel('Выживаемость (дни)')
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_by_grade.png", dpi=300)
    plt.close()
    print("✅ survival_by_grade.png")

# 6. Корреляционная матрица
numeric_for_corr = df.select_dtypes(include=[np.number]).columns
if len(numeric_for_corr) > 1:
    corr_data = df[numeric_for_corr].copy()
    # Заполняем пропуски для корреляции
    for col in corr_data.columns:
        corr_data[col].fillna(corr_data[col].median(), inplace=True)

    plt.figure(figsize=(12, 10))
    corr_matrix = corr_data.corr()
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, mask=mask, annot=True, cmap='coolwarm',
                center=0, fmt='.2f', square=True)
    plt.title('Корреляционная матрица', fontsize=14)
    plt.tight_layout()
    plt.savefig(figures_dir / "correlation_matrix.png", dpi=300)
    plt.close()
    print("✅ correlation_matrix.png")

    # Сохраняем корреляционную матрицу в CSV
    corr_matrix.to_csv(tables_dir / "correlation_matrix.csv")
    print(f"✅ correlation_matrix.csv сохранен в {tables_dir}")

# 7. Pie chart по статусу пациентов
if 'Censor (1=dead,0=alive)' in df.columns:
    status = df['Censor (1=dead,0=alive)']
    status_counts = {
        'Живы': (status == 0).sum(),
        'Умерли': (status == 1).sum(),
        'Неизвестно': status.isna().sum()
    }
    status_counts = {k: v for k, v in status_counts.items() if v > 0}

    plt.figure(figsize=(8, 8))
    plt.pie(status_counts.values(), labels=status_counts.keys(),
            autopct='%1.1f%%', startangle=90)
    plt.title('Статус пациентов', fontsize=14)
    plt.tight_layout()
    plt.savefig(figures_dir / "patient_status_pie.png", dpi=300)
    plt.close()
    print("✅ patient_status_pie.png")

# 8. Возраст vs Выживаемость (scatter)
if 'Age' in df.columns and 'OS' in df.columns:
    df_age_os = df[['Age', 'OS']].dropna()

    plt.figure(figsize=(10, 6))
    plt.scatter(df_age_os['Age'], df_age_os['OS'], alpha=0.6)
    plt.title('Возраст vs Выживаемость', fontsize=14)
    plt.xlabel('Возраст (лет)')
    plt.ylabel('Выживаемость (дни)')

    # Добавляем линию тренда
    z = np.polyfit(df_age_os['Age'], df_age_os['OS'], 1)
    p = np.poly1d(z)
    plt.plot(df_age_os['Age'].sort_values(),
             p(df_age_os['Age'].sort_values()),
             color='red', linestyle='--', label='Тренд')
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures_dir / "age_vs_survival.png", dpi=300)
    plt.close()
    print("✅ age_vs_survival.png")

# 9. Выживаемость по TCGA подтипам (если есть)
if 'TCGA subtype' in df.columns and 'OS' in df.columns:
    df_subtype = df[['TCGA subtype', 'OS']].dropna()
    subtype_order = df_subtype.groupby('TCGA subtype')['OS'].median().sort_values(ascending=False).index

    plt.figure(figsize=(10, 6))
    df_subtype.boxplot(column='OS', by='TCGA subtype')
    plt.title('Выживаемость по TCGA подтипам', fontsize=14)
    plt.suptitle('')
    plt.xlabel('TCGA подтип')
    plt.ylabel('Выживаемость (дни)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(figures_dir / "survival_by_subtype.png", dpi=300)
    plt.close()
    print("✅ survival_by_subtype.png")

# 10. Обзорный дашборд с несколькими графиками
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Возраст
if 'Age' in df.columns:
    ax = axes[0, 0]
    df['Age'].dropna().hist(ax=ax, bins=20, edgecolor='black', alpha=0.7)
    ax.set_title('Распределение возраста')
    ax.set_xlabel('Возраст (лет)')
    ax.set_ylabel('Количество')
    ax.axvline(df['Age'].median(), color='red', linestyle='--',
               label=f'Медиана: {df["Age"].median():.0f}')
    ax.legend()

# Выживаемость
if 'OS' in df.columns:
    ax = axes[0, 1]
    df['OS'].dropna().hist(ax=ax, bins=30, edgecolor='black', alpha=0.7)
    ax.set_title('Распределение выживаемости')
    ax.set_xlabel('Время (дни)')
    ax.set_ylabel('Количество')
    ax.axvline(df['OS'].median(), color='red', linestyle='--',
               label=f'Медиана: {df["OS"].median():.0f}')
    ax.legend()

# Мутации
if mutation_stats:
    ax = axes[1, 0]
    bars = ax.bar(mutation_df['gene'], mutation_df['mutated_percent'],
                  color='lightgreen', edgecolor='black')
    ax.set_title('Частота мутаций')
    ax.set_ylabel('Мутировавшие (%)')
    for bar, pct in zip(bars, mutation_df['mutated_percent']):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f'{pct:.1f}%', ha='center', va='bottom')

# Пол
if 'Gender' in df.columns:
    ax = axes[1, 1]
    gender_counts = df['Gender'].value_counts()
    ax.pie(gender_counts.values, labels=[f'Мужчины ({gender_counts.get("M", 0)})',
                                         f'Женщины ({gender_counts.get("F", 0)})'],
           autopct='%1.1f%%', startangle=90)
    ax.set_title('Распределение по полу')

plt.suptitle('Обзорный дашборд данных', fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(figures_dir / "dashboard_overview.png", dpi=300, bbox_inches='tight')
plt.close()
print("✅ dashboard_overview.png")

# 7. ИТОГОВЫЙ ОТЧЕТ
print("\n" + "=" * 80)
print("📊 ИТОГОВЫЙ ОТЧЕТ")
print("=" * 80)

print(f"\n✅ Все таблицы сохранены в: {tables_dir}")
print(f"✅ Все графики сохранены в: {figures_dir}")

print("\nСозданные файлы:")
print("-" * 60)

# Таблицы
print("\n📄 Таблицы:")
for file in sorted(tables_dir.glob("*.csv")):
    print(f"  • {file.name}")

# Графики
print("\n📊 Графики:")
for file in sorted(figures_dir.glob("*.png")):
    print(f"  • {file.name}")

print("\n" + "=" * 80)
print("✅ EDA АНАЛИЗ ЗАВЕРШЕН")
print("=" * 80)