import pandas as pd
import numpy as np
from research.config import DATA_FILE


def load_data():
    """Загрузка данных"""
    return pd.read_csv(DATA_FILE)


def clean_missing_values(df):
    """Очистка пропусков"""
    df = df.copy()

    df = df.dropna(subset=['OS', 'Censor (1=dead,0=alive)'])

    categorical_cols = ['Gender', 'Grade', 'TCGA subtype', 'Radio', 'Chemo', 'Histology']
    for col in categorical_cols:
        if col in df.columns:
            mode_value = df[col].mode()[0]
            df[col] = df[col].fillna(mode_value)

    binary_cols = ['IDH(DNA&RNA)', 'TP53', 'PTEN', 'EGFR', 'ATRX', 'IDH1-DNA', 'IDH2-DNA', 'EZH2']
    for col in binary_cols:
        if col in df.columns:
            if df[col].dtype in ['int64', 'float64']:
                df[col] = df[col].fillna(0)
            else:
                df[col] = df[col].fillna('0')

    return df


def encode_categorical(df):
    """Кодирование категориальных признаков"""
    df = df.copy()

    df['Gender'] = df['Gender'].map({'M': 0, 'F': 1})
    df['Radio'] = df['Radio'].map({'Yes': 1, 'No': 0, 1: 1, 0: 0})
    df['Chemo'] = df['Chemo'].map({'Yes': 1, 'No': 0, 1: 1, 0: 0})
    df['Grade'] = df['Grade'].astype(float)

    return df


def get_numeric_cols(df):
    """Получение числовых столбцов"""
    return df.select_dtypes(include=[np.number]).columns.tolist()


def get_categorical_cols(df):
    """Получение категориальных столбцов"""
    return df.select_dtypes(include=['object', 'string']).columns.tolist()