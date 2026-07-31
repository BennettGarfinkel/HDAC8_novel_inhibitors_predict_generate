"""
Consolidates HDAC_Docking_Inhibition.csv (73,467 rows: HDAC8/HDAC1/HDAC6 IC50 +
docking data from patents, ChEMBL, and internal FragBreed/HTVS runs) with the
existing chembl_hdac8_ic50_raw.csv. Produces clean per-target training sets.
"""
import warnings
warnings.filterwarnings('ignore')
import re
import numpy as np
import pandas as pd
from rdkit import Chem

RAW_PATH = 'HDAC_Docking_Inhibition.csv'

HDAC8_NAMES = ['HDAC8', 'HDAC8 - histone deacetylase 8 (human)', 'Histone deacetylase 8',
               'Histone deacetylase 8 (human)', 'Chain A, Histone Deacetylase 8 (human)']


def strip_stereo_annotation(smi):
    """Some rows carry SMILES extension annotations like ' |w:8.7,2.1|' or ' |r|'
    (enhanced stereochemistry markup) that RDKit's MolFromSmiles chokes on unless
    using the V3000/extended parser. Strip anything from the first ' |' onward."""
    if not isinstance(smi, str):
        return None
    idx = smi.find(' |')
    return smi[:idx] if idx != -1 else smi


def canon(smi):
    s = strip_stereo_annotation(smi)
    if s is None:
        return None
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m is not None else None


def build_ic50_set(df, target_names, label):
    sub = df[(df['target'].isin(target_names)) & (df['activity_type'] == 'IC50') &
             (df['activity_units'] == 'nM')].copy()
    sub['activity_value'] = pd.to_numeric(sub['activity_value'], errors='coerce')
    sub = sub.dropna(subset=['activity_value', 'smiles'])
    sub = sub[sub['activity_value'] > 0]
    sub['canonical_smiles'] = sub['smiles'].apply(canon)
    sub = sub.dropna(subset=['canonical_smiles'])
    sub['pIC50'] = 9 - np.log10(sub['activity_value'])
    agg = sub.groupby('canonical_smiles').agg(
        pIC50=('pIC50', 'mean'), n_measurements=('pIC50', 'count'), pIC50_std=('pIC50', 'std')
    ).reset_index()
    print(f"{label}: {len(sub)} raw rows -> {len(agg)} unique compounds "
          f"(pIC50 {agg['pIC50'].min():.2f}-{agg['pIC50'].max():.2f})")
    return agg


def build_docking_set(df, target_names, label):
    sub = df[(df['target'].isin(target_names)) & (df['activity_type'] == 'docking_score')].copy()
    sub['activity_value'] = pd.to_numeric(sub['activity_value'], errors='coerce')
    sub = sub.dropna(subset=['activity_value', 'smiles'])
    sub['canonical_smiles'] = sub['smiles'].apply(canon)
    sub = sub.dropna(subset=['canonical_smiles'])
    agg = sub.groupby('canonical_smiles').agg(
        r_i_docking_score=('activity_value', 'mean'), n_measurements=('activity_value', 'count')
    ).reset_index()
    print(f"{label}: {len(sub)} raw rows -> {len(agg)} unique compounds "
          f"(score {agg['r_i_docking_score'].min():.2f} to {agg['r_i_docking_score'].max():.2f})")
    return agg


if __name__ == '__main__':
    df = pd.read_csv(RAW_PATH)

    print("=== HDAC8 IC50: merging new data with existing ChEMBL clean set ===")
    hdac8_new = build_ic50_set(df, HDAC8_NAMES, 'HDAC8 (new file)')
    hdac8_old = pd.read_csv('hdac8_ic50_clean_VERIFIED.csv')[['canonical_smiles', 'pIC50', 'n_measurements']]
    combined = pd.concat([hdac8_old, hdac8_new[['canonical_smiles', 'pIC50', 'n_measurements']]])
    # Re-aggregate in case the same compound appears in both sources (weighted mean
    # by original measurement count so a compound with 5 real measurements isn't
    # diluted 1:1 against one with a single measurement).
    combined['weighted_pIC50'] = combined['pIC50'] * combined['n_measurements']
    merged = combined.groupby('canonical_smiles').agg(
        weighted_sum=('weighted_pIC50', 'sum'), n_measurements=('n_measurements', 'sum')
    ).reset_index()
    merged['pIC50'] = merged['weighted_sum'] / merged['n_measurements']
    merged = merged[['canonical_smiles', 'pIC50', 'n_measurements']]
    print(f"Combined HDAC8 IC50 set: {len(hdac8_old)} (old) + {len(hdac8_new)} (new) "
          f"-> {len(merged)} unique compounds after merge")
    merged.to_csv('hdac8_ic50_clean_MERGED.csv', index=False)

    print("\n=== HDAC1 / HDAC6 IC50 (real off-target data, new this session) ===")
    hdac1_ic50 = build_ic50_set(df, ['HDAC1'], 'HDAC1')
    hdac6_ic50 = build_ic50_set(df, ['HDAC6'], 'HDAC6')
    hdac1_ic50.to_csv('hdac1_ic50_clean.csv', index=False)
    hdac6_ic50.to_csv('hdac6_ic50_clean.csv', index=False)

    print("\n=== HDAC8 / HDAC1 / HDAC6 docking scores (kcal/mol) ===")
    hdac8_dock = build_docking_set(df, HDAC8_NAMES + ['HTVS', 'FragBreed'], 'HDAC8 docking')
    hdac1_dock = build_docking_set(df, ['HDAC1'], 'HDAC1 docking')
    hdac6_dock = build_docking_set(df, ['HDAC6'], 'HDAC6 docking')
    hdac8_dock.to_csv('hdac8_dock_clean_MERGED.csv', index=False)
    hdac1_dock.to_csv('hdac1_dock_clean.csv', index=False)
    hdac6_dock.to_csv('hdac6_dock_clean.csv', index=False)

    print("\nAll consolidated datasets saved.")
