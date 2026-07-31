"""
Retrain ALL SIX models from the CSVs in data/ and rebuild models/run_state.pkl.

This is the "retrain from scratch" entry point. It is OPTIONAL -- the shipped
run_state.pkl already contains trained models, and run_docking_selective.py /
the notebook reuse them without retraining. Run this only when you want to
regenerate the models (e.g. after adding data to the CSVs).

Models trained (all scaffold-grouped 5-fold CV, HistGradientBoosting):
  ic50_bundle        <- data/hdac8_ic50_clean_MERGED.csv   (pIC50)
  dock_bundle        <- data/hdac8_dock_clean_MERGED.csv   (r_i_docking_score)
  hdac1_bundle       <- data/hdac1_ic50_clean.csv          (pIC50)
  hdac6_bundle       <- data/hdac6_ic50_clean.csv          (pIC50)
  hdac1_dock_bundle  <- data/hdac1_dock_clean.csv          (r_i_docking_score)
  hdac6_dock_bundle  <- data/hdac6_dock_clean.csv          (r_i_docking_score)

The rebuilt run_state.pkl also carries ic50_df (the HDAC8 IC50 training frame with
an is_hydroxamate column) because the GA seed selection reads it. Existing derived
frames (hits_df etc.) are NOT recomputed here -- run_docking_selective.py / the
notebook regenerate those from the fresh models.

Usage (from scripts/, with ../data present):
    python3 train_all_models.py
Or from anywhere:
    python3 train_all_models.py --data-dir /path/to/data --out /path/to/run_state.pkl
"""
import os
import sys
import pickle
import argparse
import warnings
warnings.filterwarnings('ignore')
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline import train_model_from_df, is_hydroxamate
from rdkit import Chem


TARGETS = [
    # bundle_key,        filename,                      value_col
    ('ic50_bundle',      'hdac8_ic50_clean_MERGED.csv', 'pIC50'),
    ('dock_bundle',      'hdac8_dock_clean_MERGED.csv', 'r_i_docking_score'),
    ('hdac1_bundle',     'hdac1_ic50_clean.csv',        'pIC50'),
    ('hdac6_bundle',     'hdac6_ic50_clean.csv',        'pIC50'),
    ('hdac1_dock_bundle', 'hdac1_dock_clean.csv',       'r_i_docking_score'),
    ('hdac6_dock_bundle', 'hdac6_dock_clean.csv',       'r_i_docking_score'),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', default='../data', help='directory holding the training CSVs')
    ap.add_argument('--out', default='../models/run_state.pkl', help='output pickle path')
    args = ap.parse_args()

    state = {}
    print("Retraining all six models (scaffold-grouped 5-fold CV):")
    for key, fname, value_col in TARGETS:
        path = os.path.join(args.data_dir, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing training file: {path}")
        df = pd.read_csv(path)
        state[key] = train_model_from_df(df, value_col=value_col, label=key)

    # Build ic50_df (HDAC8 IC50 frame + is_hydroxamate) for GA seed selection.
    ic50_df = pd.read_csv(os.path.join(args.data_dir, 'hdac8_ic50_clean_MERGED.csv'))
    ic50_df['is_hydroxamate'] = ic50_df['canonical_smiles'].apply(
        lambda s: (lambda m: is_hydroxamate(m) if m is not None else False)(Chem.MolFromSmiles(s)))
    state['ic50_df'] = ic50_df

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'wb') as f:
        pickle.dump(state, f)

    print("\nSummary (scaffold-CV R2):")
    for key, _, _ in TARGETS:
        r2 = state[key]['cv_r2']
        n = len(state[key]['train_fps'])
        print(f"  {key:18s} R2={r2:.3f}  n={n:,}")
    print(f"\nSaved {args.out}  (ic50_df: {len(ic50_df):,} HDAC8 IC50 rows)")


if __name__ == '__main__':
    main()
