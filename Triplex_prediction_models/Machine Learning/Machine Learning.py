# Task: Scalable Machine Learning Baseline Benchmarking on Large-scale eccDNA Datasets.
# Previous Objectives & Upgrades:
# - Baseline hyperparameter optimization using RandomizedSearchCV locked to PR-AUC (AUPRC).
# - Extract 4-mer features with frequency normalization.
# - Track comprehensive metrics including MCC and Specificity.
# - Export prediction probabilities and true labels to .npz for publication-grade ROC/PR plotting.

import os
import time
import json
import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import joblib
from tqdm import tqdm
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score,
                             roc_auc_score, average_precision_score, confusion_matrix, matthews_corrcoef)
import xgboost as xgb
from scipy.stats import uniform, randint


# ==========================================
# 1. Core Utility Functions
# ==========================================
def load_data_from_csv(file_path):
    """
    Load sequence data and labels from a CSV file.
    Filters out sequences containing unknown 'N' bases.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"❌ File not found: {file_path}")
    print(f"Loading data from: {file_path} ...")
    df = pd.read_csv(file_path)

    # Identify sequence and label columns dynamically
    text_col = 'sequence' if 'sequence' in df.columns else 'text' if 'text' in df.columns else df.columns[0]
    label_col = 'label' if 'label' in df.columns else df.columns[1]

    # Filter sequences containing 'N' and convert to uppercase
    df = df[~df[text_col].str.contains('N', case=False, na=False)]
    texts = df[text_col].str.upper().tolist()
    labels = df[label_col].astype(int).tolist()
    print(f"Successfully loaded {len(texts)} valid sequences.")
    return texts, labels


def extract_kmer_features(train_texts, test_texts, k=4):
    """
    Extract k-mer count features and perform frequency normalization.
    Returns 256-dimensional sparse matrices for 4-mer features.
    """
    print(f"Extracting {k}-mer feature matrices ({4 ** k}-dimensional)...")
    vectorizer = CountVectorizer(analyzer='char', ngram_range=(k, k))
    X_train_counts = vectorizer.fit_transform(train_texts).astype(np.float32)
    X_test_counts = vectorizer.transform(test_texts).astype(np.float32)

    # Frequency normalization with zero-division safeguard
    train_sums = np.asarray(X_train_counts.sum(axis=1))
    test_sums = np.asarray(X_test_counts.sum(axis=1))

    train_sums = np.maximum(train_sums, 1e-9)
    test_sums = np.maximum(test_sums, 1e-9)

    X_train = X_train_counts.multiply(1.0 / train_sums)
    X_test = X_test_counts.multiply(1.0 / test_sums)
    return X_train, X_test, vectorizer


# ==========================================
# 2. Main Benchmarking Logic
# ==========================================
def run_baseline_benchmark(X_train, y_train, X_test, y_test, output_dir, vectorizer):
    """
    Execute RandomizedSearchCV hyperparameter optimization across baseline models.
    Evaluates test sets and exports metrics, best models, and plot data (.npz).
    """
    os.makedirs(output_dir, exist_ok=True)
    joblib.dump(vectorizer, os.path.join(output_dir, "kmer_vectorizer.pkl"))

    # Configured hyperparameter search spaces
    models_config = {
        "Logistic Regression": {
            "estimator": LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42),
            "params": {'C': uniform(0.1, 10.0)},
            "n_iter": 5
        },
        "Random Forest": {
            "estimator": RandomForestClassifier(class_weight='balanced', random_state=42, n_jobs=-1),
            "params": {
                'n_estimators': [100, 200],
                'max_depth': [10, 20, 30],
                'min_samples_leaf': [1, 2, 4]
            },
            "n_iter": 5  # Restrict iterations for scale (>600k rows)
        },
        "XGBoost": {
            "estimator": xgb.XGBClassifier(objective='binary:logistic', eval_metric='aucpr',
                                           tree_method='hist', device='cuda', random_state=42),
            "params": {
                'n_estimators': randint(100, 300),
                'max_depth': randint(3, 8),
                'learning_rate': uniform(0.01, 0.1)
            },
            "n_iter": 10  # Hardware accelerated via GPU
        }
    }

    results = []
    plot_data = {'y_test': y_test}
    params_log = {}

    print("\n" + "=" * 70)
    print("🚀 Starting Baseline Benchmark Pipeline")
    print("=" * 70)

    pbar = tqdm(models_config.items(), desc="Overall Progress", unit="model")

    for name, config in pbar:
        pbar.set_description(f"Optimizing: {name}")
        start_t = time.time()

        # Hyperparameter search with PR-AUC target metric
        search = RandomizedSearchCV(
            estimator=config["estimator"],
            param_distributions=config["params"],
            n_iter=config["n_iter"],
            scoring='average_precision',
            cv=3,
            n_jobs=-1 if "XGBoost" not in name else 1,  # Keep single thread for GPU XGBoost
            random_state=42,
            verbose=1
        )

        search.fit(X_train, y_train)
        best_model = search.best_estimator_
        params_log[name] = {k: str(v) for k, v in search.best_params_.items()}

        # Save trained estimator weights
        joblib.dump(best_model, os.path.join(output_dir, f"{name.replace(' ', '_')}_model.pkl"))

        # Predict confidence scores/probabilities
        if hasattr(best_model, "predict_proba"):
            y_prob = best_model.predict_proba(X_test)[:, 1]
        else:
            y_prob = best_model.decision_function(X_test)

        plot_data[f"{name}_prob"] = y_prob
        y_pred = (y_prob >= 0.5).astype(int)

        # Calculate academic classification metrics
        acc = accuracy_score(y_test, y_pred)
        mcc = matthews_corrcoef(y_test, y_pred)
        auroc = roc_auc_score(y_test, y_prob)
        auprc = average_precision_score(y_test, y_prob)
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0

        elapsed = time.time() - start_t
        tqdm.write(f"✅ [{name}] Completed! AUPRC: {auprc:.4f} | Runtime: {elapsed / 60:.2f} min")

        results.append({
            "Model": name, "ACC": acc, "F1": f1_score(y_test, y_pred),
            "MCC": mcc, "Spec": spec, "AUROC": auroc, "AUPRC": auprc
        })

    # Save benchmark metrics and plotting arrays
    df_res = pd.DataFrame(results).set_index("Model")
    df_res.to_csv(os.path.join(output_dir, "metrics.csv"))
    np.savez(os.path.join(output_dir, "plot_data.npz"), **plot_data)
    with open(os.path.join(output_dir, "params.json"), "w") as f:
        json.dump(params_log, f, indent=4)

    print("\n" + "=" * 80)
    print(f"🎉 Benchmark execution completed. Artifacts saved in: {output_dir}")
    print("=" * 80)
    print(df_res.to_markdown(floatfmt=".4f"))


# ==========================================
# 3. CLI Entry Point
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run machine learning baseline benchmark on eccDNA datasets.")
    parser.add_argument("--train_path", type=str, required=True, help="Path to training CSV file.")
    parser.add_argument("--test_path", type=str, required=True, help="Path to testing CSV file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save output models and metrics.")
    parser.add_argument("--kmer", type=int, default=4, help="K-mer length for feature extraction (Default: 4).")

    args = parser.parse_args()

    # Create timestamped subdirectory inside designated output directory
    timestamp = datetime.now().strftime('%m%d_%H%M')
    final_output_dir = os.path.join(args.output_dir, f"baseline_v3_{timestamp}")

    try:
        # 1. Load Data
        train_x_raw, train_y = load_data_from_csv(args.train_path)
        test_x_raw, test_y = load_data_from_csv(args.test_path)

        # 2. Extract Features
        X_tr, X_te, vec = extract_kmer_features(train_x_raw, test_x_raw, k=args.kmer)

        # 3. Run Benchmarking
        run_baseline_benchmark(X_tr, np.array(train_y), X_te, np.array(test_y), final_output_dir, vec)

    except Exception as e:
        print(f"\n❌ Fatal Error during execution: {e}")