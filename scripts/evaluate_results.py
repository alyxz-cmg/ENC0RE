import sys
import os
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from m6a_motif_project import config
from m6a_motif_project.motif.gibbs_sampler import BaselineGibbsSampler

def calculate_gelman_rubin(histories):
    chains = np.array(histories)
    half_point = chains.shape[1] // 2
    chains = chains[:, half_point:]
    M, N = chains.shape
    chain_means = np.mean(chains, axis=1)
    global_mean = np.mean(chain_means)
    B = (N / (M - 1)) * np.sum((chain_means - global_mean) ** 2)
    chain_vars = np.var(chains, axis=1, ddof=1)
    W = np.mean(chain_vars)
    var_plus = ((N - 1) / N) * W + (1 / N) * B
    return 1.0 if W == 0 else np.sqrt(var_plus / W)

def calculate_pairwise_stability(consensus_list):
    total_pairs, matching_pairs = 0, 0
    n = len(consensus_list)
    for i in range(n):
        for j in range(i + 1, n):
            total_pairs += 1
            if consensus_list[i] == consensus_list[j]:
                matching_pairs += 1
    return (matching_pairs / total_pairs) * 100 if total_pairs > 0 else 0.0

def score_sequence_with_pwm(seq, pwm, bg_freqs, char_to_idx):
    """Scores a sequence string by finding its maximum log-likelihood ratio window."""
    W = pwm.shape[1]
    max_score = -np.inf
    # Transcribe U to T if present
    seq = seq.replace('U', 'T')
    
    for i in range(len(seq) - W + 1):
        window = [char_to_idx.get(c, 0) for c in seq[i:i+W]]
        ll_motif = np.sum(np.log(pwm[window, np.arange(W)]))
        ll_bg = np.sum(np.log(bg_freqs[window]))
        score = ll_motif - ll_bg
        if score > max_score:
            max_score = score
    return max_score

def compute_auroc_evaluation(best_consensus, data_model_type, train_indices, pos_df, neg_df):
    """Evaluates the predictive capability of the model mode on a 20% held-out split."""
    # 1. Identify which rows belong to the 20% held-out test split
    all_indices = pos_df.index
    test_indices = all_indices.difference(train_indices)
    
    test_pos = pos_df.loc[test_indices, 'sequence'].tolist()
    test_neg = neg_df.loc[test_indices, 'sequence'].tolist()
    
    # 2. Reconstruct a proxy PWM from sequences that match our top consensus mode
    alphabet = ['A', 'C', 'G', 'T']
    char_to_idx = {c: i for i, c in enumerate(alphabet)}
    
    # Simple empirical alignment estimation using the target consensus
    matched_windows = []
    for seq in pos_df.loc[train_indices, 'sequence']:
        seq_clean = seq.replace('U', 'T')
        idx = seq_clean.find(best_consensus)
        if idx != -1:
            matched_windows.append([char_to_idx.get(c, 0) for c in seq_clean[idx:idx+5]])
            
    if len(matched_windows) == 0: # Fallback if manual tracking window misaligns
        return 0.5000
        
    counts = np.zeros((4, 5))
    for win in matched_windows:
        for pos, base in enumerate(win):
            counts[base, pos] += 1
            
    pwm = (counts + 1.0) / (counts + 1.0).sum(axis=0, keepdims=True)
    bg_freqs = np.array([0.25, 0.25, 0.25, 0.25]) # Standard background assumption
    
    # 3. Score test cohort
    y_true = [1] * len(test_pos) + [0] * len(test_neg)
    y_scores = [score_sequence_with_pwm(s, pwm, bg_freqs, char_to_idx) for s in test_pos + test_neg]
    
    return roc_auc_score(y_true, y_scores)

def main():
    print("--- Executing Model Evaluation Suite ---")
    log_file = os.path.join("results", "logs", "experiment_logs.json")
    pos_csv = config.INTERIM_POSITIVE_CSV
    neg_csv = os.path.join("data", "interim", "negatives", "negative_sequences.csv")
    
    if not os.path.exists(log_file) or not os.path.exists(pos_csv) or not os.path.exists(neg_csv):
        print("Error: Required scoring logs or datasets are missing.")
        return
        
    with open(log_file, "r") as f:
        data = json.load(f)
        
    pos_df = pd.read_csv(pos_csv)
    neg_df = pd.read_csv(neg_csv)
    
    # Reconstruct the exact training index mapping used in scripts/run_structure_aware.py
    train_indices = pos_df.sample(n=1000, random_state=config.RANDOM_SEED).index
    
    for model_name in ["baseline", "structure"]:
        print(f"\n================ Metrics For: {model_name.upper()} ================")
        consensus_motifs = data[model_name]["consensus"]
        ll_histories = data[model_name]["ll_history"]
        
        r_hat = calculate_gelman_rubin(ll_histories)
        stability_pct = calculate_pairwise_stability(consensus_motifs)
        
        # Find top mode
        modes = {}
        for m in consensus_motifs:
            modes[m] = modes.get(m, 0) + 1
        top_mode = max(modes, key=modes.get)
        
        # Run classification analysis on held out validation subset
        test_auroc = compute_auroc_evaluation(top_mode, model_name, train_indices, pos_df, neg_df)
        
        print(f"Gelman-Rubin Convergence (R-hat): {r_hat:.4f}")
        print(f"Pairwise Consensuses Stability:   {stability_pct:.1f}%")
        print(f"Held-Out 20% Split Test AUROC:    {test_auroc:.4f}")
        print("Captured Motif Distribution Modes:")
        for motif, count in modes.items():
            print(f"  - {motif}: {count} chains out of {len(consensus_motifs)}")

if __name__ == "__main__":
    main()