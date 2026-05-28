import sys
import os
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from m6a_motif_project import config

def find_best_shift(pwm_ref, pwm_query, max_shift=2):
    """Return (shift, corr) where shift maximizes Pearson correlation between
    pwm_ref and pwm_query over the overlapping columns. Robust to register slip."""
    W = pwm_ref.shape[1]
    best_corr = -np.inf
    best_shift = 0

    for shift in range(-max_shift, max_shift + 1):
        if shift >= 0:
            ref_slice = pwm_ref[:, shift:]
            qry_slice = pwm_query[:, :W - shift]
        else:
            ref_slice = pwm_ref[:, :W + shift]
            qry_slice = pwm_query[:, -shift:]

        if ref_slice.shape[1] < 2:
            continue

        corr = np.corrcoef(ref_slice.flatten(), qry_slice.flatten())[0, 1]
        if not np.isnan(corr) and corr > best_corr:
            best_corr = corr
            best_shift = shift

    return best_shift, best_corr


def shift_pwm(pwm, shift, bg_freqs):
    """Apply a positional shift to a PWM, padding with background freqs."""
    W = pwm.shape[1]
    aligned = np.tile(bg_freqs[:, None], (1, W))
    if shift >= 0:
        if W - shift > 0:
            aligned[:, shift:] = pwm[:, :W - shift]
    else:
        if W + shift > 0:
            aligned[:, :W + shift] = pwm[:, -shift:]
    return aligned


def align_all_pwms(pwms, ref_idx, max_shift=2):
    """Align every PWM to the reference via best Pearson shift.
    Returns aligned PWMs and per-chain shift values."""
    bg_freqs = np.array([0.25, 0.25, 0.25, 0.25])
    ref = np.array(pwms[ref_idx])
    aligned, shifts = [], []
    for p in pwms:
        p = np.array(p)
        shift, _ = find_best_shift(ref, p, max_shift=max_shift)
        aligned.append(shift_pwm(p, shift, bg_freqs))
        shifts.append(shift)
    return aligned, shifts


def pwm_to_consensus(pwm, alphabet=('A', 'C', 'G', 'T')):
    return "".join(alphabet[i] for i in np.argmax(pwm, axis=0))


def column_information_content(pwm, eps=1e-12):
    """Per-column information content in bits (4-letter alphabet)."""
    return np.log2(4) + np.sum(pwm * np.log2(pwm + eps), axis=0)

def calculate_gelman_rubin_ll(ll_histories):
    """Standard Gelman-Rubin on log-likelihood traces (post burn-in second half)."""
    chains = np.array(ll_histories)
    half_point = chains.shape[1] // 2
    chains = chains[:, half_point:]
    M, N = chains.shape
    chain_means = np.mean(chains, axis=1)
    global_mean = np.mean(chain_means)
    B = (N / (M - 1)) * np.sum((chain_means - global_mean) ** 2)
    chain_vars = np.var(chains, axis=1, ddof=1)
    W = np.mean(chain_vars)
    var_plus = ((N - 1) / N) * W + (1 / N) * B
    return 1.0 if W == 0 else float(np.sqrt(var_plus / W))


def calculate_strict_stability(consensus_list):
    """Original exact-string stability (kept for comparison/reporting)."""
    n = len(consensus_list)
    total_pairs, matching_pairs = 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            total_pairs += 1
            if consensus_list[i] == consensus_list[j]:
                matching_pairs += 1
    return (matching_pairs / total_pairs) * 100 if total_pairs > 0 else 0.0


def calculate_aligned_stability(pwms, max_shift=2):
    """Alignment-aware stability via PWM Pearson correlation across chain pairs.
    Returns dict with mean/min correlation and % of pairs at thresholds."""
    n = len(pwms)
    correlations = []
    for i in range(n):
        for j in range(i + 1, n):
            _, corr = find_best_shift(np.array(pwms[i]), np.array(pwms[j]), max_shift=max_shift)
            correlations.append(corr)
    correlations = np.array(correlations)
    return {
        "mean_corr": float(np.mean(correlations)),
        "min_corr": float(np.min(correlations)),
        "pct_pairs_corr_ge_0.8": float(np.mean(correlations >= 0.8) * 100),
        "pct_pairs_corr_ge_0.9": float(np.mean(correlations >= 0.9) * 100),
    }


def pick_best_chain_index(ll_histories):
    """Chain index with the highest final log-likelihood."""
    final_lls = [h[-1] for h in ll_histories]
    return int(np.argmax(final_lls))

def score_sequence_with_pwm(seq, log_pwm, log_bg, char_to_idx):
    W = log_pwm.shape[1]
    seq = seq.replace('U', 'T')
    L = len(seq)
    if L < W:
        return -np.inf
    seq_arr = np.array([char_to_idx.get(c, 0) for c in seq], dtype=int)
    cols = np.arange(W)
    max_score = -np.inf
    for i in range(L - W + 1):
        win = seq_arr[i:i + W]
        score = np.sum(log_pwm[win, cols]) - np.sum(log_bg[win])
        if score > max_score:
            max_score = score
    return max_score


def compute_auroc_with_learned_pwm(pwm, train_indices, pos_df, neg_df):
    """Score the held-out (non-training) split using the model's learned PWM."""
    alphabet = ['A', 'C', 'G', 'T']
    char_to_idx = {c: i for i, c in enumerate(alphabet)}
    bg_freqs = np.array([0.25, 0.25, 0.25, 0.25])

    log_pwm = np.log(pwm)
    log_bg = np.log(bg_freqs)

    test_indices = pos_df.index.difference(train_indices)
    test_pos = pos_df.loc[test_indices, 'sequence'].tolist()
    test_neg = neg_df.loc[test_indices, 'sequence'].tolist()

    y_true = [1] * len(test_pos) + [0] * len(test_neg)
    y_scores = [
        score_sequence_with_pwm(s, log_pwm, log_bg, char_to_idx)
        for s in test_pos + test_neg
    ]
    return float(roc_auc_score(y_true, y_scores))

def main():
    print("--- Executing Model Evaluation Suite (corrected) ---")
    log_file = os.path.join("results", "logs", "experiment_logs.json")
    pos_csv = config.INTERIM_POSITIVE_CSV
    neg_csv = os.path.join("data", "interim", "negatives", "negative_sequences.csv")

    if not (os.path.exists(log_file) and os.path.exists(pos_csv) and os.path.exists(neg_csv)):
        print("Error: Required scoring logs or datasets are missing.")
        return

    with open(log_file, "r") as f:
        data = json.load(f)

    if "pwms" not in data.get("baseline", {}):
        print("Error: experiment_logs.json does not contain saved PWMs.")
        print("       Re-run scripts/run_structure_aware.py with the updated version first.")
        return

    pos_df = pd.read_csv(pos_csv)
    neg_df = pd.read_csv(neg_csv)

    # Reproduce exact training split used during sampling
    if "meta" in data and "sample_indices" in data["meta"]:
        train_indices = pd.Index(data["meta"]["sample_indices"])
    else:
        train_indices = pos_df.sample(n=1000, random_state=config.RANDOM_SEED).index

    for model_name in ["baseline", "structure"]:
        print(f"\n================ Metrics For: {model_name.upper()} ================")
        consensus_motifs = data[model_name]["consensus"]
        ll_histories = data[model_name]["ll_history"]
        pwms_raw = [np.array(p) for p in data[model_name]["pwms"]]

        # ---- Convergence ---- #
        r_hat_ll = calculate_gelman_rubin_ll(ll_histories)

        # ---- Stability (strict + alignment-aware) ---- #
        strict_stab = calculate_strict_stability(consensus_motifs)
        aligned_stab = calculate_aligned_stability(pwms_raw, max_shift=2)

        # ---- Reference chain + alignment for downstream consensus ---- #
        ref_idx = pick_best_chain_index(ll_histories)
        aligned_pwms, shifts = align_all_pwms(pwms_raw, ref_idx=ref_idx, max_shift=2)
        mean_aligned_pwm = np.mean(aligned_pwms, axis=0)
        mean_aligned_pwm /= mean_aligned_pwm.sum(axis=0, keepdims=True)
        aligned_consensus = pwm_to_consensus(mean_aligned_pwm)

        # ---- AUROC: per chain (using each chain's actual PWM) ---- #
        chain_aurocs = [
            compute_auroc_with_learned_pwm(pwm, train_indices, pos_df, neg_df)
            for pwm in pwms_raw
        ]
        mean_auroc = float(np.mean(chain_aurocs))
        std_auroc = float(np.std(chain_aurocs))
        best_chain_auroc = chain_aurocs[ref_idx]

        # ---- AUROC: aligned-mean PWM (single number for headline) ---- #
        consensus_pwm_auroc = compute_auroc_with_learned_pwm(
            mean_aligned_pwm, train_indices, pos_df, neg_df
        )

        # ---- Mode reporting ---- #
        modes = {}
        for m in consensus_motifs:
            modes[m] = modes.get(m, 0) + 1

        print(f"Reference chain (highest final LL):     Chain {ref_idx} (consensus: {consensus_motifs[ref_idx]})")
        print(f"Aligned-mean PWM consensus:             {aligned_consensus}")
        print(f"")
        print(f"Gelman-Rubin (R-hat on log-likelihood): {r_hat_ll:.4f}")
        print(f"Strict stability (exact-string match):  {strict_stab:.1f}%")
        print(f"Aligned PWM mean pairwise correlation:  {aligned_stab['mean_corr']:.3f}")
        print(f"Aligned PWM min  pairwise correlation:  {aligned_stab['min_corr']:.3f}")
        print(f"Aligned stability (corr >= 0.8):        {aligned_stab['pct_pairs_corr_ge_0.8']:.1f}%")
        print(f"Aligned stability (corr >= 0.9):        {aligned_stab['pct_pairs_corr_ge_0.9']:.1f}%")
        print(f"")
        print(f"AUROC per chain (mean ± std):           {mean_auroc:.4f} ± {std_auroc:.4f}")
        print(f"AUROC of best chain (highest LL):       {best_chain_auroc:.4f}")
        print(f"AUROC of aligned-mean PWM:              {consensus_pwm_auroc:.4f}")
        print(f"")
        print(f"Per-chain shifts vs reference:          {shifts}")
        print(f"Captured Motif Distribution Modes:")
        for motif, count in sorted(modes.items(), key=lambda x: -x[1]):
            print(f"  - {motif}: {count} chains out of {len(consensus_motifs)}")


if __name__ == "__main__":
    main()