import sys
import os
import json
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from m6a_motif_project import config

def calculate_gelman_rubin(histories):
    """Computes the Gelman-Rubin (R-hat) convergence diagnostic on log-likelihood tracks."""
    # Using the second half of the MCMC trajectories to account for burn-in
    chains = np.array(histories)
    half_point = chains.shape[1] // 2
    chains = chains[:, half_point:]
    
    M, N = chains.shape # M chains, N iterations
    
    # 1. Calculate Mean of each chain and overall mean
    chain_means = np.mean(chains, axis=1)
    global_mean = np.mean(chain_means)
    
    # 2. Calculate Between-chain variance (B)
    B = (N / (M - 1)) * np.sum((chain_means - global_mean) ** 2)
    
    # 3. Calculate Within-chain variance (W)
    chain_vars = np.var(chains, axis=1, ddof=1)
    W = np.mean(chain_vars)
    
    # 4. Estimate target variance and calculate R-hat
    var_plus = ((N - 1) / N) * W + (1 / N) * B
    if W == 0:
        return 1.0
    r_hat = np.sqrt(var_plus / W)
    return r_hat

def calculate_pairwise_stability(consensus_list):
    """Measures the percentage of chains that converged to identical consensus text outputs."""
    total_pairs = 0
    matching_pairs = 0
    n = len(consensus_list)
    
    for i in range(n):
        for j in range(i + 1, n):
            total_pairs += 1
            if consensus_list[i] == consensus_list[j]:
                matching_pairs += 1
                
    return (matching_pairs / total_pairs) * 100 if total_pairs > 0 else 0.0

def main():
    print("--- Executing Model Evaluation Suite ---")
    log_file = os.path.join("results", "logs", "experiment_logs.json")
    
    if not os.path.exists(log_file):
        print(f"Error: Log file {log_file} not found. Run experiments first.")
        return
        
    with open(log_file, "r") as f:
        data = json.load(f)
        
    for model_name in ["baseline", "structure"]:
        print(f"\n================ Metrics For: {model_name.upper()} ================")
        consensus_motifs = data[model_name]["consensus"]
        ll_histories = data[model_name]["ll_history"]
        
        # 1. Convergence Diagnostic
        r_hat = calculate_gelman_rubin(ll_histories)
        
        # 2. Pairwise Chain Stability
        stability_pct = calculate_pairwise_stability(consensus_motifs)
        
        # 3. Mode Distribution Tracking
        modes = {}
        for m in consensus_motifs:
            modes[m] = modes.get(m, 0) + 1
            
        print(f"Gelman-Rubin Convergence (R-hat): {r_hat:.4f}")
        print(f"Pairwise Consensuses Stability:   {stability_pct:.1f}%")
        print("Captured Motif Distribution Modes:")
        for motif, count in modes.items():
            print(f"  - {motif}: {count} chains out of {len(consensus_motifs)}")

if __name__ == "__main__":
    main()