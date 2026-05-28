import sys
import os
import pandas as pd
import numpy as np
import json

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from m6a_motif_project.motif.gibbs_sampler import BaselineGibbsSampler, StructureAwareGibbsSampler
from m6a_motif_project import config

def main():
    print("--- 10-Seed Convergence & Stability Experiments ---")
    pos_csv = config.INTERIM_POSITIVE_CSV
    struct_npy = os.path.join("data", "interim", "structure", "positive_unpaired_probs.npy")
    log_dir = os.path.join("results", "logs")
    os.makedirs(log_dir, exist_ok=True)
    
    if not os.path.exists(pos_csv) or not os.path.exists(struct_npy):
        print("Error: Required interim data files are missing. Run prepare_dataset.py first.")
        return
        
    print("Loading data tracks into memory...")
    df = pd.read_csv(pos_csv)
    
    # Subsample 1,000 items with a fixed seed for reproducible execution across samplers
    sample_indices = df.sample(n=1000, random_state=config.RANDOM_SEED).index
    sequences = df.loc[sample_indices, 'sequence'].tolist()
    
    full_struct_matrix = np.load(struct_npy)
    struct_matrix = full_struct_matrix[sample_indices]
    
    steps = 150
    num_chains = 10
    
    results = {
        "baseline": {"consensus": [], "ll_history": []},
        "structure": {"consensus": [], "ll_history": []}
    }
    
    for seed in range(num_chains):
        print(f"\n--- Running Experimental Seed Cluster {seed + 1}/{num_chains} (Seed: {seed}) ---")
        
        # 1. Run Baseline Sequence Chain
        print("  Evaluating Baseline Sampler Chain...")
        base_sampler = BaselineGibbsSampler(sequences, motif_width=config.MOTIF_WIDTH, seed=seed)
        for _ in range(steps):
            base_sampler.step()
        results["baseline"]["consensus"].append(base_sampler.get_consensus())
        results["baseline"]["ll_history"].append(base_sampler.ll_history)
        
        # 2. Run Structure-Aware Thermodynamic Chain
        print("  Evaluating Structure-Aware Sampler Chain...")
        struct_sampler = StructureAwareGibbsSampler(
            sequences, struct_matrix, alpha=config.ALPHA_PRIOR, 
            motif_width=config.MOTIF_WIDTH, seed=seed
        )
        for _ in range(steps):
            struct_sampler.step()
        results["structure"]["consensus"].append(struct_sampler.get_consensus())
        results["structure"]["ll_history"].append(struct_sampler.ll_history)

        print(f"  -> Seed {seed} Baseline Final Consensus: {base_sampler.get_consensus()}")
        print(f"  -> Seed {seed} Structure Final Consensus: {struct_sampler.get_consensus()}")

    # Output logs into results directory tree for downstream evaluation metrics
    log_file = os.path.join(log_dir, "experiment_logs.json")
    with open(log_file, "w") as f:
        json.dump(results, f)
        
    print(f"\n✅ All 10 chains finished cleanly. Trjectory metrics saved to: {log_file}")

if __name__ == "__main__":
    main()