import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from m6a_motif_project.motif.gibbs_sampler import StructureAwareGibbsSampler
from m6a_motif_project import config

def main():
    print("--- Running Structure-Aware Gibbs Sampler ---")
    pos_csv = config.INTERIM_POSITIVE_CSV
    struct_npy = os.path.join("data", "interim", "structure", "positive_unpaired_probs.npy")
    
    if not os.path.exists(pos_csv) or not os.path.exists(struct_npy):
        print("Error: Missing data files. Run prepare_dataset.py first.")
        return
        
    print("Loading RNA sequences and Structural Probabilities...")
    df = pd.read_csv(pos_csv)
    
    # Subsample to keep test runtime manageable, exactly as we did for the baseline
    sample_indices = df.sample(n=1000, random_state=config.RANDOM_SEED).index
    sequences = df.loc[sample_indices, 'sequence'].tolist()
    
    # Load and subset the corresponding structure probabilities
    full_struct_matrix = np.load(struct_npy)
    struct_matrix = full_struct_matrix[sample_indices]
    
    steps = 150
    num_chains = 5
    
    for seed in range(num_chains):
        print(f"\n--- Starting Structure-Aware Chain {seed + 1}/{num_chains} (Seed: {seed}) ---")
        sampler = StructureAwareGibbsSampler(
            sequences=sequences, 
            unpaired_probs=struct_matrix,
            alpha=config.ALPHA_PRIOR,
            motif_width=config.MOTIF_WIDTH, 
            seed=seed
        )
        
        for step in range(1, steps + 1):
            sampler.step()
            if step % 50 == 0:
                print(f"  Step {step:03d} | Current Consensus: {sampler.get_consensus()}")
                
        print(f"-> Chain {seed + 1} Final Consensus: {sampler.get_consensus()}")

if __name__ == "__main__":
    main()