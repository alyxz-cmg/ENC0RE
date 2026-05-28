import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from m6a_motif_project.motif.gibbs_sampler import BaselineGibbsSampler
from m6a_motif_project import config

def main():
    print("--- Running Baseline Gibbs Sampler on Real ENCORE Data ---")
    pos_csv = config.INTERIM_POSITIVE_CSV
    
    if not os.path.exists(pos_csv):
        print(f"Error: Could not find {pos_csv}. Run prepare_dataset.py first.")
        return
        
    print(f"Loading true RNA sequences from {pos_csv}...")
    df = pd.read_csv(pos_csv)
    
    sample_df = df.sample(n=1000, random_state=config.RANDOM_SEED)
    sequences = sample_df['sequence'].tolist()
    
    steps = 150
    num_chains = 5
    
    for seed in range(num_chains):
        print(f"\n--- Starting Chain {seed + 1}/{num_chains} (Seed: {seed}) ---")
        sampler = BaselineGibbsSampler(sequences, motif_width=config.MOTIF_WIDTH, seed=seed)
        
        for step in range(1, steps + 1):
            sampler.step()
            if step % 50 == 0:
                print(f"  Step {step:03d} | Current Consensus: {sampler.get_consensus()}")
                
        print(f"-> Chain {seed + 1} Final Consensus: {sampler.get_consensus()}")

if __name__ == "__main__":
    main()