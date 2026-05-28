import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from m6a_motif_project import config
from m6a_motif_project.motif.gibbs_sampler import BaselineGibbsSampler

def generate_synthetic_data(num_seqs=200, seq_len=51, planted_motif="GGACT"):
    """Creates synthetic text windows containing an explicitly embedded motif target."""
    rng = np.random.default_rng(config.RANDOM_SEED)
    alphabet = ['A', 'C', 'G', 'T']
    
    sequences = []
    true_starts = []
    
    for _ in range(num_seqs):
        seq_array = rng.choice(alphabet, size=seq_len)
        start_idx = rng.integers(0, seq_len - len(planted_motif) + 1)
        seq_array[start_idx:start_idx+len(planted_motif)] = list(planted_motif)
        
        sequences.append("".join(seq_array))
        true_starts.append(start_idx)
        
    return sequences, true_starts

def main():
    print("--- Starting Checklist Step 4: Synthetic Validation Baseline ---")
    print("Generating background sequences with a planted 'GGACT' motif...")
    sequences, true_starts = generate_synthetic_data()
    
    sampler = BaselineGibbsSampler(sequences, motif_width=5, seed=config.RANDOM_SEED)
    print(f"Initial Random Consensus: {sampler.get_consensus()}")
    
    steps = 100
    print(f"Running MCMC updates across {steps} steps...")
    for step in range(1, steps + 1):
        sampler.step()
        if step % 20 == 0 or step == steps:
            acc = np.mean(sampler.z == true_starts) * 100
            print(f"Step {step:03d} | Consensus: {sampler.get_consensus()} | Position Accuracy: {acc:.1f}%")
            
    if sampler.get_consensus() == "GGACT":
        print("\n✅ SUCCESS: Baseline Gibbs sampler passed validation check.")
    else:
        print("\n❌ FAILED: Convergence error. Review likelihood weights.")

if __name__ == "__main__":
    main()