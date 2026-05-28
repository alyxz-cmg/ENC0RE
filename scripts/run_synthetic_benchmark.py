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
    print("--- Starting Synthetic Validation Baseline ---")
    print("Generating background sequences with a planted 'GGACT' motif...\n")
    sequences, true_starts = generate_synthetic_data()
    
    steps = 100
    num_chains = 5
    
    success = False
    
    for seed in range(num_chains):
        print(f"--- Running Chain {seed + 1}/{num_chains} (Seed: {seed}) ---")
        sampler = BaselineGibbsSampler(sequences, motif_width=5, seed=seed)
        
        for step in range(1, steps + 1):
            sampler.step()
            
        consensus = sampler.get_consensus()
        acc = np.mean(sampler.z == true_starts) * 100
        
        # Check for phase shifts to provide better debugging feedback
        shift_minus_1 = np.mean(sampler.z == (np.array(true_starts) - 1)) * 100
        shift_plus_1 = np.mean(sampler.z == (np.array(true_starts) + 1)) * 100
        
        print(f"Final Consensus: {consensus}")
        print(f"Exact Match Accuracy: {acc:.1f}%")
        if shift_minus_1 > 50:
            print(f"⚠️ Phase Shift Detected: Chain got trapped -1 position to the left ({shift_minus_1:.1f}% shifted)")
        elif shift_plus_1 > 50:
            print(f"⚠️ Phase Shift Detected: Chain got trapped +1 position to the right ({shift_plus_1:.1f}% shifted)")
            
        print("-" * 40)
        
        if consensus == "GGACT" and acc > 90:
            success = True

    if success:
        print("\n✅ SUCCESS: At least one chain cleanly recovered the global optimum planted motif.")
    else:
        print("\n❌ FAILED: All chains diverged or got trapped in local optima.")

if __name__ == "__main__":
    main()