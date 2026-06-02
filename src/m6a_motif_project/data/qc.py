import os
import sys
import pandas as pd
import numpy as np

# Add the parent directory to the path so config can be imported smoothly if run directly
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))
from m6a_motif_project import config

def calculate_gc_content(sequence):
    """Calculates the percentage of G and C bases in a sequence string."""
    seq_upper = sequence.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    if len(seq_upper) == 0:
        return 0.0
    return (gc_count / len(seq_upper)) * 100

def verify_confounders():
    """Checks for GC-content matching between positive and negative cohorts."""
    print("--- Running Confounder Verification Checks ---")
    
    pos_csv = config.INTERIM_POSITIVE_CSV
    # Derived relative path for the negatives based on project structure
    neg_csv = os.path.join("data", "interim", "negatives", "negative_sequences.csv")
    
    if not os.path.exists(pos_csv) or not os.path.exists(neg_csv):
        print(f"Error: Datasets missing. Ensure {pos_csv} and {neg_csv} exist.")
        return False
        
    pos_df = pd.read_csv(pos_csv)
    neg_df = pd.read_csv(neg_csv)
    
    pos_gc = pos_df['sequence'].apply(calculate_gc_content).to_numpy()
    neg_gc = neg_df['sequence'].apply(calculate_gc_content).to_numpy()
    
    mean_pos_gc = np.mean(pos_gc)
    mean_neg_gc = np.mean(neg_gc)
    
    print(f"Positive Set Mean GC: {mean_pos_gc:.2f}%")
    print(f"Negative Set Mean GC: {mean_neg_gc:.2f}%")
    
    gc_delta = abs(mean_pos_gc - mean_neg_gc)
    print(f"Absolute GC Difference: {gc_delta:.2f}%")
    
    if gc_delta < 1.0:
        print("✅ SUCCESS: GC-content is perfectly controlled between groups (delta < 1%).")
        return True
    else:
        print("⚠️ WARNING: GC bias detected between cohorts. Verify your shuffle script step.")
        return False

if __name__ == "__main__":
    verify_confounders()