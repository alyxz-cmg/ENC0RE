import numpy as np
import pandas as pd
from ushuffle import shuffle
from m6a_motif_project import config

def generate_dinucleotide_negatives(positive_df):
    """Generates a 1:1 dinucleotide-shuffled negative control set from positive sequences."""
    print("Generating dinucleotide-shuffled negative control set...")
    
    negative_records = []
    # Seed ushuffle's internal random seed for strict pipeline reproducibility
    np.random.seed(config.RANDOM_SEED)
    
    for idx, row in positive_df.iterrows():
        seq_bytes = row['sequence'].encode('utf-8')
        
        # ushuffle works with string bytes and preserves exact dinucleotide frequency distributions
        shuffled_bytes = shuffle(seq_bytes, 2) 
        shuffled_seq = shuffled_bytes.decode('utf-8')
        
        negative_records.append({
            'chrom': row['chrom'],
            'pos': row['pos'],
            'strand': row['strand'],
            'sequence': shuffled_seq
        })
        
    negative_df = pd.DataFrame(negative_records)
    return negative_df

def verify_base_compositions(pos_df, neg_df):
    """Sanity check: Verifies that mononucleotide content distributions are perfectly matched."""
    print("\n--- Running Pipeline Sanity Check ---")
    
    def get_mono_freqs(df):
        all_seqs = "".join(df['sequence'].tolist())
        total_len = len(all_seqs)
        return {base: all_seqs.count(base) / total_len for base in ['A', 'C', 'G', 'U', 'T']}

    pos_freqs = get_mono_freqs(pos_df)
    neg_freqs = get_mono_freqs(neg_df)
    
    print(f"Positive Set Mononucleotide Frequencies: {pos_freqs}")
    print(f"Negative Set Mononucleotide Frequencies: {neg_freqs}")
    
    # Assert tolerance threshold to ensure the shuffling algorithm didn't mutate characters
    for base in ['A', 'C', 'G', 'U', 'T']:
        diff = abs(pos_freqs[base] - neg_freqs[base])
        assert diff < 1e-5, f"Composition mismatch detected on base {base}! Check shuffle routine."
        
    print("Success: Positive and negative background sequences are perfectly composition-matched.")