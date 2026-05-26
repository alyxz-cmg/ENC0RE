import os
import sys
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from m6a_motif_project import config
from m6a_motif_project.data.parse_encore import load_and_filter_encore
from m6a_motif_project.data.extract_windows import extract_rna_windows
from m6a_motif_project.data.make_negatives import generate_dinucleotide_negatives, verify_base_compositions

def main():
    raw_bed = config.RAW_BED_PATH
    ref_fasta = config.RAW_FASTA_PATH
    pos_output_csv = config.INTERIM_POSITIVE_CSV
    
    # Configure path target for the new negative output matching your tree layout
    neg_output_csv = os.path.join("data", "interim", "negatives", "negative_sequences.csv")
    
    if not os.path.exists(raw_bed) or not os.path.exists(ref_fasta):
        print(f"Error: Raw resource assets missing. Check your download logs.")
        return
        
    # 1. Step 1: Positives Generation
    filtered_data = load_and_filter_encore(raw_bed)
    positive_df = extract_rna_windows(filtered_data, ref_fasta)
    
    os.makedirs(os.path.dirname(pos_output_csv), exist_ok=True)
    positive_df.to_csv(pos_output_csv, index=False)
    print(f"Saved positive entries to: {pos_output_csv}")
    
    # 2. Step 2: Negatives Generation & Quality Verification
    negative_df = generate_dinucleotide_negatives(positive_df)
    verify_base_compositions(positive_df, negative_df)
    
    os.makedirs(os.path.dirname(neg_output_csv), exist_ok=True)
    negative_df.to_csv(neg_output_csv, index=False)
    print(f"Saved matched negative control entries to: {neg_output_csv}")
    print("\nData preparation pipeline phase complete.")

if __name__ == "__main__":
    main()