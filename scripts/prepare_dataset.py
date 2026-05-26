import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from m6a_motif_project.data.parse_encore import load_and_filter_encore
from m6a_motif_project.data.extract_windows import extract_rna_windows

def main():
    RAW_BED = "data/raw/encore/hg38_m6A_site.bed"
    REF_FASTA = "data/raw/reference/GRCh38.primary_assembly.genome.fa"
    OUTPUT_CSV = "data/interim/positives/positive_sequences.csv"
    
    # Safety Check
    if not os.path.exists(RAW_BED) or not os.path.exists(REF_FASTA):
        print(f"Error: Ensure files exist:\n - {RAW_BED}\n - {REF_FASTA}")
        return
        
    # Execute step-by-step modular pipeline
    filtered_data = load_and_filter_encore(RAW_BED)
    processed_windows = extract_rna_windows(filtered_data, REF_FASTA)
    
    # Ensure directory exists and save
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    processed_windows.to_csv(OUTPUT_CSV, index=False)
    print(f"Pipeline complete! Saved positive entries to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()