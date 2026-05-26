import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from m6a_motif_project import config
from m6a_motif_project.data.parse_encore import load_and_filter_encore
from m6a_motif_project.data.extract_windows import extract_rna_windows

def main():
    # Track paths directly from the production config layout
    raw_bed = config.RAW_BED_PATH
    ref_fasta = config.RAW_FASTA_PATH
    output_csv = config.INTERIM_POSITIVE_CSV
    
    # Safety Check to guide execution
    if not os.path.exists(raw_bed):
        print(f"Error: Missing raw ENCORE BED file at: {raw_bed}")
        print("Please download it or gunzip it into place.")
        return
        
    if not os.path.exists(ref_fasta):
        print(f"Error: Missing reference genome FASTA at: {ref_fasta}")
        print("Please ensure it is gunzipped and matches the name exactly.")
        return
        
    # Execute modular pipeline functions
    filtered_data = load_and_filter_encore(raw_bed)
    processed_windows = extract_rna_windows(filtered_data, ref_fasta)
    
    # Ensure nested target directories are created safely, then write output
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    processed_windows.to_csv(output_csv, index=False)
    print(f"Pipeline complete! Saved processed positive entries to: {output_csv}")

if __name__ == "__main__":
    main()