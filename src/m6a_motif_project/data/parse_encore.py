import pandas as pd

def load_and_filter_encore(bed_path):
    """Loads raw ENCORE m6A data and filters for high-confidence single-nucleotide sites."""
    print(f"Parsing raw ENCORE bed file: {bed_path}")
    
    # BED6 standard columns
    columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    df = pd.read_csv(bed_path, sep='\t', names=columns, comment='#')
    
    # QC: Ensure single-nucleotide resolution (End - Start == 1)
    filtered_df = df[(df['chromEnd'] - df['chromStart']) == 1]
    
    print(f"Filtered down to {len(filtered_df)} high-confidence single-nucleotide sites.")
    return filtered_df