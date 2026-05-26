import pandas as pd
from pyfaidx import Fasta
from m6a_motif_project import config

def extract_rna_windows(df, fasta_path):
    """Extracts ±NN nt windows around the m6A modification site, respecting strand orientation."""
    print(f"Opening reference genome: {fasta_path}")
    genome = Fasta(fasta_path)
    
    sequences = []
    w = config.FLANK_WINDOW
    
    for idx, row in df.iterrows():
        chrom = row['chrom']
        mod_pos = int(row['chromStart']) 
        strand = row['strand']
        
        # 0-based indexing boundaries for pyfaidx slicing
        start = mod_pos - w
        end = mod_pos + w + 1
        
        try:
            seq_obj = genome[chrom][start:end]
            seq_str = seq_obj.seq.upper()
            
            # Handle minus strand transcripts to reflect real RNA sequences
            if strand == '-':
                seq_str = seq_obj.reverse.complement.seq.upper()
                
            sequences.append({
                'chrom': chrom,
                'pos': mod_pos,
                'strand': strand,
                'sequence': seq_str
            })
        except KeyError:
            # Skip positions flanking chromosome boundaries
            continue
            
    output_df = pd.DataFrame(sequences)
    print(f"Extracted {len(output_df)} RNA window sequences.")
    return output_df