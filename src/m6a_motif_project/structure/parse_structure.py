import os
import numpy as np

def parse_lunp_file(file_path, seq_len=51):
    """Parses a single ViennaRNA _lunp file to extract l=1 unpaired probabilities."""
    p_unpaired = np.zeros(seq_len)
    
    if not os.path.exists(file_path):
        return p_unpaired
        
    with open(file_path, 'r') as f:
        for line in f:
            line_stripped = line.strip()
            # Skip empty lines, comment lines, or structural metadata headers
            if not line_stripped or line_stripped.startswith('#'):
                continue
                
            parts = line_stripped.split()
            if len(parts) >= 2:
                if not parts[0].isdigit():
                    continue
                    
                # 1-based index from file -> 0-based index for numpy array
                pos_idx = int(parts[0]) - 1 
                
                # Column index 1 contains the l=1 single-nucleotide unpairing probability
                try:
                    prob = float(parts[1])
                    if 0 <= pos_idx < seq_len:
                        p_unpaired[pos_idx] = prob
                except ValueError:
                    # Guard against non-numeric data in the probability columns
                    continue
                    
    return p_unpaired