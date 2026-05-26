import numpy as np

# Reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Data Pipeline Settings
FLANK_WINDOW = 25  # ±25 nt yields a 51 nt total window

# Source URLs for Documentation/Scripts
RMBASE_M6A_URL = "http://bioinformaticsscience.cn/rmbase/download/hg38/hg38_m6A_site.bed.gz"
GENCODE_FASTA_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
GENCODE_GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"

# MCMC Modeling Settings
MOTIF_WIDTH = 5   # Canonical DRACH length
ALPHA_PRIOR = 1.0  # Structure prior weight hyperparameter
NUM_CHAINS = 10
MCMC_STEPS = 5000