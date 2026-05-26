import os
import numpy as np

# Reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Data Pipeline Settings
FLANK_WINDOW = 25  # ±25 nt yields a 51 nt total window

# Project-Specific Source URLs for Documentation/Scripts
RMBASE_M6A_URL = "https://raw.githubusercontent.com/mevers/RNAModR/master/inst/extdata/miCLIP_m6A_Linder2015_hg38.bed"
GENCODE_FASTA_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
GENCODE_GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"

# Concrete Tree File Paths
RAW_BED_PATH = os.path.join("data", "raw", "encore", "hg38_m6A_site.bed")
RAW_FASTA_PATH = os.path.join("data", "raw", "reference", "GRCh38.primary_assembly.genome.fa")
RAW_GTF_PATH = os.path.join("data", "raw", "reference", "gencode.v44.primary_assembly.annotation.gtf")

INTERIM_POSITIVE_CSV = os.path.join("data", "interim", "positives", "positive_sequences.csv")

# MCMC Modeling Settings
MOTIF_WIDTH = 5   # Canonical DRACH length
ALPHA_PRIOR = 1.0  # Structure prior weight hyperparameter
NUM_CHAINS = 10
MCMC_STEPS = 5000