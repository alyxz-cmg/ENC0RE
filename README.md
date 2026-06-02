# ENC0RE: Structure-Aware Gibbs Sampling for m6A Motif Discovery

ENC0RE is a computational genomics pipeline designed to discover $m^6A$ (RNA methylation) sequence motifs by pairing structural thermodynamic constraints with a Bayesian collapsed Gibbs sampling framework. 

By integrating per-nucleotide RNA unpairing probabilities calculated via `RNAplfold`, this project explicitly demonstrates how structural priors smooth out rugged sequence-only search landscapes, rescuing Markov Chain Monte Carlo (MCMC) samplers from phase-shifted local optima and drastically stabilizing convergence.

---

## 📊 Core Empirical Findings

* **Landscape Smoothing:** The sequence-only baseline frequently fractures across multiple phase-shifted local traps (`GACTG`, `CTGGA`). Introducing the structural prior eliminates background wander, increasing the strict mode recovery of the true biological global optimum (`GGACT`) from **40% to 50%**.
* **Convergence Stability Boost:** Pairwise correlation stability among independent random-seeded matrices jumps from **53.3%** in the Baseline up to **80.0%** in the Structure-Aware sampler.
* **Predictive Ceiling Consistency:** Both models reach an identical maximum held-out test split **AUROC of 0.6645**, demonstrating that structural priors function as an explicit *optimization accelerator* rather than altering the core information content of the discovered sequence model.
* **Confounder Control:** Strict dinucleotide-shuffling verified an absolute **0.00% GC-content delta** between the positive and negative training sets, shielding the downstream evaluator from compositional biases.

---

## 📁 Repository Architecture

```text
.
├── data/
│   ├── interim/               # Extracted window sets and structural arrays
│   └── raw/                   # Raw hg38 BED tracks & GENCODE references
├── src/
│   └── m6a_motif_project/
│       ├── data/              # Parsing, negative generation, and QC verification
│       ├── eval/              # Stability, AUROC, and Gelman-Rubin diagnostics
│       ├── motif/             # Gibbs engines, PWM processing, and priors
│       ├── structure/         # RNAplfold execution wrappers
│       └── viz/               # Publication logos and trajectory plotters
├── scripts/                   # Executable terminal entry points
├── results/
│   ├── figures/               # Convergence, stability, and logo exports
│   └── logs/                  # JSON traces of MCMC histories
├── reports/                   # Written manuscripts and final conclusions
├── environment.yml            # Conda environment definition
└── requirements.txt           # Pip dependency stack
```

---

## 🚀 Getting Started & Quickstart
### 1. Environment Setup

Clone this repository and instantiate the pre-configured environment:
```Bash

git clone [https://github.com/YOUR_USERNAME/ENC0RE.git](https://github.com/YOUR_USERNAME/ENC0RE.git)
cd ENC0RE

# Option A: Using Pip
pip install -r requirements.txt

# Option B: Using Conda
conda env create -f environment.yml
conda activate encore
```

### 2. Complete Execution Pipeline

Run the pipeline sequentially from raw data down to publication visuals using the following execution track:
```Bash

# 1. Download primary GENCODE assembly and miCLIP tracks
PYTHONPATH=src python scripts/download_data.py

# 2. Extract context windows and generate 0.00% delta GC-matched negative sets
PYTHONPATH=src python scripts/prepare_dataset.py

# 3. Run Quality Control to verify confounder shielding
PYTHONPATH=src python m6a_motif_project.data.qc

# 4. Spin up 10 independent MCMC chains for both sampling architectures
PYTHONPATH=src python scripts/run_baseline.py
PYTHONPATH=src python scripts/run_structure_aware.py

# 5. Compute Gelman-Rubin, matrix similarity matrices, and test split AUROCs
PYTHONPATH=src python scripts/evaluate_results.py

# 6. Export publication-ready convergence plots, heatmaps, and logos
PYTHONPATH=src python scripts/make_figures.py
```
---

## 📈 Visual Intermediate Artifacts

All experimental results are logged natively to results/logs/experiment_logs.json. Re-running the figure generation script writes the following primary visual assets back into results/figures/:

    convergence_trajectories.png: Tracks the step-by-step log-likelihood stabilization across all 10 initializations.

    stability_heatmaps.png: Illustrates the structural compaction of final state discrepancies.

    main_comparison_logos.png: Displays information-content scaled biological logos derived directly from the true learned PWM spaces.

---

## References

- Linder et al., 2015. *Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome.* Nature Methods. (canonical DRACH motif)
- Lorenz et al., 2011. *ViennaRNA Package 2.0.* Algorithms for Molecular Biology. (RNAfold / RNAplfold)
- Bernhart et al., 2006. *Local RNA base pairing probabilities in large sequences.* Bioinformatics. (RNAplfold algorithm)
- Lawrence et al., 1993. *Detecting subtle sequence signals: a Gibbs sampling strategy for multiple alignment.* Science. (foundational Gibbs motif sampler)
- Heller et al., 2017. *ssHMM: extracting intuitive sequence-structure motifs from high-throughput RNA-binding protein data.* Nucleic Acids Research. (related structure-aware approach)
- ENCORE / RMBase v3.0 documentation.

---