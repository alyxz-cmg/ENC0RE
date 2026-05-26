# Project: Does RNA Secondary Structure Improve MCMC Motif Recovery? A Case Study on Human m6A Sites from ENCORE
## 1. Overview

This project investigates whether incorporating local RNA secondary structure information enhances the ability of a Bayesian MCMC motif discovery algorithm to recover the canonical m6A motif (DRACH) from human transcript data. The core comparison is between a sequence-only Gibbs sampler and a structure-aware variant that uses unpaired probabilities derived from RNAplfold as a prior. The goal is to produce a reproducible pipeline, quantitative analysis, and insights into the role of RNA structure in motif detection.
## 2. Objectives

    Extract high-confidence human m6A sites from ENCORE.
    Generate matched negative control sequences via dinucleotide shuffling.
    Implement a baseline sequence-only MCMC motif finder (Gibbs sampler).
    Implement a structure-aware MCMC variant incorporating RNA secondary structure priors.
    Compare motif recovery, stability, and convergence between the two methods.
    Document the process, results, and biological interpretation.

## 3. Data Sources
Data	Source	Notes
Human m6A sites	ENCORE	Filter to high-confidence single-nucleotide-resolution entries
Reference genome	GRCh38 via GENCODE v44	Used for sequence extraction around sites
Transcript annotations	GENCODE v44 GTF	For coordinate mapping
Canonical DRACH PWM	Linder et al. 2015 (literature)	Ground truth for evaluation
RNA structure prediction	ViennaRNA's RNAplfold	Local unpaired probabilities

All sources are publicly available and free to use.
## 4. Tools & Libraries
Python environment

    Python 3.10+
    numpy, scipy — numerical computation
    pandas — data manipulation
    biopython or pyfaidx — sequence extraction from FASTA
    ushuffle (or equivalent) — dinucleotide-preserving shuffling for negatives
    logomaker — sequence logo generation
    matplotlib, seaborn — figures

External CLI tools

    ViennaRNA package, specifically RNAplfold for local base-pair probabilities
        Recommended parameters: -W 80 -L 40 -u 8

Recommended hardware

    MacBook Air M2 is sufficient for the full MVP pipeline.
    RTX 2070 Super is optional and only useful for parallelizing chains across CUDA streams as a stretch goal.

## 5. Repository Structure (suggested)

mcmc-m6a-structure/  
├── README.md  
├── HandDown.md                # this file  
├── environment.yml            # conda environment spec  
├── data/  
│   ├── raw/                   # original ENCORE downloads (gitignored)  
│   ├── processed/             # extracted sequences, negatives  
│   └── structure/             # RNAplfold outputs  
├── src/  
│   ├── data_pipeline.py       # ENCORE → sequences → negatives  
│   ├── structure.py           # RNAplfold wrapper  
│   ├── mcmc_baseline.py       # vanilla Gibbs sampler  
│   ├── mcmc_structure.py      # structure-aware variant  
│   ├── evaluation.py          # PWM similarity, AUROC, stability  
│   └── plotting.py            # logos, heatmaps, convergence  
├── notebooks/  
│   ├── 01_data_exploration.ipynb  
│   ├── 02_synthetic_validation.ipynb  
│   ├── 03_real_data_run.ipynb  
│   └── 04_results_figures.ipynb  
├── results/  
│   ├── figures/  
│   └── tables/  
└── report/  
    └── final_report.pdf  

## 6. Pipeline / Workflow
Step 1: Data acquisition

    Download human m6A sites from ENCORE (BED-like format).
    Filter to high-confidence single-nucleotide sites.
    Map each site to GRCh38 coordinates and extract ±25 nt flanking sequence (51 nt windows).

Step 2: Negative set construction

    Generate one dinucleotide-shuffled sequence per positive (1:1 ratio).
    Verify base-composition matching as sanity check.

Step 3: Structure prediction

    Run RNAplfold on each positive sequence with parameters -W 80 -L 40 -u 8.
    Parse output to obtain per-position P(unpaired) values.

Step 4: Sequence-only MCMC (baseline)

    Standard collapsed Gibbs sampler over motif start positions.
    PWM with Dirichlet pseudocounts on each column.
    Uniform prior on motif start position.
    Motif width fixed at 5 nt (DRACH length).
    Validate on synthetic planted-motif data before running on real data.

Step 5: Structure-aware MCMC

    Identical to baseline except the prior on motif start position is proportional to mean P(unpaired) over the candidate window.
    One additional hyperparameter: prior weight α (held fixed at a sensible default for MVP; sweep as stretch goal).

Step 6: Evaluation

    PWM similarity: Pearson correlation between recovered PWM and canonical DRACH PWM.
    Stability: pairwise PWM correlation across ≥10 chains with different random seeds.
    Convergence: Gelman-Rubin diagnostic on chain log-likelihoods.
    Information content: per-position bits.
    Held-out AUROC: PWM-score positives vs negatives on a 20% held-out set.

## 7. Key Design Decisions (and Why)
Decision	Rationale
m6A only (not pseudouridine, m5C, etc.)	DRACH provides clean ground truth; minimizes biological prerequisites for the team
±25 nt window (51 nt total)	Captures local structural context without bloating RNAplfold runtime
Dinucleotide-shuffled negatives only	Computational standard; avoids region-stratification complexity
RNAplfold (not RNAfold)	Local windowed structure is biologically appropriate for short motif scale and runs faster
Position-level structure prior (Option A)	Cleanly slots into Gibbs sampler full conditional with one extra term; column-wise priors (Option B) deferred to stretch
Fixed motif width = 5 (DRACH length)	Avoids motif-width search; can extend to 6–7 as sensitivity check
Multiple random seeds (≥10)	Stability metric is the most likely place where structure prior shows a clean win

## 8. Validation Strategy

Two-tier sanity check before reporting any result:

    Synthetic data: Plant a known motif (e.g., GGACU) into random sequences. Both samplers must recover it. If they don't, there is a bug — fix before proceeding.
    Real-data baseline: The sequence-only sampler must recover DRACH-like output on real m6A data before the structure-aware comparison is meaningful.

Only after both checks pass should comparative claims (structure-aware vs sequence-only) be made.

## 9. Known Issues & Gotchas

    ENCORE format quirks: Coordinate systems (0-based vs 1-based) and strand conventions vary by source. Always confirm by manually checking a few sites against IGV or the UCSC Browser.
    RNAplfold runtime: Scales with sequence length. Run on extracted 51 nt windows rather than full transcripts; do not feed long transcripts.
    MCMC convergence: Multiple short chains beat one long chain for diagnostic purposes. Always log seeds.
    DRACH "looks right" trap: Visual inspection of a motif logo is not sufficient evidence. Always quantify (PWM correlation + AUROC + stability).
    Confounder check: Make sure positives and negatives have similar GC content distributions; if they don't, the structure prior may pick up GC bias rather than true accessibility signal.
    Reverse complement: RNA is single-stranded; do not search the reverse complement strand. This is a common mistake when adapting DNA motif tools.

## 10. Reproducibility Checklist

    Random seeds logged for every MCMC chain run
    Conda environment file (environment.yml) committed
    Raw data download script (or explicit URLs + dates) in repo
    All hyperparameters in a single config file (not hardcoded)
    README includes step-by-step run instructions
    Figures regenerable from saved intermediate results

## 11. Team Roles
Bio-focused member (~30% of total work)

    Data filtering decisions and biological motivation
    DRACH PWM extraction from literature
    Interpretation of recovered motifs
    Most of the report's background and discussion sections
    Bio prereq: ~30 min of reading on m6A, DRACH, and METTL3 substrate preference

CS/ML-focused member (~70% of total work)

    All code (pipeline, samplers, structure integration, evaluation, plotting)
    Hyperparameter and validation experiments
    Repository structure and reproducibility infrastructure
    Methods and results sections of the report

## 12. Timeline Snapshot

    Days 1–2: Data download, sequence extraction, negative generation
    Days 3–5: Baseline MCMC implementation + synthetic validation
    Days 6–7: RNAplfold integration + structure-aware variant
    Days 8–10: Multi-chain runs + evaluation
    Days 11–12: Figures + writeup
    Day 13: Polish + optional Streamlit demo
    Day 14: Buffer / submission

## 13. Success Criteria

MVP (must hit):

    Working baseline that recovers DRACH on real data
    Working structure-aware variant
    One quantitative comparison figure
    Final written report

Strong project (target):

    ≥10 chains per method with stability metrics
    Held-out AUROC for both methods
    Directional answer to the research question with confidence intervals
    Polished, reproducible repo

Stretch:

    α-sweep over structure prior weight
    Streamlit web demo
    Position-specific structure prior (Option B)
    Pseudouridine as secondary modification

## 14. Future Work / Open Questions

    Region-stratified motifs: Does DRACH differ in 5′UTR vs CDS vs 3′UTR? (Was descoped from this project but worth revisiting.)
    Other modifications: Apply the same structure-aware framework to pseudouridine (ψ), m5C, or 2′-O-methylation.
    Position-specific structural priors: Some motif positions may prefer paired vs unpaired contexts (Option B).
    Comparison with existing tools: Benchmark against MEME, RNAcontext, ssHMM.
    Cell-type-specific motifs: Are recovered motifs identical across tissues/cell lines, or context-dependent?
    Writer attribution: Compare METTL3-dependent vs METTL16-dependent sites; do they yield different motifs?

## 15. References

    Linder et al., 2015. Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. Nature Methods.
    Lorenz et al., 2011. ViennaRNA Package 2.0. Algorithms for Molecular Biology.
    Bernhart et al., 2006. Local RNA base pairing probabilities in large sequences. Bioinformatics.
    Lawrence et al., 1993. Detecting subtle sequence signals: a Gibbs sampling strategy for multiple alignment. Science.
    Heller et al., 2017. ssHMM: extracting intuitive sequence-structure motifs from high-throughput RNA-binding protein data. Nucleic Acids Research.
    ENCORE / RMBase v3.0 documentation.
