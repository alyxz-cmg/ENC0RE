**Project Title:** Does RNA Secondary Structure Improve MCMC Motif Recovery? A Case Study on Human m6A Sites from ENCORE

**Document Version:** 1.0
**Date:** May 21, 2026
**Duration:** 2 weeks
**Compute Resources:** MacBook Air (M2), Desktop (RTX 2070 Super, 8GB VRAM)

---

## 1. Executive Summary

This project investigates whether incorporating local RNA secondary structure information improves the recovery of sequence motifs around N6-methyladenosine (m6A) modification sites in human transcripts. We will implement two variants of an MCMC-based motif discovery algorithm — a standard sequence-only Gibbs sampler and a structure-aware variant that incorporates per-position unpaired probabilities as a Bayesian prior — and compare their ability to recover the canonical DRACH consensus motif from m6A sites curated by the ENCORE epitranscriptome resource. The output is a methodological case study with quantitative evidence on whether structural priors meaningfully aid de novo motif discovery in the human epitranscriptome.

---

## 2. Background & Motivation

RNA modifications, particularly m6A, are central to post-transcriptional gene regulation, influencing mRNA stability, translation efficiency, and splicing. m6A is deposited by the METTL3/METTL14 writer complex and shows a well-characterized consensus motif known as DRACH (D = A/G/U; R = A/G; A = the modified adenosine; C; H = A/C/U).

Despite the strength of the DRACH consensus, two open questions motivate this project:

1. **Sequence is not the whole story.** The METTL3/14 complex is biochemically known to prefer single-stranded (accessible) substrate, suggesting that local RNA folding biases where m6A is actually deposited beyond what sequence alone predicts.
2. **De novo motif discovery rarely uses structure.** Most widely used motif discovery tools (MEME, Gibbs sampler variants) ignore RNA secondary structure entirely. The few structure-aware tools (RNAcontext, ssHMM) are specialized and not commonly benchmarked against modern Bayesian samplers in the modification-site setting.

This creates a clean methodological question with a known ground-truth motif (DRACH) for evaluation, making it an ideal case study.

---

## 3. Research Question & Hypothesis

**Primary research question:**
Does incorporating local RNA secondary structure (per-position unpaired probabilities) as a prior in MCMC motif discovery improve the recovery, stability, and confidence of the m6A sequence motif compared to a sequence-only sampler?

**Hypothesis:**
A structure-aware MCMC sampler will recover DRACH-like motifs with (a) higher similarity to the canonical PWM, (b) greater stability across random seeds, and (c) faster convergence than a sequence-only sampler — because m6A sites are biologically enriched in single-stranded contexts.

**Null result is also informative:** If structure priors do not improve motif recovery, this suggests that m6A motif identity is sufficiently encoded in primary sequence and that structural context affects binding strength rather than motif-level sequence preference.

---

## 4. Objectives

1. Build a reproducible pipeline for extracting m6A site sequences and matched negative controls from ENCORE.
2. Implement a vanilla Gibbs/MCMC motif discovery sampler from scratch.
3. Implement a structure-aware variant that uses RNAplfold-derived unpaired probabilities as an informative prior on motif start positions.
4. Quantitatively compare both samplers on motif recovery (DRACH similarity), stability (across random seeds), and convergence diagnostics.
5. Produce a written report and reproducible codebase suitable as a course final project.

---

## 5. Scope

### In Scope
- One RNA modification: **m6A**
- One organism: **human**
- One data source: **ENCORE**, filtered to high-confidence single-nucleotide-resolution sites
- One window size: ±25 nt centered on each modification site (51 nt total)
- One negative control type: dinucleotide-shuffled positives
- Two methods: sequence-only MCMC vs structure-aware MCMC
- One structure prior formulation: motif-start position weighted by mean P(unpaired) over the window (Option A from project planning)

### Out of Scope
- Multiple modification types (e.g., pseudouridine, m5C)
- Region-stratified analysis (5′UTR vs CDS vs 3′UTR)
- Cross-cell-line or cross-tissue comparisons
- Position-specific structural priors (Option B / column-wise)
- Comparison against existing tools (MEME, RNAcontext, ssHMM) — mentioned as stretch goal only
- Genome-wide scanning or de novo prediction of new m6A sites

---

## 6. Success Criteria

### Minimum Viable Project (must hit all)
- Working sequence-only MCMC sampler that recovers a DRACH-like motif from real m6A data
- Working structure-aware variant integrated with RNAplfold output
- At least one quantitative comparison figure between the two methods
- Written report covering motivation, methods, results, and interpretation

### Strong Project (target)
- Multi-chain runs (≥10 random seeds per method) with stability metrics
- Held-out classification AUROC for both methods
- Clear directional answer to the research question, with confidence intervals
- Reproducible repository with README and environment file

### Stretch (if time permits)
- Streamlit or simple web demo
- Comparison against MEME or another standard tool
- Pseudouridine as a secondary modification for cross-modification comparison

---

## 7. Methodology

### 7.1 Data Pipeline
1. Download high-confidence human m6A sites from ENCORE
2. Filter to single-nucleotide-resolution sites with confidence above a defined threshold
3. Map sites to the GRCh38 / GENCODE v44 reference; extract ±25 nt flanking sequence per site
4. Generate dinucleotide-shuffled negatives (1:1 ratio with positives) using a standard library (e.g., `ushuffle`)

### 7.2 Structure Prediction
- Run **RNAplfold** on each positive sequence with parameters `-W 80 -L 40 -u 8`
- Output per-position unpaired probabilities for use in the structure-aware prior

### 7.3 Sequence-Only MCMC (Baseline)
- Standard collapsed Gibbs sampler over motif start positions
- PWM with Dirichlet pseudocounts on each column
- Uniform prior on motif start position
- Motif width fixed at 5 nt (DRACH length); optionally explored at 6–7 nt

### 7.4 Structure-Aware MCMC
- Identical to baseline, except prior on motif start position is proportional to mean P(unpaired) across the candidate window
- One additional hyperparameter (prior weight α) to control structure influence; held fixed at a reasonable default for MVP, swept as a stretch goal

### 7.5 Validation
- **Synthetic test:** plant a known motif into random sequences, verify both samplers recover it
- **Sanity test:** confirm baseline recovers DRACH on real m6A data before claiming any structure-aware improvement

### 7.6 Evaluation Metrics
- **Motif similarity:** Pearson correlation between recovered PWM columns and the canonical DRACH PWM; optionally Tomtom-style E-value
- **Stability:** pairwise PWM correlation across 10+ MCMC chains with different random seeds
- **Convergence:** Gelman-Rubin diagnostic on chain log-likelihoods
- **Information content:** per-position bits of the recovered motif
- **Held-out AUROC:** PWM-score positives vs negatives on a held-out 20% test set

---

## 8. Data Requirements

| Data | Source | Approximate Size |
|---|---|---|
| Human m6A sites (high-confidence) | ENCORE | ~50k–500k sites pre-filter |
| Reference genome (GRCh38) | GENCODE / Ensembl | ~3 GB |
| Transcript annotations | GENCODE v44 GTF | ~1 GB |
| Canonical DRACH PWM (for evaluation) | Literature (Linder et al. 2015 or equivalent) | Negligible |

All data is publicly available and free.

---

## 9. Technical Requirements

### Languages & Core Libraries
- **Python 3.10+**
- `numpy`, `scipy` — numerical and statistical computation
- `pandas` — data manipulation
- `biopython` or `pyfaidx` — sequence extraction
- `ushuffle` (or equivalent) — dinucleotide shuffling
- `logomaker` — sequence logo plotting
- `matplotlib` / `seaborn` — figures

### External Tools
- **ViennaRNA / RNAplfold** — RNA structure prediction (CLI tool, called via subprocess)

### Hardware
- Primary development on **MacBook Air M2** (sufficient for MCMC and RNAplfold runs)
- **RTX 2070 Super** used optionally for parallelizing chains across CUDA streams (stretch goal; not required for MVP)

### Reproducibility
- All code in a single Git repository
- `environment.yml` (conda) or `requirements.txt`
- Random seeds logged for every MCMC chain
- README with data download steps and run instructions

---

## 10. Deliverables

1. **Codebase** (Git repository)
   - Data preprocessing scripts (ENCORE site extraction, sequence retrieval, negative generation)
   - Vanilla MCMC motif sampler (Python module)
   - Structure-aware MCMC variant
   - RNAplfold integration wrapper
   - Evaluation and plotting scripts

2. **Figures**
   - Sequence logos: sequence-only vs structure-aware recovered motifs, side-by-side with canonical DRACH
   - Stability heatmap: pairwise PWM correlation across chains for each method
   - Convergence plot: log-likelihood traces and Gelman-Rubin diagnostic
   - AUROC bar chart for held-out classification
   - Structure-context plot: average unpaired probability profile around motif occurrences

3. **Final report** (written document)
   - Background and motivation
   - Methods (data, sampler, structure prior, evaluation)
   - Results with figures
   - Discussion: when does structure help, when doesn't it
   - Limitations and future work

4. **Optional (extra credit)**
   - Streamlit web demo for interactive motif discovery on user-supplied sequences
   - Tutorial notebook walking through the full pipeline

---

## 11. Timeline & Milestones

| Day | Bio-focused member | CS/ML-focused member | Milestone |
|---|---|---|---|
| 1 | Download ENCORE m6A sites; basic filtering | Set up repo, environment, RNAplfold install | Data acquired |
| 2 | Sequence extraction QC; document filtering choices | Sequence extraction script; dinucleotide shuffling | Positive + negative sets ready |
| 3 | Begin background section of report | Implement vanilla Gibbs sampler | Baseline sampler skeleton |
| 4 | Background draft; pull DRACH PWM from literature | Validate baseline on synthetic planted motifs | Sampler validated on synthetic data |
| 5 | Light review of m6A writer biology | Run baseline on real m6A; confirm DRACH recovery | Baseline working on real data |
| 6 | Help interpret early baseline results | Run RNAplfold across all positives | Structure features computed |
| 7 | Methods section draft | Implement structure-aware prior | Structure-aware sampler skeleton |
| 8 | — | Validate structure-aware version on synthetic data | Both samplers operational |
| 9 | — | Multi-chain runs (10+ seeds per method) | Stability data collected |
| 10 | — | AUROC evaluation; convergence diagnostics | Quantitative comparison done |
| 11 | Results section draft; help interpret motifs | Generate all final figures | Figures finalized |
| 12 | Discussion section draft | Code cleanup, README, environment file | Repo polished |
| 13 | Polish writeup | Optional: Streamlit demo | Report near-final |
| 14 | Final pass on writeup | Buffer / final figure tweaks | **Submission** |

---

## 12. Team & Responsibilities

### Bio-focused member (~30% of total work)
- Data filtering decisions (which ENCORE confidence threshold; which transcript types to include)
- Biological motivation and background section of the report
- Interpretation of recovered motifs (does it look like DRACH? are deviations sensible?)
- Sanity-checking that the analysis remains biologically defensible
- **Required bio knowledge upfront:** ~3 sentences on m6A; the DRACH consensus; the fact that METTL3 prefers single-stranded substrate. Estimated ~30 min of reading.

### CS/ML-focused member (~70% of total work)
- All code: pipeline, samplers, structure integration, evaluation
- Hyperparameter choices and validation experiments
- Figure generation
- Repository structure and reproducibility

---

## 13. Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Baseline MCMC fails to recover DRACH | Medium | High | Validate first on synthetic planted-motif data; if real data fails, debug with smaller cleaner subset before adding structure prior |
| ENCORE data format harder to parse than expected | Medium | Medium | Have m6A-Atlas / RMBase as backup sources; sites are largely overlapping |
| Structure prior makes no difference | Medium | Low (still a finding) | Frame as a legitimate negative result; report exactly when and why structure does/doesn't help |
| RNAplfold runtime explodes on long transcripts | Low | Medium | Use local windows (`-W 80`); run on extracted ±25 nt sequences rather than whole transcripts |
| MCMC convergence is slow | Medium | Medium | Use multiple short chains rather than one long chain; subsample sequences if needed |
| Team member time conflicts | Medium | High | Front-load data pipeline so CS member is unblocked; bio member's deliverables are mostly writing |

---

## 14. Stretch Goals (if MVP completed early)

- Sweep over the structure-prior weight α to characterize how much weight is optimal
- Add a position-specific structural prior (Option B from planning)
- Compare against MEME or another standard motif tool
- Extend to pseudouridine (ψ) sites as a contrasting case where no canonical motif is known
- Streamlit web demo: user pastes RNA sequences, gets recovered motif and structure context

---

## 15. References

- Linder et al., 2015. *Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome.* Nature Methods. (canonical DRACH motif)
- Lorenz et al., 2011. *ViennaRNA Package 2.0.* Algorithms for Molecular Biology. (RNAfold / RNAplfold)
- Bernhart et al., 2006. *Local RNA base pairing probabilities in large sequences.* Bioinformatics. (RNAplfold algorithm)
- Lawrence et al., 1993. *Detecting subtle sequence signals: a Gibbs sampling strategy for multiple alignment.* Science. (foundational Gibbs motif sampler)
- Heller et al., 2017. *ssHMM: extracting intuitive sequence-structure motifs from high-throughput RNA-binding protein data.* Nucleic Acids Research. (related structure-aware approach)
- ENCORE / RMBase v3.0 documentation.

---