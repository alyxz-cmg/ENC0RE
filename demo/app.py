import sys
import os
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Run with: PYTHONPATH=src streamlit run demo/app.py

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from m6a_motif_project.motif.gibbs_sampler import BaselineGibbsSampler
from m6a_motif_project.viz.logos import plot_sequence_logo

st.set_page_config(page_title="ENC0RE Motif Discoverer", page_icon="🧬", layout="wide")

st.title("🧬 ENC0RE: $m^6A$ Motif Discovery Dashboard")
st.markdown("""
Welcome to the ENC0RE interactive demo! This tool uses a **Bayesian Gibbs Sampler** to discover hidden RNA methylation motifs from unaligned sequence data.
""")

st.sidebar.header("⚙️ Model Hyperparameters")
iterations = st.sidebar.slider("MCMC Iteration Steps", min_value=10, max_value=500, value=100, step=10)
motif_width = st.sidebar.number_input("Motif Width", min_value=3, max_value=10, value=5)
random_seed = st.sidebar.number_input("Random Seed", min_value=0, max_value=9999, value=42)

# --- Sequence Input Section ---
st.subheader("1. Sequence Input Space")
st.markdown("Upload a `.fasta` file, or proceed with the default provided demo sequences.")

fasta_file = st.file_uploader("Upload FASTA sequences", type=["fa", "fasta"])

@st.cache_data
def read_fasta(file_content_bytes):
    """Parses a FASTA file into a list of pure sequences."""
    lines = file_content_bytes.decode("utf-8").splitlines()
    seqs = []
    for line in lines:
        if not line.startswith(">") and line.strip():
            seqs.append(line.strip().upper())
    return seqs

sequences = []
if fasta_file is not None:
    sequences = read_fasta(fasta_file.read())
    st.success(f"Successfully loaded {len(sequences)} sequences from upload.")
else:
    demo_path = os.path.join(os.path.dirname(__file__), "demo_sequences.fa")
    if os.path.exists(demo_path):
        with open(demo_path, "rb") as f:
            sequences = read_fasta(f.read())
        st.info(f"Loaded {len(sequences)} sequences from default `demo_sequences.fa`.")
    else:
        st.warning("No sequences loaded. Please upload a FASTA file.")

# --- Execution Section ---
st.markdown("---")
st.subheader("2. Run Sequence Motif Discovery")

if st.button("🚀 Execute Gibbs Sampler", type="primary"):
    if not sequences:
        st.error("Cannot run model: No sequence data provided.")
    else:
        with st.spinner("Initializing chains and running MCMC sampling..."):
            sampler = BaselineGibbsSampler(sequences, motif_width=int(motif_width), seed=int(random_seed))
            
            progress_bar = st.progress(0)
            for i in range(iterations):
                sampler.step()
                if i % max(1, iterations // 100) == 0:
                    progress_bar.progress((i + 1) / iterations)
            progress_bar.progress(1.0)
            
            consensus = sampler.get_consensus()
            pwm = sampler.pwm
            
        st.success("Sampling Complete!")
        
        # --- Results Columns ---
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown(f"### Discovered Consensus: **`{consensus}`**")
            fig, ax = plt.subplots(figsize=(6, 4))
            plot_sequence_logo(pwm, ax=ax, title="Learned Biological Logo")
            st.pyplot(fig)
            
        with col2:
            st.markdown("### Position Weight Matrix (PWM)")
            # Format PWM for clean display
            df_pwm = pd.DataFrame(pwm, index=['A', 'C', 'G', 'T'], columns=[f"Pos {i+1}" for i in range(motif_width)])
            st.dataframe(df_pwm.style.background_gradient(cmap='Blues', axis=None))

# --- Project Showcase Section ---
st.markdown("---")
st.header("📊 Final Project Results: The Power of Structural Priors")
st.markdown("""
While the live demo above uses the sequence-only baseline model for speed, our full research pipeline proved that integrating 
**thermodynamic unpairing probabilities (`RNAplfold`)** acts as a powerful optimization prior. 

* **Stabilizes Convergence:** Structural priors rescue chains from phase-shifted background traps.
* **Consolidates Top Modes:** Eliminates heavy sequence variance across random initialization seeds.
""")

col_fig1, col_fig2 = st.columns(2)
try:
    with col_fig1:
        st.image(os.path.join("results", "figures", "main_comparison_logos.png"), 
                 caption="True Biological MCMC PWMs")
    with col_fig2:
        st.image(os.path.join("results", "figures", "convergence_trajectories.png"), 
                 caption="Log-Likelihood Stabilization via Structure")
except Exception:
    st.info("Run `scripts/make_figures.py` to generate the comparison visual artifacts for the dashboard.")