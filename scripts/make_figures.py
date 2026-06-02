import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from m6a_motif_project.viz.logos import plot_sequence_logo
from m6a_motif_project import config

def plot_convergence(data, output_dir):
    """Generates comparative log-likelihood convergence plots for all 10 seeds."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    steps = np.arange(len(data["baseline"]["ll_history"][0]))
    
    # Plot Baseline Chains
    for i, history in enumerate(data["baseline"]["ll_history"]):
        axes[0].plot(steps, history, alpha=0.6, label=f"Seed {i}" if i < 3 else "")
    axes[0].set_title("Baseline Gibbs Sampler Convergence", fontsize=12, fontweight='bold')
    axes[0].set_xlabel("MCMC Iteration Step")
    axes[0].set_ylabel("Log-Likelihood")
    axes[0].grid(True, linestyle="--", alpha=0.5)
    
    # Plot Structure-Aware Chains
    for i, history in enumerate(data["structure"]["ll_history"]):
        axes[1].plot(steps, history, alpha=0.6, color='crimson')
    axes[1].set_title("Structure-Aware Gibbs Sampler Convergence", fontsize=12, fontweight='bold')
    axes[1].set_xlabel("MCMC Iteration Step")
    axes[1].grid(True, linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "convergence_trajectories.png"), dpi=300)
    plt.close()
    print("  -> Generated: convergence_trajectories.png")

def plot_stability_heatmap(data, output_dir):
    """Builds a heatmap showing why the structure prior stabilizes cross-run outcomes."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, model in enumerate(["baseline", "structure"]):
        ll_finals = [history[-1] for history in data[model]["ll_history"]]
        # Generate a synthetic distance matrix based on absolute log-likelihood discrepancies
        matrix = np.zeros((10, 10))
        for i in range(10):
            for j in range(10):
                matrix[i, j] = abs(ll_finals[i] - ll_finals[j])
                
        sns.heatmap(matrix, ax=axes[idx], cmap="rocket_r", annot=False, cbar=True)
        axes[idx].set_title(f"{model.upper()} Final State Discrepancy Matrix", fontsize=11, fontweight='bold')
        axes[idx].set_xlabel("Chain Seed Index")
        axes[idx].set_ylabel("Chain Seed Index")
        
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "stability_heatmaps.png"), dpi=300)
    plt.close()
    print("  -> Generated: stability_heatmaps.png")

def plot_main_comparison_logos(data, output_dir, pos_csv):
    """Generates side-by-side sequence logos for the best Baseline and Structure chains."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    pos_df = pd.read_csv(pos_csv)
    alphabet = ['A', 'C', 'G', 'T']
    char_to_idx = {c: i for i, c in enumerate(alphabet)}
    
    for idx, model in enumerate(["baseline", "structure"]):
        ll_finals = [history[-1] for history in data[model]["ll_history"]]
        best_chain_idx = np.argmax(ll_finals)
        best_consensus = data[model]["consensus"][best_chain_idx]

        matched_windows = []
        for seq in pos_df['sequence']:
            seq_clean = seq.replace('U', 'T')
            match_idx = seq_clean.find(best_consensus)
            if match_idx != -1:
                matched_windows.append([char_to_idx.get(c, 0) for c in seq_clean[match_idx:match_idx+5]])
                
        counts = np.zeros((4, 5))
        for win in matched_windows:
            for pos, base in enumerate(win):
                counts[base, pos] += 1
                
        pwm = (counts + 1.0) / (counts + 1.0).sum(axis=0, keepdims=True)
        
        title = f"{model.upper()} Top Motif\n(Chain {best_chain_idx}, LL: {ll_finals[best_chain_idx]:.0f})"
        plot_sequence_logo(pwm, ax=axes[idx], title=title)
        
    plt.tight_layout()
    logo_path = os.path.join(output_dir, "main_comparison_logos.png")
    plt.savefig(logo_path, dpi=300)
    plt.close()
    print(f"  -> Generated: main_comparison_logos.png")

def main():
    print("--- Transforming Trajectories Into Visual Figures ---")
    log_file = os.path.join("results", "logs", "experiment_logs.json")
    pos_csv = config.INTERIM_POSITIVE_CSV
    output_dir = os.path.join("results", "figures")
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(log_file):
        print(f"Error: Log track data {log_file} not found.")
        return
        
    with open(log_file, "r") as f:
        data = json.load(f)
        
    plot_convergence(data, output_dir)
    plot_stability_heatmap(data, output_dir)
    
    plot_main_comparison_logos(data, output_dir, pos_csv)
    
    print(f"✅ All publication visual elements exported successfully to: {output_dir}")

if __name__ == "__main__":
    main()