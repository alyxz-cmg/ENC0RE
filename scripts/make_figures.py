import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

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

def main():
    print("--- Checklist Step 9: Transforming Trajectories Into Visual Figures ---")
    log_file = os.path.join("results", "logs", "experiment_logs.json")
    output_dir = os.path.join("results", "figures")
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(log_file):
        print(f"Error: Log track data {log_file} not found. Run evaluation or experiments first.")
        return
        
    with open(log_file, "r") as f:
        data = json.load(f)
        
    plot_convergence(data, output_dir)
    plot_stability_heatmap(data, output_dir)
    print(f"✅ All publication visual elements exported successfully to: {output_dir}")

if __name__ == "__main__":
    main()