import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import logomaker
except ImportError:
    logomaker = None
    print("Warning: 'logomaker' library not found. Sequence logos require it. (pip install logomaker)")

def compute_information_content(pwm):
    """
    Converts a Position Weight Matrix (PWM) to an Information Content (IC) matrix.
    
    Args:
        pwm (np.ndarray): Array of shape (4, L) representing frequencies for A, C, G, T.
        
    Returns:
        pd.DataFrame: A formatted dataframe suitable for logomaker.
    """
    # Clip to avoid log2(0) mathematical errors
    eps = 1e-10
    pwm_safe = np.clip(pwm, eps, 1.0)
    
    # Calculate Shannon entropy per position
    entropy = -np.sum(pwm_safe * np.log2(pwm_safe), axis=0)
    
    # Maximum entropy for 4 bases is log2(4) = 2 bits
    ic = 2.0 - entropy
    
    # Scale frequencies by Information Content
    ic_matrix = pwm * ic[np.newaxis, :]
    
    # Convert to DataFrame with explicit base columns
    df = pd.DataFrame(ic_matrix.T, columns=['A', 'C', 'G', 'T'])
    return df

def plot_sequence_logo(pwm, ax=None, title="Sequence Logo"):
    """
    Draws a standard biological sequence logo from a PWM.
    
    Args:
        pwm (np.ndarray): 4xL probability matrix.
        ax (matplotlib.axes.Axes, optional): Axis to plot on.
        title (str): Title for the subplot.
        
    Returns:
        matplotlib.axes.Axes: The generated axis.
    """
    if logomaker is None:
        raise ImportError("Cannot plot sequence logos without logomaker installed.")
        
    df = compute_information_content(pwm)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
        
    logo = logomaker.Logo(df, ax=ax, color_scheme='classic', vpad=.1, width=.8)
    
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    ax.set_ylabel('Information (bits)', fontweight='bold')
    ax.set_title(title, fontweight='bold', pad=10)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels([str(i+1) for i in range(len(df))])
    
    return ax